#' Δ prevalence vs pooled prevalence (per peptide)
#'
#' @description
#' Build a static ggplot showing the per-peptide shift in prevalence
#' (Δ = group2 − group1) as a function of pooled prevalence ((group1 + group2)/2).
#' Intended for inputs like `make_prev_tbl()` (long format with columns
#' `peptide_id`, `Group`, `Prevalence`).
#'
#' @param prev_tbl Long table with at least three columns: `peptide_id`,
#'   `group_col` (defaults to `"Group"`), and `value_col` (defaults to `"Prevalence"`).
#' @param group_col Name of the column with the group labels. Default `"Group"`.
#' @param value_col Name of the column with prevalence values. Default `"Prevalence"`.
#' @param peptide_col Name of the peptide id column. Default `"peptide_id"`.
#' @param group_order Optional character vector of length 2 giving the order
#'   (c(group1, group2)). If `NULL`, the function uses the unique values found
#'   in `group_col` (sorted).
#' @param jitter_width,jitter_height Jitter amounts for the points. Defaults 0.005.
#' @param point_alpha Point transparency. Default 0.25.
#' @param point_size Point size. Default 0.6.
#' @param smooth Add a GAM smooth curve (`mgcv`). Default `TRUE`.
#' @param smooth_k Basis dimension `k` for the smooth. Default 5.
#' @param arrow_color Color for the directional arrows and labels. Default `"red"`.
#' @param arrow_length_mm Arrow length in mm. Default 4.
#' @param arrow_pos_x Where to place arrows horizontally as a fraction of the
#'   max pooled prevalence in the data. Default 0.97 (near the right edge).
#' @param title,subtitle,xlab,ylab Optional plot labels.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' # prev_tbl <- make_prev_tbl(x_kid1, rank="is_flagellum", g1="T2", g2="T8", feature_value=TRUE)
#' # p <- deltaplot_prevalence(prev_tbl, group_order = c("B","M12"))
#' # p
#'
#' @export
deltaplot_prevalence <- function(
  prev_tbl,
  group_col = "Group",
  value_col = "Prevalence",
  peptide_col = "peptide_id",
  group_order = NULL,
  jitter_width = 0.005,
  jitter_height = 0.005,
  point_alpha = 0.25,
  point_size = 0.6,
  smooth = TRUE,
  smooth_k = 5,
  arrow_color = "red",
  arrow_length_mm = 4,
  arrow_pos_x = 0.97,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL
) {
  # --- basic checks -----------------------------------------------------------
  need <- c(peptide_col, group_col, value_col)
  miss <- setdiff(need, names(prev_tbl))
  if (length(miss)) {
    stop(
      "deltaplot_prevalence: missing required columns: ",
      paste(miss, collapse = ", ")
    )
  }

  # --- wide table with exactly two groups ------------------------------------
  wide <- prev_tbl |>
    dplyr::select(dplyr::all_of(c(peptide_col, group_col, value_col))) |>
    tidyr::pivot_wider(
      names_from = !!rlang::sym(group_col),
      values_from = !!rlang::sym(value_col)
    )

  grp_names <- setdiff(names(wide), peptide_col)

  if (length(grp_names) != 2L && is.null(group_order)) {
    stop(
      "deltaplot_prevalence: expected exactly two groups in `", group_col,
      "`, found: ", paste(grp_names, collapse = ", ")
    )
  }

  if (!is.null(group_order)) {
    if (length(group_order) != 2) stop("group_order must have length 2.")
    if (!all(group_order %in% grp_names)) {
      stop(
        "group_order values must be present in ", group_col,
        ". Found groups: ", paste(grp_names, collapse = ", ")
      )
    }
    g1 <- group_order[1]
    g2 <- group_order[2]
  } else {
    # stable (alphabetical) order
    g1 <- sort(grp_names)[1]
    g2 <- sort(grp_names)[2]
  }

  # --- pooled & delta ---------------------------------------------------------
  w <- wide |>
    dplyr::mutate(
      pooled = (.data[[g1]] + .data[[g2]]) / 2,
      delta = .data[[g2]] - .data[[g1]]
    )

  w_clean <- w |>
    dplyr::filter(is.finite(pooled), is.finite(delta)) |>
    dplyr::mutate(pooled_clip = pmin(pmax(as.numeric(pooled), 1e-6), 1 - 1e-6))

  if (!nrow(w_clean)) {
    stop("deltaplot_prevalence: no finite rows to plot after computing pooled & delta.")
  }

  # Arrow position near the right edge of the data
  arrow_x <- max(w_clean$pooled_clip, na.rm = TRUE) * arrow_pos_x

  # --- build plot -------------------------------------------------------------
  p <- ggplot2::ggplot(w_clean, ggplot2::aes(pooled_clip, delta)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_jitter(
      alpha = point_alpha,
      height = jitter_height,
      width = jitter_width,
      size = point_size
    )

  if (isTRUE(smooth)) {
    # mgcv::gam; ok to depend on mgcv indirectly via ggplot2
    p <- p + ggplot2::geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = smooth_k),
      se = FALSE
    )
  }

  # directional arrows & labels
  p <- p +
    ggplot2::annotate(
      "segment",
      x = arrow_x, xend = arrow_x, y = 0, yend = 0.06,
      colour = arrow_color, linewidth = 0.6,
      arrow = grid::arrow(length = grid::unit(arrow_length_mm, "mm"))
    ) +
    ggplot2::annotate(
      "text",
      x = arrow_x, y = 0.065,
      label = paste0("More in ", g2),
      colour = arrow_color, fontface = "bold", vjust = 0
    ) +
    ggplot2::annotate(
      "segment",
      x = arrow_x, xend = arrow_x, y = 0, yend = -0.06,
      colour = arrow_color, linewidth = 0.6,
      arrow = grid::arrow(length = grid::unit(arrow_length_mm, "mm"))
    ) +
    ggplot2::annotate(
      "text",
      x = arrow_x, y = -0.065,
      label = paste0("More in ", g1),
      colour = arrow_color, fontface = "bold", vjust = 1
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    ggplot2::labs(
      title = title %||% paste0("Per-peptide shift vs pooled prevalence (", g2, " - ", g1, ")"),
      subtitle = subtitle,
      x = xlab %||% paste0("Pooled prevalence (", g1, " & ", g2, ")"),
      y = ylab %||% paste0("\u0394 prevalence (", g2, " - ", g1, ")")
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(10, 20, 10, 10),
      text = ggplot2::element_text(size = 20)
    )

  p
}

#' Build per-peptide prevalence for a chosen rank/feature and two groups (standalone)
#'
#' @description
#' Construct a long per-peptide prevalence table for a specific `rank` value
#' (e.g., a species or a boolean annotation) comparing two groups in
#' `x$data_long`. Works directly with `phip_data` that carries a
#' `peptide_library`, without calling any other helper.
#'
#' @param x A `phip_data` object with `data_long` and `peptide_library`.
#' @param rank Peptide-library column to filter on (e.g. `"species"` or `"is_flagellum"`).
#' @param feature_value Value(s) to keep in `rank` (e.g. `"Acadevirus PM93"`,
#'   `TRUE`, or a character vector). Character/logical matching is case-insensitive.
#' @param group_col Column in `x$data_long` defining groups (default `"timepoint_factor"`).
#' @param groups Character(2) with the two groups to compare, e.g. `c("T2","T8")`.
#' @param labels Character(2) display labels for those groups (default `c("B","M12")`).
#' @param filter Optional named list of simple equality filters applied to `x$data_long`
#'   before computing prevalences, e.g. `list(big_group = "kid_serum")` or
#'   `list(site = c("gut","oral"))`. (Uses `%in%`.)
#' @param eps Laplace numerator epsilon (default `0.5`).
#' @param den_mult Denominator epsilon multiplier (default `2`).
#'
#' @return A tibble with columns: `peptide_id`, `Group` (display label),
#'   and `Prevalence` (smoothed, in 0–1).
#'
#' @examples
#' # Species example (case-insensitive match), kid_serum, T2 vs T8 labeled B/M12:
#' prev_tbl <- ph_make_prev_tbl_standalone(
#'   x = ps_merged_box_bin,
#'   rank = "species",
#'   feature_value = "Acadevirus PM93",
#'   group_col = "timepoint_factor",
#'   groups = c("T2", "T8"),
#'   labels = c("B", "M12"),
#'   filter = list(big_group = "kid_serum")
#' )
#' # Now plot with your delta plotter:
#' # deltaplot_prevalence(prev_tbl, group_order = c("B","M12"))
#'
#' @export
ph_make_prev_tbl_standalone <- function(
  x,
  rank = "is_flagellum",
  feature_value = TRUE,
  group_col = "timepoint_factor",
  groups = c("T2", "T8"),
  labels = c("B", "M12"),
  filter = list(),
  eps = 0.5,
  den_mult = 2
) {
  stopifnot(length(groups) == 2L, length(labels) == 2L)
  if (!inherits(x, "phip_data")) {
    stop("`x` must be a phip_data with $data_long and $peptide_library.")
  }

  dl <- x$data_long
  need_cols <- c("sample_id", "peptide_id", group_col, "exist")
  miss <- setdiff(need_cols, colnames(dl))
  if (length(miss)) stop("x$data_long is missing: ", paste(miss, collapse = ", "))

  # 1) apply equality filters safely (SQL-friendly)
  if (length(filter)) {
    for (nm in names(filter)) {
      if (!nm %in% colnames(dl)) {
        stop("Filter column '", nm, "' not in x$data_long.")
      }
      vals <- filter[[nm]]
      dl <- dplyr::filter(dl, .data[[nm]] %in% !!vals)
    }
  }

  # target groups only
  dl <- dplyr::filter(dl, .data[[group_col]] %in% !!groups)

  # 2) ensure peptide_library is on same connection, join rank
  lib_src <- x$peptide_library %>%
    dplyr::select(peptide_id, !!rlang::sym(rank)) %>%
    dplyr::distinct()

  con_dl <- tryCatch(dbplyr::remote_con(dl), error = function(...) NULL)
  con_lib <- tryCatch(dbplyr::remote_con(lib_src), error = function(...) NULL)
  same_con <- !is.null(con_dl) && !is.null(con_lib) && identical(con_dl, con_lib)

  if (!same_con && !is.null(con_dl)) {
    lib_df <- dplyr::collect(lib_src)
    tmp <- paste0("ph_tmp_lib_", as.integer(runif(1) * 1e9))
    DBI::dbWriteTable(con_dl, tmp, tibble::as_tibble(lib_df), temporary = TRUE, overwrite = TRUE)
    lib_tbl <- dplyr::tbl(con_dl, tmp)
  } else {
    lib_tbl <- lib_src
  }

  dl <- dl %>%
    dplyr::left_join(lib_tbl, by = "peptide_id") %>%
    dplyr::mutate(
      rank_chr = dplyr::if_else(is.na(.data[[rank]]), NA_character_, as.character(.data[[rank]])),
      rank_ci  = tolower(rank_chr)
    )

  # 3) filter by feature_value (case-insensitive) on SQL
  fv_ci <- tolower(as.character(feature_value))
  dl <- dplyr::filter(dl, .data$rank_ci %in% !!fv_ci)

  # 4) subject fallback id (still SQL)
  dl <- dl %>%
    dplyr::mutate(
      subj_chr   = dplyr::if_else(!is.na(.data$subject_id), as.character(.data$subject_id), NA_character_),
      sample_chr = paste0("S", as.character(.data$sample_id)),
      id         = dplyr::coalesce(subj_chr, sample_chr)
    )

  # 5) compute x, n, prevalence on SQL ONLY — no labeling yet
  prev_sql <- dl %>%
    dplyr::group_by(peptide_id, !!rlang::sym(group_col)) %>%
    dplyr::summarise(
      x = sum(exist > 0L),
      n = dplyr::n_distinct(id),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Prevalence = (x + eps) / (n + den_mult * eps)) %>%
    dplyr::rename(Group_raw = !!rlang::sym(group_col)) %>%
    dplyr::select(peptide_id, Group_raw, Prevalence)

  # 6) now materialize, then map to labels in R (avoids dbplyr CASE/CAST issues)
  prev <- dplyr::collect(prev_sql)

  prev <- prev %>%
    dplyr::mutate(
      Group = dplyr::case_when(
        Group_raw == groups[1] ~ labels[1],
        Group_raw == groups[2] ~ labels[2],
        TRUE ~ as.character(Group_raw)
      )
    ) %>%
    dplyr::select(peptide_id, Group, Prevalence)

  # sanity: warn if unexpected groups made it through
  bad <- setdiff(unique(prev$Group), labels)
  if (length(bad)) {
    warning("Unexpected groups in output: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  prev
}

# -------------------------------------------------------------------------
# forest_delta: top/bottom forest plot for DELTA/Stouffer results (RAW T_obs)
# -------------------------------------------------------------------------

#' Forest plot of top/bottom raw Stouffer T by rank (with arrows)
#'
#' @description
#' Builds a forest plot for the most extreme DELTA/Stouffer results within a
#' chosen rank using the **raw** `T_obs` (no division by any SD). Rows can be
#' restricted to BH-significant results (per rank). Colors are purely aesthetic
#' and do not change geometry.
#'
#' @param x A data frame/tibble with at least:
#'   `rank, feature, group1, group2, design, T_obs, p_perm, p_adj_rank, category_rank_bh`.
#' @param rank_of_interest Character scalar, e.g. `"species"`.
#' @param n_each Integer; how many items from negative and positive tails (default 15).
#' @param color Logical; if `TRUE`, lines/points are shaded blue (left, T<0) to
#'   red (right, T>0) with higher contrast; otherwise monochrome.
#' @param bh_only Logical; if `TRUE` (default) keep only rows with
#'   `category_rank_bh == "significant (BH, per rank)"`.
#' @param add_signed_z_from_p Logical; if `TRUE`, computes a signed Z from
#'   permutation p (`sign(T_obs) * qnorm(1 - p_perm/2)`) and returns it in data.
#' @param left_label Text for the left arrow/side label.
#' @param right_label Text for the right arrow/side label.
#' @param arrow_frac Fraction of max |T_obs| used as half-length of arrows.
#' @param text_size Arrow-end label size (ggplot2 size units).
#' @param text_gap_frac Horizontal gap for arrow labels beyond arrow tips (fraction of max |T_obs|).
#' @param y_pad Vertical padding above the top category for arrows/labels (in y units).
#' @param text_y_offset Additional vertical offset for arrow-end labels (in y units).
#' @param text_vjust Passed to `annotate("text", vjust=...)` (default -0.30).
#' @param arrow_color Arrow/label color.
#' @param arrow_lwd Arrow line width.
#' @param arrow_head_cm Arrow head length in cm.
#'
#' @return A list with:
#' \itemize{
#'   \item `data` – tibble used for plotting (selected top/bottom items),
#'   \item `plot` – ggplot object.
#' }
#'
#' @details
#' - `T_obs` here is plotted **as is** (raw Stouffer statistic from your shift).
#' - For calibrated inference or cross-figure comparability, prefer permutation
#'   p-values (optionally expose `add_signed_z_from_p=TRUE` for reporting).
#'
#' @examples
#' \dontrun{
#' out <- forest_delta(res_shift,
#'   rank_of_interest = "species",
#'   n_each = 15, color = TRUE
#' )
#' out$plot
#' }
#'
#' @importFrom dplyr filter arrange slice_head bind_rows mutate case_when desc any_of
#' @importFrom ggplot2 ggplot aes geom_vline geom_segment geom_point labs theme
#' @importFrom ggplot2 annotate scale_color_gradientn coord_cartesian theme_void
#' @importFrom forcats fct_reorder
#' @importFrom stats qnorm
#' @export
forest_delta <- function(
  x,
  rank_of_interest,
  n_each = 15,
  color = FALSE,
  bh_only = TRUE,
  add_signed_z_from_p = FALSE,
  # arrow/label params
  left_label = "More in group1",
  right_label = "More in group2",
  arrow_frac = 0.35,
  text_size = 3.8,
  text_gap_frac = 0.06,
  y_pad = 0.6,
  text_y_offset = 0.00,
  text_vjust = -0.30,
  arrow_color = "red",
  arrow_lwd = 0.6,
  arrow_head_cm = 0.30
) {
  # ---- basic checks -----------------------------------------------------------
  if (!is.data.frame(x)) {
    if (exists(".ph_abort", mode = "function")) .ph_abort("`x` must be a data.frame/tibble.")
    stop("`x` must be a data.frame/tibble.")
  }
  need_cols <- c("rank", "feature", "group1", "group2", "design", "T_obs", "p_perm", "p_adj_rank", "category_rank_bh")
  miss <- setdiff(need_cols, colnames(x))
  if (length(miss)) {
    msg <- paste("`x` is missing required columns:", paste(miss, collapse = ", "))
    if (exists(".ph_abort", mode = "function")) .ph_abort(msg) else stop(msg)
  }

  # ---- simple timer + closing log --------------------------------------------
  t0 <- Sys.time()
  on.exit(
    {
      elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 3)
      if (exists(".ph_log_ok", mode = "function")) {
        # no 'step' param → nic nie dotyka nzchar(step)
        try(
          .ph_log_ok(
            headline = "forest_delta finished",
            bullets = sprintf("elapsed: %s s", elapsed)
          ),
          silent = TRUE
        )
      } else {
        message(sprintf("[forest_delta] elapsed: %0.3f s", elapsed))
      }
    },
    add = TRUE
  )

  # ---- optional signed Z from permutation p ----------------------------------
  if (isTRUE(add_signed_z_from_p)) {
    x$Z_signed_from_p_2s <- sign(x$T_obs) * stats::qnorm(1 - x$p_perm / 2)
    if (exists(".ph_log_info", mode = "function")) {
      try(
        .ph_log_info(
          headline = "added signed Z from permutation p",
          bullets = "column: Z_signed_from_p_2s"
        ),
        silent = TRUE
      )
    }
  }

  # ---- filter rank & BH if requested -----------------------------------------
  df_rk <- dplyr::filter(x, .data$rank == rank_of_interest)
  if (isTRUE(bh_only)) {
    df_rk <- dplyr::filter(df_rk, .data$category_rank_bh == "significant (BH, per rank)")
  }

  if (nrow(df_rk) == 0L) {
    if (exists(".ph_log_info", mode = "function")) {
      try(
        .ph_log_info(
          headline = "no rows left for plotting",
          bullets = c(
            paste0("rank: ", rank_of_interest),
            paste0("bh_only: ", bh_only)
          )
        ),
        silent = TRUE
      )
    }
    return(list(
      data = df_rk,
      plot = .empty_placeholder_plot(sprintf("No rows to plot (rank='%s')", rank_of_interest))
    ))
  }

  # ---- select tails on RAW T_obs ---------------------------------------------
  n_pos <- sum(df_rk$T_obs > 0, na.rm = TRUE)
  n_neg <- sum(df_rk$T_obs < 0, na.rm = TRUE)

  top_pos <- df_rk |>
    dplyr::arrange(dplyr::desc(.data$T_obs)) |>
    dplyr::slice_head(n = min(n_each, n_pos))

  top_neg <- df_rk |>
    dplyr::arrange(.data$T_obs) |>
    dplyr::slice_head(n = min(n_each, n_neg))

  df_plot <- dplyr::bind_rows(top_neg, top_pos) |>
    dplyr::mutate(
      species_label = .data$feature,
      species_label = forcats::fct_reorder(.data$species_label, .data$T_obs)
    )

  # ---- subtitle context -------------------------------------------------------
  g1 <- if (nrow(df_plot)) df_plot$group1[1] else ""
  g2 <- if (nrow(df_plot)) df_plot$group2[1] else ""
  des <- if (nrow(df_plot)) df_plot$design[1] else ""

  # ---- color normalization (cosmetic only) -----------------------------------
  max_neg <- max(abs(df_plot$T_obs[df_plot$T_obs < 0]), na.rm = TRUE)
  max_pos <- max(df_plot$T_obs[df_plot$T_obs > 0], na.rm = TRUE)
  if (!is.finite(max_neg) || max_neg == 0) max_neg <- max(abs(df_plot$T_obs), 1, na.rm = TRUE)
  if (!is.finite(max_pos) || max_pos == 0) max_pos <- max(abs(df_plot$T_obs), 1, na.rm = TRUE)

  df_plot$T_col <- dplyr::case_when(
    df_plot$T_obs < 0 ~ -abs(df_plot$T_obs) / max_neg,
    df_plot$T_obs > 0 ~ df_plot$T_obs / max_pos,
    TRUE ~ 0
  )
  gamma <- 0.85
  df_plot$T_col <- sign(df_plot$T_col) * (abs(df_plot$T_col))^gamma

  # ---- base plot --------------------------------------------------------------
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$T_obs, y = .data$species_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
    {
      if (isTRUE(color)) {
        ggplot2::geom_segment(
          ggplot2::aes(
            x = 0, xend = .data$T_obs, y = .data$species_label, yend = .data$species_label,
            color = .data$T_col
          ),
          linewidth = 1, alpha = 0.9, show.legend = FALSE
        )
      } else {
        ggplot2::geom_segment(
          ggplot2::aes(x = 0, xend = .data$T_obs, y = .data$species_label, yend = .data$species_label),
          linewidth = 0.8, alpha = 0.85
        )
      }
    } +
    {
      if (isTRUE(color)) {
        ggplot2::geom_point(ggplot2::aes(color = .data$T_col), size = 3.6, show.legend = FALSE)
      } else {
        ggplot2::geom_point(size = 2.2)
      }
    } +
    ggplot2::labs(
      title = sprintf("Top/Bottom raw Stouffer T - rank: %s", rank_of_interest),
      subtitle = sprintf(
        "Contrast: %s vs %s (%s); shown: %d neg + %d pos",
        g1, g2, des, nrow(top_neg), nrow(top_pos)
      ),
      x = "Stouffer T (raw)",
      y = NULL
    ) +
    theme_phip()

  if (isTRUE(color)) {
    p <- p + ggplot2::scale_color_gradientn(
      colors = c("#1f4e79", "#6f94c2", "#eeeeee", "#e7a39c", "#8e1b10"),
      values = c(0, 0.25, 0.5, 0.75, 1),
      limits = c(-1, 1)
    )
  }

  # ---- arrows sized by RAW |T_obs| -------------------------------------------
  lvl_n <- length(levels(df_plot$species_label))
  if (lvl_n == 0L) lvl_n <- length(unique(df_plot$species_label))
  y_top <- lvl_n + y_pad

  maxabs <- max(abs(df_plot$T_obs), na.rm = TRUE)
  if (!is.finite(maxabs) || maxabs == 0) maxabs <- 1

  x_off <- arrow_frac * maxabs
  x_left_text <- -x_off - text_gap_frac * maxabs
  x_right_text <- x_off + text_gap_frac * maxabs

  # rewrite subtitle header to reflect left/right labels
  if (!is.null(p$labels$subtitle)) {
    s <- p$labels$subtitle
    if (grepl("^Contrast:", s)) {
      s <- sub("^Contrast:\\s*.*?;", paste0("Contrast: ", left_label, " vs ", right_label, ";"), s)
    } else {
      s <- paste0("Contrast: ", left_label, " vs ", right_label, "; ", s)
    }
    p$labels$subtitle <- s
  }

  p <- p +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0.05, 0.10))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::annotate("segment",
      x = 0, xend = -x_off, y = y_top, yend = y_top,
      colour = arrow_color, linewidth = arrow_lwd,
      arrow = ggplot2::arrow(
        length = grid::unit(arrow_head_cm, "cm"),
        type = "closed", ends = "last"
      )
    ) +
    ggplot2::annotate("text",
      x = x_left_text, y = y_top + text_y_offset, label = left_label,
      size = text_size, colour = arrow_color, fontface = "bold",
      hjust = 1, vjust = text_vjust
    ) +
    ggplot2::annotate("segment",
      x = 0, xend = x_off, y = y_top, yend = y_top,
      colour = arrow_color, linewidth = arrow_lwd,
      arrow = ggplot2::arrow(
        length = grid::unit(arrow_head_cm, "cm"),
        type = "closed", ends = "last"
      )
    ) +
    ggplot2::annotate("text",
      x = x_right_text, y = y_top + text_y_offset, label = right_label,
      size = text_size, colour = arrow_color, fontface = "bold",
      hjust = 0, vjust = text_vjust
    )

  list(data = df_plot, plot = p)
}

# --- small internal helper: white placeholder with red label ------------------
.empty_placeholder_plot <- function(label = "No significant results") {
  ggplot2::ggplot() +
    ggplot2::geom_text(
      data = data.frame(x = 0.5, y = 0.5),
      ggplot2::aes(x = x, y = y),
      label = label, color = "red", fontface = "bold", size = 6
    ) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = NA))
}
