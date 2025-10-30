# ==============================================================================
#' Static scatterplot of percent1 vs percent2 from `ph_prevalence_compare()`
#'
#' @description
#' A ggplot2 clone of `scatter_interactive()` with the same filtering behavior
#' and default coloring rules:
#' - If `color_by = NULL`, points are colored by significance category.
#'   If `category_rank_wbh` exists, it is used directly.
#'   Otherwise, categories are derived with priority:
#'     p_adj_rank_wbh -> p_adj_rank -> p_raw (threshold = `alpha`).
#'   If only one category remains, it auto-falls back to p-value bins for
#'   informative coloring.
#' - If `color_by` is provided, it joins peptide metadata (like interactive).
#'
#' @return A ggplot object.
#' @export
# ==============================================================================
scatter_static <- function(df,
                           pair = NULL,
                           rank = NULL,
                           universe = NULL,
                           features = NULL,
                           features_regex = FALSE,
                           universe_regex = FALSE,
                           xlab = NULL,
                           ylab = NULL,
                           alpha = 0.05,
                           prefer_flags = TRUE,
                           color_by = NULL,
                           color_title = NULL) {

  # -- keep original so we can access prev_meta$peptide_library ---------------
  df_original <- df

  # -- ph_prev_result: same filtering behavior as scatter_interactive ----------
  if (inherits(df, "ph_prev_result")) {
    if (!is.null(pair)) {
      df <- prev_filter_pairs(
        df,
        gA = pair[1], gB = pair[2],
        ranks = rank,
        features = features, features_regex = features_regex,
        group_universe = universe, universe_regex = universe_regex,
        drop_na = TRUE
      )
    } else {
      df <- as.data.frame(df)
      if (!is.null(rank) && "rank" %in% names(df)) {
        df <- df[df$rank %in% rank, , drop = FALSE]
      }
      if (!is.null(universe) && "group_col" %in% names(df)) {
        gvec <- as.character(df$group_col); gvec[is.na(gvec)] <- ""
        if (isTRUE(universe_regex)) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(universe))
            keep <- keep | grepl(p, gvec, ignore.case = TRUE, perl = TRUE)
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(gvec) %in% tolower(as.character(universe)), , drop = FALSE]
        }
      }
      if (!is.null(features) && "feature" %in% names(df)) {
        fvec <- as.character(df$feature); fvec[is.na(fvec)] <- ""
        if (isTRUE(features_regex)) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(features))
            keep <- keep | grepl(p, fvec, ignore.case = TRUE, perl = TRUE)
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(fvec) %in% tolower(as.character(features)), , drop = FALSE]
        }
      }
    }
  } else {
    df <- as.data.frame(df)
  }

  # -- empty guard -------------------------------------------------------------
  if (is.null(df) || !nrow(df)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "no data for this contrast") +
        ggplot2::theme_void()
    )
  }

  # -- default axis labels from pair ------------------------------------------
  if (!is.null(pair)) {
    if (is.null(xlab)) xlab <- pair[1]
    if (is.null(ylab)) ylab <- pair[2]
  }
  if (is.null(xlab)) xlab <- "group a"
  if (is.null(ylab)) ylab <- "group b"

  # ---------- (a) significance category when color_by is NULL ----------------
  color_var <- NULL
  color_val_label <- NULL
  if (is.null(color_by)) {
    # Use existing category if present; otherwise derive from p-values
    has_cat <- "category_rank_wbh" %in% names(df) &&
      any(nzchar(trimws(as.character(df$category_rank_wbh %||% ""))))

    if (has_cat) {
      df <- df %>%
        dplyr::mutate(
          category_raw = tolower(trimws(as.character(.data$category_rank_wbh %||% NA_character_))),
          category = dplyr::case_when(
            grepl("wbh", category_raw) & grepl("significant", category_raw) ~ "significant (wBH, per rank)",
            grepl("nominal", category_raw)                                  ~ "nominal only",
            grepl("not significant", category_raw)                          ~ "not significant",
            TRUE                                                            ~ "not significant"
          )
        )
    } else {
      df <- df %>%
        dplyr::mutate(
          category = dplyr::case_when(
            !is.na(.data$p_adj_rank_wbh) & .data$p_adj_rank_wbh <= alpha ~ "significant (wBH, per rank)",
            !is.na(.data$p_raw)          & .data$p_raw          <= alpha ~ "nominal only",
            TRUE ~ "not significant"
          )
        )
    }

    # If there is only a single category, fall back to p-value bins for color
    cat_levels <- c("significant (wBH, per rank)",
                    "significant (BH, per rank)",
                    "nominal only",
                    "not significant")
    df$category <- factor(df$category, levels = cat_levels)

    if (length(stats::na.omit(unique(df$category))) <= 1L) {
      base_p <- if ("p_adj_rank_wbh" %in% names(df)) df$p_adj_rank_wbh else df$p_raw
      df$p_bin <- cut(base_p,
                      breaks = c(0, 1e-3, 1e-2, 5e-2, 1, Inf),
                      labels = c("≤1e-3", "(1e-3,1e-2]", "(1e-2,0.05]", ">0.05", "NA"),
                      include.lowest = TRUE, right = TRUE)
      color_var <- "p_bin"
      color_val_label <- "p-value"
    } else {
      color_var <- "category"
      color_val_label <- NULL
    }
  }

  # ---------- (b) peptide-level join for color_by using saved library ---------
  if (!is.null(color_by)) {
    # derive peptide key
    if ("peptide_id" %in% names(df)) {
      df <- dplyr::mutate(df, pep_key = as.character(.data$peptide_id))
    } else if (all(c("feature", "rank") %in% names(df))) {
      df <- dplyr::mutate(df, pep_key = dplyr::if_else(.data$rank == "peptide_id",
                                                       as.character(.data$feature),
                                                       NA_character_))
    } else if ("feature" %in% names(df)) {
      df <- dplyr::mutate(df, pep_key = as.character(.data$feature))
    } else {
      df <- dplyr::mutate(df, pep_key = NA_character_)
    }

    has_any_keys <- any(!is.na(df$pep_key))
    if (has_any_keys) {
      lib_handle <- tryCatch({
        if (inherits(df_original, "ph_prev_result"))
          attr(df_original, "prev_meta")$peptide_library
        else NULL
      }, error = function(...) NULL)

      pm <- if (!is.null(lib_handle)) {
        lib_handle %>%
          dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
          dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
          dplyr::collect()
      } else {
        get_peptide_meta() %>%
          dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
          dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
          dplyr::collect()
      }

      df <- df %>%
        dplyr::left_join(pm, by = dplyr::join_by(pep_key == peptide_id))

      if (color_by %in% names(df) && is.list(df[[color_by]])) {
        df[[color_by]] <- vapply(
          df[[color_by]],
          function(z) {
            if (length(z) == 0) return(NA_character_)
            if (is.atomic(z)) return(as.character(z)[1])
            as.character(z[[1]])
          },
          character(1)
        )
      }
    } else {
      df[[color_by]] <- NA
    }

    if (is.logical(df[[color_by]])) {
      df[[color_by]] <- dplyr::case_when(
        df[[color_by]] %in% TRUE  ~ "yes",
        df[[color_by]] %in% FALSE ~ "no",
        is.na(df[[color_by]])     ~ "na"
      )
    }
    df[[color_by]] <- as.factor(as.character(df[[color_by]]))

    color_var <- color_by
    if (is.null(color_title)) color_title <- color_by
    color_val_label <- color_title
  }

  # --- ensure ratio exists: derive from prop1/prop2 or create NA column -------
  if (!"ratio" %in% names(df)) {
    if (all(c("prop1","prop2") %in% names(df))) {
      prop1_eps <- ifelse(is.na(df$prop1) | df$prop1 <= 0, df$prop1 + 1e-12, df$prop1)
      prop2_eps <- ifelse(is.na(df$prop2) | df$prop2 <= 0, df$prop2 + 1e-12, df$prop2)
      df$ratio  <- prop1_eps / prop2_eps
    } else {
      df$ratio <- NA_real_
    }
  }

  # ---------- (c) build plot --------------------------------------------------
  p <- ggplot2::ggplot(df, ggplot2::aes(x = percent1, y = percent2))

  # point layer
  if (identical(color_var, "category")) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data$category), alpha = 0.85, size = 1.9) +
      ggplot2::scale_color_manual(
        values = c(
          "significant (wBH, per rank)" = "#e31a1c",
          "nominal only"                = "#1b9e77",
          "not significant"             = "#386cb0"
        ),
        name = NULL, drop = FALSE
      )
  } else if (!is.null(color_var)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color_var]]), alpha = 0.85, size = 1.9) +
      ggplot2::labs(color = color_val_label)
  } else {
    p <- p + ggplot2::geom_point(alpha = 0.85, size = 1.9, color = "grey35")
  }

  # diagonal and limits
  rng <- range(c(df$percent1, df$percent2), na.rm = TRUE)
  p <- p +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey60") +
    ggplot2::expand_limits(x = rng, y = rng) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(x = xlab, y = ylab)

  # aesthetics
  p <- p +
    theme_phip() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}

#' Volcano plots (log2 ratio vs -log10 p) from ph_prevalence_compare()
#'
#' @param comparison_tbl tibble or DuckDB lazy tbl from ph_prevalence_compare()
#' @param ranks character vector of ranks to plot (default: all present)
#' @param facet_pairs logical; TRUE (default) facets all (group1,group2) pairs per rank,
#'   FALSE returns multiple plots per pair
#' @param facet_nrow,facet_ncol integers controlling facet_wrap layout (optional)
#' @param sparse_rank character(1) rank whose "not significant" rows you want to downsample
#' @param sparse_drop_pct numeric in [0,1]; fraction of "not significant" rows to drop
#'   within sparse_rank (default 0 = keep all). All significant rows are kept.
#' @param sparse_seed integer or NULL; set to make downsampling reproducible
#' @param fc_cut numeric; vertical cutoff for |log2(ratio)| (default 1)
#' @param p_cut  numeric; horizontal cutoff for p-value (default 0.05)
#' @param significant_colors named vector for category colors
#' @param interactive TRUE (default) for plotly, FALSE for static ggplot
#'
#' @return A named list of plots. If facet_pairs=TRUE: one plot per rank.
#'         If facet_pairs=FALSE: for each rank a sub-list of plots per (group1 vs group2).
#' @export
ph_volcano <- function(
    comparison_tbl,
    ranks = NULL,
    facet_pairs = TRUE,
    facet_nrow = NULL,
    facet_ncol = NULL,
    sparse_rank = NULL,
    sparse_drop_pct = 0,
    sparse_seed = NULL,
    fc_cut = 1,
    p_cut = 0.05,
    significant_colors = c(
      "not significant"                 = "#386cb0",
      "significant prior correction"    = "#1b9e77",
      "significant post FDR correction" = "#e31a1c"
    ),
    interactive = TRUE) {
  .ph_with_timing(
    headline = "make_interactive_volcano",
    expr = {
      # required + optional
      required_cols <- c("rank", "feature", "group_col", "group1", "group2", "category")
      optional_cols <- c("ratio", "prop1", "prop2", "p_raw", "p_adj", "pvals_not_adj")

      # columns present (works for lazy)
      have <- tryCatch(colnames(comparison_tbl), error = function(...) character())
      if (length(have) == 0 && inherits(comparison_tbl, "tbl_sql")) have <- colnames(comparison_tbl)

      miss <- setdiff(required_cols, have)
      if (length(miss)) {
        .ph_abort("Comparison table is missing required columns", bullets = paste("-", miss))
      }

      # filter ranks (lazy-safe)
      if (!is.null(ranks) && length(ranks)) {
        comparison_tbl <- dplyr::filter(comparison_tbl, rank %in% ranks)
      }

      sel_cols <- c(required_cols, intersect(optional_cols, have))
      df <- comparison_tbl |>
        dplyr::select(tidyselect::any_of(sel_cols)) |>
        dplyr::collect()

      if (!nrow(df)) {
        .ph_warn("No rows after filtering; returning empty list.")
        return(invisible(list()))
      }

      df <- df |>
        dplyr::mutate(
          rank = as.character(rank),
          group_col = as.character(group_col),
          group1 = as.character(group1),
          group2 = as.character(group2),
          category = as.character(category),
          pair_label = paste0(group1, " vs ", group2)
        )

      # ratio fallback
      if (!"ratio" %in% colnames(df)) {
        if (!all(c("prop1", "prop2") %in% colnames(df))) {
          .ph_abort("Need either 'ratio' or both 'prop1' and 'prop2' to build volcano.")
        }
        df <- df |>
          dplyr::mutate(
            prop1_eps = dplyr::if_else(is.na(prop1) | prop1 <= 0, prop1 + 1e-12, prop1),
            prop2_eps = dplyr::if_else(is.na(prop2) | prop2 <= 0, prop2 + 1e-12, prop2),
            ratio     = prop1_eps / prop2_eps
          )
      }

      # p column selection; clamp at machine eps to avoid Inf
      p_col <- intersect(c("p_raw", "pvals_not_adj", "p_adj"), colnames(df))[1]
      if (is.na(p_col)) .ph_abort("No p-value column found (need one of: p_raw, pvals_not_adj, p_adj).")

      df <- df |>
        dplyr::mutate(
          log2ratio = log2(ratio),
          p_use_raw = .data[[p_col]],
          p_use     = pmax(p_use_raw, .Machine$double.xmin),
          nlog10p   = -log10(p_use)
        )

      if (!any(is.finite(df$log2ratio)) || !any(is.finite(df$nlog10p))) {
        .ph_abort("Non-finite values for log2(ratio) or -log10(p). Check input.")
      }

      # downsample not significant in chosen rank
      if (!is.null(sparse_rank) && sparse_drop_pct > 0) {
        if (!sparse_rank %in% unique(df$rank)) {
          .ph_warn(paste0("sparse_rank '", sparse_rank, "' not found; skipping downsampling"))
        } else {
          if (!is.null(sparse_seed)) set.seed(as.integer(sparse_seed))
          pool_idx <- which(df$rank == sparse_rank & df$category == "not significant")
          n_pool <- length(pool_idx)
          if (n_pool > 0) {
            n_drop <- floor(sparse_drop_pct * n_pool)
            if (n_drop > 0) {
              drop_idx <- sample(pool_idx, size = n_drop, replace = FALSE)
              df <- df[-drop_idx, , drop = FALSE]
              .ph_log_info("Downsampled 'not significant' rows",
                           bullets = c(
                             paste0("rank: ", sparse_rank),
                             paste0("pool: ", n_pool),
                             paste0("dropped: ", n_drop),
                             paste0("kept: ", n_pool - n_drop)
                           )
              )
            }
          }
        }
      }

      # tooltips
      df <- df |>
        dplyr::mutate(
          tooltip_txt = paste0(
            "Feature: ", feature, "<br>",
            "Rank: ", rank, "<br>",
            "log2(ratio): ", signif(log2ratio, 3), "<br>",
            "-log10(p): ", signif(nlog10p, 3), "<br>",
            "Pair: ", pair_label, "<br>",
            "Category: ", category
          )
        )

      # color scale helper
      add_color_scale <- function() {
        if (is.null(significant_colors)) {
          scale_color_phip(name = NULL)
        } else {
          ggplot2::scale_color_manual(values = significant_colors, name = NULL)
        }
      }

      # builder per rank
      build_rank_plot <- function(dat_rank) {
        base <- ggplot2::ggplot(dat_rank, ggplot2::aes(x = log2ratio, y = nlog10p)) +
          {
            if (interactive) {
              ggplot2::geom_point(ggplot2::aes(color = category, text = tooltip_txt), alpha = 0.7)
            } else {
              ggplot2::geom_point(ggplot2::aes(color = category), alpha = 0.7)
            }
          } +
          add_color_scale() +
          ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "gray50") +
          ggplot2::geom_vline(xintercept = c(fc_cut, -fc_cut), linetype = "dashed", color = "gray50") +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.6) +
          ggplot2::labs(
            x = "log₂ ratio (group1 / group2)",
            y = "-log₁₀(p-value)",
            title = paste0("Rank: ", unique(dat_rank$rank)[1])
          ) +
          theme_phip(base_size = 12) +
          ggplot2::theme(legend.position = "none")

        if (facet_pairs) {
          base <- base + ggplot2::facet_wrap(~pair_label,
                                             scales = "fixed",
                                             nrow = facet_nrow, ncol = facet_ncol
          )
          if (interactive) {
            p <- plotly::ggplotly(base, tooltip = c("text", "x", "y"))
            p <- p |>
              plotly::layout(
                showlegend = FALSE,
                xaxis = list(automargin = TRUE),
                yaxis = list(automargin = TRUE, title = list(standoff = 10)),
                margin = list(l = 60, r = 40, b = 60, t = 60)
              )
            return(p)
          } else {
            return(base)
          }
        } else {
          pieces <- split(dat_rank, dat_rank$pair_label, drop = TRUE)
          out <- lapply(names(pieces), function(lbl) {
            subd <- pieces[[lbl]]
            p <- base %+% subd + ggplot2::ggtitle(paste0("Rank: ", unique(subd$rank)[1], " • ", lbl))
            if (interactive) {
              plt <- plotly::ggplotly(p, tooltip = c("text", "x", "y"))
              plt <- plt |>
                plotly::layout(
                  showlegend = FALSE,
                  xaxis = list(automargin = TRUE),
                  yaxis = list(automargin = TRUE, title = list(standoff = 10))
                )
              plt
            } else {
              p
            }
          })
          names(out) <- names(pieces)
          out
        }
      }

      .ph_log_info("Building volcano plots per rank")
      by_rank <- split(df, df$rank, drop = TRUE)
      plots <- lapply(names(by_rank), function(rk) {
        dat <- by_rank[[rk]]
        if (!nrow(dat)) {
          .ph_warn(paste0("No rows for rank '", rk, "'"))
          return(NULL)
        }
        build_rank_plot(dat)
      })
      names(plots) <- names(by_rank)

      .ph_log_ok("Volcano plots built", bullets = c(
        paste0("ranks: ", paste(names(plots), collapse = ", ")),
        paste0("facet_pairs: ", facet_pairs),
        paste0("fc_cut: ", fc_cut, ", p_cut: ", p_cut),
        if (!is.null(sparse_rank)) paste0("sparse_rank: ", sparse_rank, " (drop=", sparse_drop_pct, ")") else NULL
      )[!vapply(c(
        paste0("ranks: ", paste(names(plots), collapse = ", ")),
        paste0("facet_pairs: ", facet_pairs),
        paste0("fc_cut: ", fc_cut, ", p_cut: ", p_cut),
        if (!is.null(sparse_rank)) paste0("sparse_rank: ", sparse_rank, " (drop=", sparse_drop_pct, ")") else NA_character_
      ), is.na, logical(1))])

      plots
    }
  )
}


# ==============================================================================
#' Interactive prevalence scatter for `ph_prev_result`
#'
#' @description
#' creates an interactive scatter (plotly) comparing prevalence in **group a**
#' vs **group b** for per-feature results produced by `ph_prevalence_compare()`.
#' accepts either a `ph_prev_result` (recommended) or a plain data.frame with the
#' same columns (`percent1`, `percent2`, `feature`, `group1`, `group2`, etc.).
#'
#' when a `ph_prev_result` is passed, you may specify a concrete pair of group
#' levels via `pair = c("A","B")`, and optionally restrict to a `rank`, a
#' `universe` (`group_col`) and/or `features`. internally, subsetting uses
#' `prev_filter_pairs()`, preserving metadata.
#'
#' color mapping:
#' - by default (when `color_by = NULL`), points are colored by `category_rank_wbh`
#'   with three levels: "significant (wBH, per rank)", "nominal only",
#'   "not significant".
#' - when `color_by` is a peptide-level meta column name (e.g. `"is_flagellum"`),
#'   the function tries to join peptide metadata using `peptide_id` (or `feature`
#'   if `rank == "peptide_id"`), and color by that variable.
#'
#' @param df a `ph_prev_result` or a `data.frame` with columns:
#'   `percent1`, `percent2`, `feature`, `group1`, `group2`,
#'   `p_raw`, `p_adj_rank`, `p_adj_rank_wbh`, `category_rank_wbh`,
#'   optionally `n_peptides`, `rank`, `peptide_id`, `group_col`.
#' @param pair optional length-2 character, e.g. `c("kid_serum::T2","kid_serum::T8")`.
#'   only used when `df` is a `ph_prev_result`; filters rows to that exact pair
#'   (order-agnostic).
#' @param rank optional single rank (character) to keep (only used with `ph_prev_result`).
#' @param universe optional `group_col` value or regex (if `universe_regex = TRUE`);
#'   only used with `ph_prev_result`.
#' @param features optional character vector or regex patterns (if `features_regex = TRUE`);
#'   only used with `ph_prev_result`.
#' @param features_regex logical; treat `features` as regex patterns (or).
#' @param universe_regex logical; treat `universe` as regex pattern(s) (or).
#' @param xlab,ylab axis labels; if missing and `pair` is provided, they default
#'   to `pair[1]` and `pair[2]`.
#' @param alpha numeric in (0,1]; used only for nominal labels; not the plotly alpha.
#' @param prefer_flags logical; reserved for future use (kept for back-compat).
#' @param color_by optional peptide-level meta column name to color by.
#' @param color_title optional legend title for `color_by`.
#'
#' @return a `plotly` object.
#'
#' @examples
#' # typical usage with ph_prev_result:
#' # p <- scatter_interactive(scatters,
#' #   pair     = c("kid_serum::T2","kid_serum::T8"),
#' #   rank     = "peptide_id",
#' #   color_by = "is_flagellum",
#' #   color_title = "Flagellum"
#' # )
#' # p
#'
#' @export
# ==============================================================================
scatter_interactive <- function(df,
                                pair = NULL,
                                rank = NULL,
                                universe = NULL,
                                features = NULL,
                                features_regex = FALSE,
                                universe_regex = FALSE,
                                xlab = NULL,
                                ylab = NULL,
                                alpha = 0.05,
                                prefer_flags = TRUE,
                                color_by = NULL,
                                color_title = NULL) {

  # -- keep original to access prev_meta (peptide library handle) -------------
  df_original <- df

  # -- ph_prev_result: optional filtering via prev_filter_pairs ----------------
  if (inherits(df, "ph_prev_result")) {
    if (!is.null(pair)) {
      df <- prev_filter_pairs(
        df,
        gA = pair[1], gB = pair[2],
        ranks = rank,
        features = features, features_regex = features_regex,
        group_universe = universe, universe_regex = universe_regex,
        drop_na = TRUE
      )
    } else {
      df <- as.data.frame(df)
      # optional structural filters without requiring a pair
      if (!is.null(rank) && "rank" %in% names(df)) {
        df <- df[df$rank %in% rank, , drop = FALSE]
      }
      if (!is.null(universe) && "group_col" %in% names(df)) {
        gvec <- as.character(df$group_col); gvec[is.na(gvec)] <- ""
        if (universe_regex) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(universe)) keep <- keep | grepl(p, gvec, ignore.case = TRUE, perl = TRUE)
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(gvec) %in% tolower(as.character(universe)), , drop = FALSE]
        }
      }
      if (!is.null(features) && "feature" %in% names(df)) {
        fvec <- as.character(df$feature); fvec[is.na(fvec)] <- ""
        if (features_regex) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(features)) keep <- keep | grepl(p, fvec, ignore.case = TRUE, perl = TRUE)
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(fvec) %in% tolower(as.character(features)), , drop = FALSE]
        }
      }
    }
  } else {
    df <- as.data.frame(df)
  }

  # -- empty guard -------------------------------------------------------------
  if (is.null(df) || !nrow(df)) {
    p <- plotly::plot_ly()
    return(plotly::layout(
      p,
      annotations = list(x = 0.5, y = 0.5, text = "no data for this contrast",
                         showarrow = FALSE, xref = "paper", yref = "paper"),
      xaxis = list(title = xlab %||% "group a", zeroline = FALSE),
      yaxis = list(title = ylab %||% "group b", zeroline = FALSE),
      legend = list(orientation = "h", y = -0.15)
    ))
  }

  # -- default axis labels from pair ------------------------------------------
  if (!is.null(pair)) {
    if (is.null(xlab)) xlab <- pair[1]
    if (is.null(ylab)) ylab <- pair[2]
  }
  if (is.null(xlab)) xlab <- "group a"
  if (is.null(ylab)) ylab <- "group b"

  # ---------- (a) significance category when color_by is null -----------------
  cat_levels <- c("significant (wBH, per rank)", "nominal only", "not significant")
  if (is.null(color_by)) {
    df <- df %>%
      dplyr::mutate(
        category_raw = tolower(trimws(as.character(.data$category_rank_wbh %||% NA_character_))),
        category = dplyr::case_when(
          grepl("wbh", category_raw) & grepl("significant", category_raw) ~ cat_levels[1],
          grepl("nominal", category_raw)                                  ~ cat_levels[2],
          grepl("not significant", category_raw)                          ~ cat_levels[3],
          TRUE                                                            ~ "other/NA"
        ),
        category = factor(category, levels = c(cat_levels, "other/NA"))
      )
  }

  # ---------- (b) peptide-level join for color_by using saved library ---------
  color_var <- NULL
  color_val_label <- NULL

  if (!is.null(color_by)) {
    # derive peptide key
    if ("peptide_id" %in% names(df)) {
      df <- dplyr::mutate(df, pep_key = as.character(.data$peptide_id))
    } else if (all(c("feature", "rank") %in% names(df))) {
      df <- dplyr::mutate(df, pep_key = dplyr::if_else(.data$rank == "peptide_id",
                                                       as.character(.data$feature),
                                                       NA_character_))
    } else if ("feature" %in% names(df)) {
      df <- dplyr::mutate(df, pep_key = as.character(.data$feature))
    } else {
      df <- dplyr::mutate(df, pep_key = NA_character_)
    }

    has_any_keys <- any(!is.na(df$pep_key))
    if (has_any_keys) {
      # try the peptide library saved in prev_meta
      lib_handle <- NULL
      if (inherits(df_original, "ph_prev_result")) {
        meta <- attr(df_original, "prev_meta")
        if (!is.null(meta) && "peptide_library" %in% names(meta)) {
          lib_handle <- meta$peptide_library
        }
      }

      if (!is.null(lib_handle)) {
        # use the saved duckdb/dbplyr handle without downloading
        pm <- lib_handle %>%
          dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
          dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
          dplyr::collect()
      } else {
        # fallback to global provider
        pm <- get_peptide_meta() %>%
          dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
          dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
          dplyr::collect()
      }

      df <- df %>%
        dplyr::left_join(pm, by = dplyr::join_by(pep_key == peptide_id))

      # flatten list-cols if any
      if (color_by %in% names(df) && is.list(df[[color_by]])) {
        df[[color_by]] <- vapply(
          df[[color_by]],
          function(z) {
            if (length(z) == 0) return(NA_character_)
            if (is.atomic(z)) return(as.character(z)[1])
            as.character(z[[1]])
          },
          character(1)
        )
      }
    } else {
      df[[color_by]] <- NA
    }

    # normalize logicals and coerce to factor (atomic)
    if (is.logical(df[[color_by]])) {
      df[[color_by]] <- dplyr::case_when(
        df[[color_by]] %in% TRUE  ~ "yes",
        df[[color_by]] %in% FALSE ~ "no",
        is.na(df[[color_by]])     ~ "na"
      )
    }
    df[[color_by]] <- as.factor(as.character(df[[color_by]]))

    color_var <- color_by
    if (is.null(color_title)) color_title <- color_by
    color_val_label <- color_title
  }

  # --- ensure ratio exists: derive from prop1/prop2 or create NA column -------
  if (!"ratio" %in% names(df)) {
    if (all(c("prop1","prop2") %in% names(df))) {
      prop1_eps <- ifelse(is.na(df$prop1) | df$prop1 <= 0, df$prop1 + 1e-12, df$prop1)
      prop2_eps <- ifelse(is.na(df$prop2) | df$prop2 <= 0, df$prop2 + 1e-12, df$prop2)
      df$ratio  <- prop1_eps / prop2_eps
    } else {
      df$ratio <- NA_real_
    }
  }

  # ---------- (c) hover fields and numerics ----------------------------------
  df <- df %>%
    dplyr::mutate(
      p1r  = round(.data$percent1, 2),
      p2r  = round(.data$percent2, 2),
      rr   = dplyr::coalesce(round(.data$ratio, 2), NA_real_),
      praw = round(.data$p_raw, 3),
      padj_wbh = round(.data$p_adj_rank_wbh, 3),
      padj_bh  = if ("p_adj_rank" %in% names(df)) round(.data$p_adj_rank, 3) else NA_real_,
      padj_bh_str = ifelse(is.na(.data$padj_bh), "NA", sprintf("%.3f", .data$padj_bh)),
      npep = dplyr::coalesce(
        if ("n_peptides" %in% names(df)) .data[["n_peptides"]] else NULL,
        if ("n_peptide"  %in% names(df)) .data[["n_peptide" ]] else NULL
      ),
      color_info = if (!is.null(color_var)) as.character(.data[[color_var]]) else NA_character_
    )

  df <- df %>%
    dplyr::mutate(
      text = sprintf(
        "<b>%s</b>\
<br>%s: %d/%d (%.2f%%)\
<br>%s: %d/%d (%.2f%%)%s\
<br>peptides: %s\
<br>ratio: %s\
<br>p: %.3f\
<br>p_adj (bh, per rank): %s\
<br>p_adj (wbh, per rank): %.3f",
        .data$feature,
        xlab, .data$n1, .data$N1, .data$p1r,
        ylab, .data$n2, .data$N2, .data$p2r,
        if (!is.null(color_var)) sprintf("<br>%s: %s", color_val_label, .data$color_info) else "",
        format(.data$npep, big.mark = "\u2009", scientific = FALSE),
        ifelse(is.na(.data$rr), "NA", sprintf("%.2f", .data$rr)),
        .data$praw, .data$padj_bh_str, .data$padj_wbh
      )
    )

  # ---------- (d) plotly scatter ---------------------------------------------
  color_mapping <- if (is.null(color_var)) {
    ~category
  } else {
    stats::as.formula(paste0("~`", color_var, "`"))
  }

  p <- plotly::plot_ly(
    df,
    x = ~percent1, y = ~percent2,
    color = color_mapping,
    type = "scatter", mode = "markers",
    text = ~text,
    hovertemplate = "%{text}<extra></extra>",
    marker = list(size = 7, opacity = 0.85)
  )

  rng <- range(c(df$percent1, df$percent2), na.rm = TRUE)
  p <- plotly::add_lines(p, x = rng, y = rng, inherit = FALSE, showlegend = FALSE)

  leg <- list(orientation = "h", y = -0.15)
  if (!is.null(color_var) && !is.null(color_title)) {
    leg$title <- list(text = color_title)
  }

  plotly::layout(
    p,
    xaxis = list(title = xlab, zeroline = FALSE),
    yaxis = list(title = ylab, zeroline = FALSE),
    legend = leg
  )
}


################################################################################
# Definitely to work on --> doesnt plot as intended. Especially in paired
# designs
# ==============================================================================
#' Interactive volcano plots from `ph_prevalence_compare()` results
#'
#' @description
#' Build volcano plot(s) of log2(ratio) vs -log10(p) from results produced by
#' `ph_prevalence_compare()`. Accepts a `ph_prev_result` (recommended) or a
#' plain `data.frame` with the standard columns. Supports convenient subsetting
#' by pair, rank, universe (`group_col`) and features. By default, colors show
#' significance category; alternatively you can color by any peptide-level
#' attribute available in the **saved peptide library** attached in
#' `attr(x, "prev_meta")$peptide_library`.
#'
#' @details
#' **Which p-value is on the Y axis**
#' - `p_mode = "raw"`  uses `p_raw`
#' - `p_mode = "bh"`   uses `p_adj_rank`
#' - `p_mode = "wbh"`  uses `p_adj_rank_wbh`
#'
#' **Categories (default coloring)**
#' If `category_rank_wbh` exists, it's used. Otherwise categories are derived:
#' - "significant post fdr correction": chosen p (`p_mode`) ≤ `p_cut` **and** |log2ratio| ≥ `fc_cut`
#' - "significant prior correction"   : raw p ≤ `p_cut` **and** |log2ratio| ≥ `fc_cut` but adjusted > `p_cut`
#' - "not significant"                : otherwise
#'
#' **Coloring by peptide metadata**
#' Set `color_by = "<column>"` to color points by a column from the saved
#' peptide library (e.g. `is_flagellum`, `kingdom`, `domain`, ...). No download
#' is attempted; the function uses the handle stored in `prev_meta`.
#'
#' **Output**
#' - If `facet_pairs = TRUE`: one plot per rank (as a named list). For
#'   interactive plots, pairs are tiled using `plotly::subplot`.
#' - If `facet_pairs = FALSE`: per rank, a named sub-list with one plot per
#'   `(group1 vs group2)`.
#'
#' @param comparison_tbl `ph_prev_result` or `data.frame` with columns:
#'   `rank`, `feature`, `group_col`, `group1`, `group2`,
#'   `prop1`, `prop2`, optional `ratio`,
#'   `p_raw`, `p_adj_rank`, `p_adj_rank_wbh`,
#'   and optionally `category_rank_wbh`.
#' @param ranks character vector of ranks to plot (default: all present)
#' @param pair optional length-2 character, e.g. `c("kid_serum::T2","kid_serum::T8")`
#' @param universe optional `group_col` value(s) or regex (with `universe_regex=TRUE`)
#' @param features optional character vector or regex patterns (with `features_regex=TRUE`)
#' @param features_regex logical; treat `features` as regex (OR over patterns)
#' @param universe_regex logical; treat `universe` as regex (OR over patterns)
#' @param color_by optional peptide-library column to color by (e.g. `"is_flagellum"`)
#' @param color_title optional legend title when `color_by` is used
#' @param facet_pairs logical; if `TRUE` (default) facet/subplot all pairs per rank
#' @param facet_nrow,facet_ncol integers controlling subplot grid (interactive only)
#' @param sparse_rank character(1) rank within which to downsample "not significant"
#' @param sparse_drop_pct numeric in [0,1]; fraction to drop within `sparse_rank`
#' @param sparse_seed integer or `NULL`; set for reproducible downsampling
#' @param fc_cut numeric; vertical cutoff for `|log2(ratio)|` (default 1)
#' @param p_cut numeric; cutoff for p/q (default 0.05)
#' @param p_mode one of `c("raw","bh","wbh")` (default `"raw"`)
#' @param significant_colors named vector used for default category colors
#' @param interactive `TRUE` (default) uses plotly; `FALSE` returns ggplot
#'
#' @return named list of plots:
#'   - if `facet_pairs=TRUE`: one plot per rank
#'   - if `facet_pairs=FALSE`: per rank a sub-list of `(group1 vs group2)` plots
#'
#' @examples
#' # p <- ph_volcano(
#' #   scatters,
#' #   ranks = "peptide_id",
#' #   pair  = c("kid_serum::T2","kid_serum::T8"),
#' #   p_mode = "wbh", fc_cut = 1, p_cut = 0.05,
#' #   color_by = "is_flagellum", color_title = "Flagellum"
#' # )
#' # p[["peptide_id"]]
#' @export
# ==============================================================================
ph_volcano <- function(
    comparison_tbl,
    ranks = NULL,
    pair = NULL,
    universe = NULL,
    features = NULL,
    features_regex = FALSE,
    universe_regex = FALSE,
    color_by = NULL,
    color_title = NULL,
    facet_pairs = TRUE,
    facet_nrow = NULL,
    facet_ncol = NULL,
    sparse_rank = NULL,
    sparse_drop_pct = 0,
    sparse_seed = NULL,
    fc_cut = 1,
    p_cut = 0.05,
    p_mode = c("raw","bh","wbh"),
    significant_colors = c(
      "not significant"                 = "#386cb0",
      "significant prior correction"    = "#1b9e77",
      "significant post fdr correction" = "#e31a1c"
    ),
    interactive = TRUE
) {
  p_mode <- match.arg(p_mode)

  # ---------- (0) normalize / subset input (mirror scatter_interactive) ------
  df_in <- comparison_tbl
  if (inherits(df_in, "ph_prev_result")) {
    if (!is.null(pair)) {
      df <- prev_filter_pairs(
        df_in,
        gA = pair[1], gB = pair[2],
        ranks = ranks,
        features = features, features_regex = features_regex,
        group_universe = universe, universe_regex = universe_regex,
        drop_na = TRUE
      )
    } else {
      df <- as.data.frame(df_in)
      if (!is.null(ranks) && "rank" %in% names(df)) {
        df <- df[df$rank %in% ranks, , drop = FALSE]
      }
      if (!is.null(universe) && "group_col" %in% names(df)) {
        gvec <- as.character(df$group_col); gvec[is.na(gvec)] <- ""
        if (universe_regex) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(universe)) keep <- keep | grepl(p, gvec, ignore.case = TRUE, perl = TRUE)
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(gvec) %in% tolower(as.character(universe)), , drop = FALSE]
        }
      }
      if (!is.null(features) && "feature" %in% names(df)) {
        fvec <- as.character(df$feature); fvec[is.na(fvec)] <- ""
        if (features_regex) {
          keep <- rep(FALSE, nrow(df))
          for (p in as.character(features)) keep <- keep | grepl(p, fvec, ignore.case = TRUE, perl = TRUE)
          df <- df[keep, , drop = FALSE]
        } else {
          df <- df[tolower(fvec) %in% tolower(as.character(features)), , drop = FALSE]
        }
      }
    }
  } else {
    df <- as.data.frame(df_in)
    if (!is.null(ranks) && "rank" %in% names(df)) df <- df[df$rank %in% ranks, , drop = FALSE]
  }

  if (is.null(df) || !nrow(df)) return(invisible(list()))

  # ---------- (1) required columns, minimal coercions ------------------------
  req <- c("rank","feature","group1","group2")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("ph_volcano: missing required columns: ", paste(miss, collapse = ", "))

  df$rank    <- as.character(df$rank)
  df$feature <- as.character(df$feature)
  df$group1  <- as.character(df$group1)
  df$group2  <- as.character(df$group2)
  if (!"group_col" %in% names(df)) df$group_col <- NA_character_
  df$pair_label <- paste0(df$group1, " vs ", df$group2)

  # ---------- (2) ensure ratio exists (fallback from prop1/prop2) ------------
  if (!"ratio" %in% names(df)) {
    if (!all(c("prop1","prop2") %in% names(df)))
      stop("ph_volcano: need either 'ratio' or both 'prop1' and 'prop2'.")
    p1 <- ifelse(is.na(df$prop1) | df$prop1 <= 0, df$prop1 + 1e-12, df$prop1)
    p2 <- ifelse(is.na(df$prop2) | df$prop2 <= 0, df$prop2 + 1e-12, df$prop2)
    df$ratio <- p1 / p2
  }

  # ---------- (3) choose p column and compute axes ---------------------------
  p_col <- switch(p_mode, raw = "p_raw", bh = "p_adj_rank", wbh = "p_adj_rank_wbh")
  if (!p_col %in% names(df)) stop("ph_volcano: requested p_mode='", p_mode, "' but column '", p_col, "' not found.")
  df$log2ratio <- log2(df$ratio)
  p_use <- pmax(df[[p_col]], .Machine$double.xmin)
  df$nlog10p <- -log10(p_use)

  # ---------- (4) category coloring (fallback if color_by is NULL) -----------
  if (is.null(color_by)) {
    if ("category_rank_wbh" %in% names(df)) {
      raw <- tolower(trimws(as.character(df$category_rank_wbh)))
      cat_use <- ifelse(grepl("significant", raw) & grepl("wbh", raw),
                        "significant post fdr correction",
                        ifelse(grepl("nominal", raw),
                               "significant prior correction",
                               "not significant"))
    } else {
      rawp  <- if ("p_raw" %in% names(df)) df$p_raw else p_use
      sigfc <- abs(df$log2ratio) >= fc_cut
      post  <- (p_use <= p_cut) & sigfc
      prior <- (rawp  <= p_cut) & sigfc & !post
      cat_use <- ifelse(post,  "significant post fdr correction",
                        ifelse(prior,"significant prior correction","not significant"))
    }
    df$category_use <- factor(cat_use,
                              levels = c("significant post fdr correction",
                                         "significant prior correction",
                                         "not significant"))
  }

  # ---------- (5) color by peptide metadata (use saved library handle) -------
  color_var <- NULL
  color_label <- NULL
  if (!is.null(color_by)) {
    if ("peptide_id" %in% names(df)) {
      df$pep_key <- as.character(df$peptide_id)
    } else if (all(c("feature","rank") %in% names(df))) {
      df$pep_key <- ifelse(df$rank == "peptide_id", as.character(df$feature), NA_character_)
    } else if ("feature" %in% names(df)) {
      df$pep_key <- as.character(df$feature)
    } else {
      df$pep_key <- NA_character_
    }

    has_keys <- any(!is.na(df$pep_key))
    if (has_keys) {
      lib_handle <- tryCatch(attr(df_in, "prev_meta")$peptide_library, error = function(...) NULL)
      if (!is.null(lib_handle)) {
        pm <- lib_handle %>%
          dplyr::select("peptide_id", tidyselect::all_of(color_by)) %>%
          dplyr::distinct(peptide_id, .keep_all = TRUE) %>%
          dplyr::collect()
        df <- dplyr::left_join(df, pm, by = dplyr::join_by(pep_key == peptide_id))

        # flatten potential list-cols and coerce logicals to friendly labels
        if (color_by %in% names(df) && is.list(df[[color_by]])) {
          df[[color_by]] <- vapply(
            df[[color_by]],
            function(z) {
              if (length(z) == 0) return(NA_character_)
              if (is.atomic(z)) return(as.character(z)[1])
              as.character(z[[1]])
            },
            character(1)
          )
        }
        if (is.logical(df[[color_by]])) {
          df[[color_by]] <- dplyr::case_when(
            df[[color_by]] %in% TRUE  ~ "Yes",
            df[[color_by]] %in% FALSE ~ "No",
            is.na(df[[color_by]])     ~ "NA"
          )
        }
      } else {
        df[[color_by]] <- NA
      }
    } else {
      df[[color_by]] <- NA
    }

    df[[color_by]] <- as.factor(as.character(df[[color_by]]))
    color_var   <- color_by
    color_label <- if (is.null(color_title)) color_by else color_title
  }

  # ---------- (6) optional downsampling of "not significant" -----------------
  if (!is.null(sparse_rank) && sparse_drop_pct > 0) {
    if (sparse_rank %in% unique(df$rank)) {
      pool_idx <- which(df$rank == sparse_rank &
                          (if (is.null(color_var)) df$category_use == "not significant" else TRUE))
      if (length(pool_idx)) {
        if (!is.null(sparse_seed)) set.seed(as.integer(sparse_seed))
        n_drop <- floor(sparse_drop_pct * length(pool_idx))
        if (n_drop > 0) {
          drop_idx <- sample(pool_idx, size = n_drop, replace = FALSE)
          df <- df[-drop_idx, , drop = FALSE]
        }
      }
    }
  }

  # ---------- (7) split per rank and builders --------------------------------
  by_rank <- split(df, df$rank, drop = TRUE)

  build_one_plotly <- function(dd, title_suffix = NULL) {
    hov <- sprintf(
      "<b>%s</b><br>rank: %s<br>pair: %s vs %s<br>log2(ratio): %.3f<br>-log10(p): %.3f%s",
      dd$feature, dd$rank, dd$group1, dd$group2, dd$log2ratio, dd$nlog10p,
      if (!is.null(color_var))
        sprintf("<br>%s: %s", color_label, as.character(dd[[color_var]]))
      else ""
    )

    # color mapping (meta or categories)
    if (!is.null(color_var)) {
      pal <- grDevices::hcl.colors(n = max(3L, length(levels(dd[[color_var]]))), palette = "Dark 3")
      pal <- setNames(pal[seq_along(levels(dd[[color_var]]))], levels(dd[[color_var]]))
      col_vec <- pal[as.character(dd[[color_var]])]
    } else {
      pal <- significant_colors
      cat_vec <- as.character(dd$category_use)
      col_vec <- pal[cat_vec]; col_vec[is.na(col_vec)] <- pal[["not significant"]]
    }

    p <- plotly::plot_ly(
      dd,
      x = ~log2ratio, y = ~nlog10p,
      type = "scatter", mode = "markers",
      text = hov,
      hovertemplate = "%{text}<extra></extra>",
      marker = list(size = 7, opacity = 0.85, color = col_vec)
    )

    # thresholds
    p <- plotly::add_segments(p, x = -Inf, xend = Inf, y = -log10(p_cut), yend = -log10(p_cut),
                              inherit = FALSE, showlegend = FALSE)
    p <- plotly::add_segments(p, x =  fc_cut, xend =  fc_cut, y = -Inf, yend = Inf,
                              inherit = FALSE, showlegend = FALSE)
    p <- plotly::add_segments(p, x = -fc_cut, xend = -fc_cut, y = -Inf, yend = Inf,
                              inherit = FALSE, showlegend = FALSE)
    p <- plotly::add_segments(p, x = 0, xend = 0, y = -Inf, yend = Inf,
                              inherit = FALSE, showlegend = FALSE)

    ttl <- paste0("rank: ", unique(dd$rank))
    if (!is.null(title_suffix)) ttl <- paste0(ttl, " • ", title_suffix)

    leg <- list(orientation = "h", y = -0.15)
    if (!is.null(color_var) && !is.null(color_label)) {
      leg$title <- list(text = color_label)
    }

    plotly::layout(
      p,
      title = ttl,
      xaxis = list(title = "log\u2082 ratio (group1 / group2)", zeroline = FALSE),
      yaxis = list(title = "-log\u2081\u2080(p)", zeroline = FALSE),
      legend = leg,
      showlegend = !is.null(color_var)
    )
  }

  build_one_gg <- function(dd, title_suffix = NULL) {
    # drop non-finite to avoid "Ignoring N observations" warnings
    dd <- dd[is.finite(dd$log2ratio) & is.finite(dd$nlog10p), , drop = FALSE]

    aes_base <- ggplot2::aes(x = log2ratio, y = nlog10p)
    if (!is.null(color_var)) {
      p <- ggplot2::ggplot(dd, aes_base + ggplot2::aes(color = .data[[color_var]])) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::labs(color = color_label)
    } else {
      p <- ggplot2::ggplot(dd, aes_base + ggplot2::aes(color = category_use)) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::scale_color_manual(values = significant_colors, name = NULL)
    }
    p +
      ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "gray50") +
      ggplot2::geom_vline(xintercept = c(-fc_cut, 0, fc_cut), linetype = "dashed", color = "gray50") +
      ggplot2::labs(
        x = "log\u2082 ratio (group1 / group2)",
        y = "-log\u2081\u2080(p)",
        title = paste0("rank: ", unique(dd$rank), if (!is.null(title_suffix)) paste0(" • ", title_suffix) else "")
      ) +
      theme_phip(base_size = 12) +
      ggplot2::theme(legend.position = if (!is.null(color_var)) "bottom" else "none")
  }

  # ---------- (8) compose output per rank ------------------------------------
  out <- vector("list", length(by_rank)); names(out) <- names(by_rank)
  for (rk in names(by_rank)) {
    dr <- by_rank[[rk]]
    if (!nrow(dr)) { out[[rk]] <- NULL; next }
    pieces <- split(dr, dr$pair_label, drop = TRUE)

    if (isTRUE(facet_pairs)) {
      plots <- lapply(
        names(pieces),
        function(lbl) if (interactive) build_one_plotly(pieces[[lbl]], lbl) else build_one_gg(pieces[[lbl]], lbl)
      )
      n <- length(plots)
      r <- facet_nrow; c <- facet_ncol
      if (is.null(r) && is.null(c)) { r <- floor(sqrt(n)); c <- ceiling(n / max(1, r)) }
      else if (is.null(r))         { r <- ceiling(n / max(1, c)) }
      else if (is.null(c))         { c <- ceiling(n / max(1, r)) }

      out[[rk]] <- if (interactive) {
        plotly::subplot(plots, nrows = r, shareX = TRUE, shareY = TRUE, titleX = TRUE, titleY = TRUE, margin = 0.03)
      } else {
        plots  # return list of ggplots; caller can arrange
      }
    } else {
      lst <- lapply(
        names(pieces),
        function(lbl) if (interactive) build_one_plotly(pieces[[lbl]], lbl) else build_one_gg(pieces[[lbl]], lbl)
      )
      names(lst) <- names(pieces)
      out[[rk]] <- lst
    }
  }

  out
}
