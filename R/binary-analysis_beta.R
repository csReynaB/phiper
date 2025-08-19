# generic
#' @export
compute_beta_diversity <- function(x, ...) UseMethod("compute_beta_diversity")

#' Compute beta diversity (PCoA, PERMANOVA, dispersion) for one or many ranks
#' from a phip_data object. Supports longitudinal time and multicore distances.
#' @rdname compute_beta_diversity
#' @export
#' Compute beta diversity (PCoA, PERMANOVA, dispersion) — categorical time only
#' @export
compute_beta_diversity.phip_data <- function(x,
                                             group_col,
                                             ranks = "peptide",
                                             present_by = c("exist", "fold_change"),
                                             fc_threshold = 0,
                                             norm_method = c("auto", "relative", "hellinger", "log", "none"),
                                             method = "bray",
                                             permutations = 999,
                                             sample_meta = NULL,
                                             filter_rank = NULL,
                                             # categorical longitudinal
                                             time_col = NULL,
                                             # carry extra columns from data_long into the output (joined by sample_id)
                                             carry_cols = NULL,
                                             # speed
                                             engine_dist = c("auto", "vegan", "parallelDist"),
                                             n_threads = max(1L, (if (rlang::is_installed("parallel")) parallel::detectCores() else 1L) - 1L)) {
  stopifnot(inherits(x, "phip_data"))
  present_by <- match.arg(present_by)
  norm_method <- match.arg(norm_method)
  engine_dist <- match.arg(engine_dist)

  .d <- rlang::.data
  sid <- rlang::sym("sample_id")
  pid <- rlang::sym("peptide_id")
  grp <- rlang::sym(group_col)

  tbl <- x$data_long

  # ----- ensure group column exists (join in from sample_meta if needed) -----
  if (!group_col %in% colnames(tbl)) {
    stopifnot(!is.null(sample_meta), all(c("sample_id", group_col) %in% colnames(sample_meta)))
    tbl <- tbl %>%
      dplyr::left_join(
        sample_meta %>% dplyr::select(sample_id, !!grp) %>% dplyr::distinct(),
        by = "sample_id"
      )
  }

  # ----- time handling: categorical only; DON'T factor() on DB backends -----
  has_time <- !is.null(time_col)
  tlev <- NULL
  if (has_time) {
    if (!time_col %in% colnames(tbl)) stop("time_col '", time_col, "' not found in data_long.")
    if (inherits(tbl, "tbl_lazy")) {
      # On DB backends we cannot reliably test type; collect levels as character
      tlev <- tbl %>%
        dplyr::distinct(!!rlang::sym(time_col)) %>%
        dplyr::arrange(!!rlang::sym(time_col)) %>%
        dplyr::collect() %>%
        dplyr::pull(1) %>%
        as.character()
      # no mutate(factor()) on DB; we will factor AFTER collect on small frames
    } else {
      # Local data frame: warn if not factor/character and treat as categorical
      if (!(is.factor(tbl[[time_col]]) || is.character(tbl[[time_col]]))) {
        warning(
          "Continuous time is disabled. Treating `", time_col,
          "` as categorical (coercing to character for analysis)."
        )
      }
      tlev <- as.character(sort(unique(tbl[[time_col]])))
    }
  }

  # ----- presence rule (stay lazy) -----
  if (present_by == "exist") {
    stopifnot("exist" %in% colnames(tbl))
    tbl <- tbl %>% dplyr::filter(!!rlang::sym("exist") == 1L)
    value_expr <- 1L
  } else {
    stopifnot("fold_change" %in% colnames(tbl))
    tbl <- tbl %>% dplyr::filter(!!rlang::sym("fold_change") > !!fc_threshold)
    value_expr <- !!rlang::sym("fold_change")
  }

  # ----- normalization helper -----
  .normalize <- function(m, how) {
    if (how == "auto") how <- "relative"
    if (how == "none") {
      return(m)
    }
    rs <- rowSums(m, na.rm = TRUE)
    rs[rs == 0] <- 1
    rel <- m / rs
    switch(how,
      relative  = rel,
      hellinger = sqrt(rel),
      log       = log1p(m),
      rel
    )
  }

  # ----- distance engine (parallelDist -> vegan fallback) -----
  .dist_fast <- function(m, method, engine, n_threads) {
    if (engine %in% c("auto", "parallelDist") && rlang::is_installed("parallelDist")) {
      out <- try(parallelDist::parDist(m, method = "bray", threads = n_threads), silent = TRUE)
      if (!inherits(out, "try-error")) {
        return(out)
      }
    }
    vegan::vegdist(m, method = method)
  }

  main_con <- dbplyr::remote_con(x$data_long)

  build_for_rank <- function(rank) {
    cat("# 1) lazy -> collect agg")
    t0 <- proc.time()
    rk <- tolower(rank)
    if (rk %in% c("peptide", "peptide_id", "peptideid")) {
      ranked <- tbl %>%
        dplyr::transmute(!!sid, !!grp, rank_val = !!pid, value = !!value_expr)
      nm <- if (norm_method == "auto" && present_by == "exist") "none" else norm_method
    } else {
      stopifnot(rank %in% colnames(x$peptide_library))
      map_df <- x$peptide_library %>%
        dplyr::select(!!pid, !!rlang::sym(rank)) %>%
        dplyr::rename(rank_val = !!rlang::sym(rank)) %>%
        dplyr::collect()
      map_tbl <- dplyr::copy_to(
        main_con, map_df,
        name = paste0("peplib_map_tmp_", as.integer(Sys.time()), "_", rank),
        temporary = TRUE, overwrite = TRUE
      )
      ranked <- tbl %>%
        dplyr::inner_join(map_tbl, by = "peptide_id") %>%
        dplyr::transmute(!!sid, !!grp, rank_val, value = !!value_expr)
      nm <- if (norm_method == "auto") "relative" else norm_method
    }

    # Filter rank (vector or function) while still lazy where possible
    if (is.function(filter_rank)) {
      # get allowed set in R, then filter on DB
      allowed <- ranked %>%
        dplyr::distinct(.data$rank_val) %>%
        dplyr::collect() %>%
        dplyr::pull()
      allowed <- allowed[as.logical(filter_rank(allowed))]
      ranked <- ranked %>% dplyr::filter(.data$rank_val %in% !!allowed)
    } else if (!is.null(filter_rank)) {
      ranked <- ranked %>% dplyr::filter(.data$rank_val %in% !!filter_rank)
    }

    # ---- aggregate on DB, then collect compact table ----
    agg <- ranked %>%
      dplyr::group_by(!!sid, !!grp, .data$rank_val) %>%
      dplyr::summarise(abund = sum(.data$value), .groups = "drop") %>%
      dplyr::collect()

    if (!nrow(agg)) {
      return(NULL)
    }

    # ---- row metadata (small): group + time + carry_cols ----
    meta_cols <- unique(c("sample_id", group_col, if (has_time) time_col, carry_cols))
    sm <- x$data_long %>%
      dplyr::distinct(dplyr::across(dplyr::all_of(intersect(meta_cols, colnames(x$data_long))))) %>%
      dplyr::collect() %>%
      dplyr::rename(group = !!grp)

    # enforce character ids; categorical time AFTER collect
    sm$sample_id <- as.character(sm$sample_id)
    if (has_time && time_col %in% names(sm)) {
      sm[[time_col]] <- factor(as.character(sm[[time_col]]), levels = tlev)
      sm$time_fac <- sm[[time_col]]
    }
    if (!rlang::is_installed("data.table")) stop("Install data.table for this path.")
    #' @importFrom data.table as.data.table dcast set setDT
    DT <- data.table::as.data.table(agg)
    # ensure character ids
    DT[, sample_id := as.character(sample_id)]
    # wide table (fill=0), first column becomes rownames
    wide <- data.table::dcast(DT, sample_id ~ rank_val, value.var = "abund", fill = 0L)
    rn <- wide[["sample_id"]]
    mat <- as.matrix(wide[, -"sample_id"])
    rownames(mat) <- rn

    # align sm
    sm <- sm[match(rownames(mat), sm$sample_id), , drop = FALSE]
    rownames(sm) <- sm$sample_id
    sm_df <- as.data.frame(sm)

    # ---- normalize & distances ----
    mat_norm <- .normalize(mat, nm)
    dist_obj <- .dist_fast(mat_norm, method = method, engine = engine_dist, n_threads = n_threads)

    # ---- PCoA ----
    pcoa <- stats::cmdscale(dist_obj, eig = TRUE, k = 2)
    eig <- pcoa$eig
    ve <- if (is.null(eig) || sum(eig) == 0) c(NA_real_, NA_real_) else round(100 * eig[1:2] / sum(eig), 2)
    print(proc.time() - t0)
    cat("\n\n")
    cat("# 5) stats: adonis2, betadisper, permutest")
    t0 <- proc.time()
    # ---- tests (categorical only) ----
    run_perm <- isTRUE(permutations > 0) && length(unique(sm$group)) > 1 && nrow(sm) >= 4
    if (run_perm && has_time && "time_fac" %in% names(sm)) {
      perm <- try(vegan::adonis2(dist_obj ~ group + time_fac + group:time_fac,
        data = sm_df, permutations = permutations
      ), silent = TRUE)
      disp_vec <- interaction(sm$group, sm$time_fac, drop = TRUE)
    } else if (run_perm) {
      perm <- try(vegan::adonis2(dist_obj ~ group, data = sm_df, permutations = permutations), silent = TRUE)
      disp_vec <- sm$group
    } else {
      perm <- NULL
      disp_vec <- sm$group
    }

    disp <- try(vegan::betadisper(dist_obj, disp_vec), silent = TRUE)
    p_disp <- if (inherits(disp, "try-error")) {
      NA_real_
    } else {
      out <- try(vegan::permutest(disp, permutations = permutations), silent = TRUE)
      if (inherits(out, "try-error")) NA_real_ else out$tab$`Pr(>F)`[1]
    }
    print(proc.time() - t0)
    # ---- tidy outputs ----
    scores <- tibble::tibble(
      sample_id = rownames(pcoa$points),
      pcoa1     = pcoa$points[, 1],
      pcoa2     = pcoa$points[, 2]
    ) %>%
      dplyr::left_join(sm, by = "sample_id") %>%
      dplyr::mutate(rank = rank, .before = 1)

    # dispersion distances table
    disp_df <- if (!is.null(disp) && !inherits(disp, "try-error")) {
      tibble::tibble(
        sample_id = names(disp$distances),
        distance = as.numeric(disp$distances),
        group_disp = as.character(disp$group),
        rank = rank
      ) %>%
        dplyr::left_join(sm %>% dplyr::select(sample_id, group, dplyr::any_of("time_fac")), by = "sample_id")
    } else {
      tibble::tibble(sample_id = character(), distance = double(), group_disp = character(), rank = character())
    }

    perm_tbl <- if (is.null(perm) || inherits(perm, "try-error")) {
      tibble::tibble(rank = rank, term = character(), p_value = numeric())
    } else {
      tibble::tibble(rank = rank, term = rownames(perm), p_value = perm$`Pr(>F)`)
    }

    list(
      rank = rank,
      scores = scores,
      var_exp = tibble::tibble(rank = rank, pcoa1 = ve[1], pcoa2 = ve[2]),
      tests = perm_tbl,
      dispersion = disp_df
    )
  }

  outs <- lapply(unique(ranks), build_for_rank)
  outs <- outs[!vapply(outs, is.null, logical(1))]
  if (!length(outs)) stop("No data available to compute beta diversity.")

  out <- list(
    pcoa          = dplyr::bind_rows(lapply(outs, `[[`, "scores")),
    var_explained = dplyr::bind_rows(lapply(outs, `[[`, "var_exp")),
    tests         = dplyr::bind_rows(lapply(outs, `[[`, "tests")),
    dispersion    = dplyr::bind_rows(lapply(outs, `[[`, "dispersion"))
  )
  class(out) <- c("phip_beta", class(out))
  out
}


#' Faceted PCoA (group base color; per-time light→dark shades; solid ellipses)
#' @export
plot_beta_pcoa <- function(beta,
                           group_col = "group",
                           custom_colors = NULL, # named by group
                           ellipse = TRUE,
                           show_legend = TRUE, # legend for Group | Time
                           facet_scales = "fixed",
                           ncol = 2,
                           time_col = NULL, # if NULL, tries 'time_fac'
                           centroid_size = 6,
                           centroid_stroke = 1.1,
                           draw_centroid_paths = TRUE,
                           lighten_from = 0.85) { # stronger lightening for clarity
  stopifnot(inherits(beta, "phip_beta"))
  df <- beta$pcoa
  ve <- beta$var_explained

  # facet strip text with explained variance
  df_lab <- dplyr::left_join(
    ve, dplyr::tibble(rank = unique(df$rank)),
    by = "rank"
  ) %>%
    dplyr::mutate(strip = sprintf("%s  (PCoA 1 %s%%, PCoA 2 %s%%)", rank, pcoa1, pcoa2))
  df <- dplyr::left_join(df, df_lab[, c("rank", "strip")], by = "rank")

  gsym <- rlang::sym(group_col)

  # resolve time (categorical only)
  if (is.null(time_col) && "time_fac" %in% names(df)) time_col <- "time_fac"
  has_time <- !is.null(time_col) && time_col %in% names(df)
  if (has_time && !(is.factor(df[[time_col]]) || is.character(df[[time_col]]))) {
    stop("plot_beta_pcoa: `time_col` must be a factor/character (continuous time disabled).")
  }
  if (has_time) df[[time_col]] <- factor(df[[time_col]])
  tsym <- if (has_time) rlang::sym(time_col) else NULL

  # base colors per group
  groups <- sort(unique(df[[group_col]]))
  if (is.null(custom_colors)) {
    base_cols <- scales::hue_pal()(length(groups))
    names(base_cols) <- groups
  } else {
    if (is.null(names(custom_colors))) {
      stopifnot(length(custom_colors) >= length(groups))
      base_cols <- custom_colors[seq_along(groups)]
      names(base_cols) <- groups
    } else {
      base_cols <- custom_colors[groups]
    }
  }

  # helper to lighten
  .make_shades <- function(base, n, from = 0.85) {
    if (n <= 1) {
      return(rep(base, n))
    }
    if (rlang::is_installed("colorspace")) {
      vapply(
        seq(from, 0, length.out = n),
        function(a) colorspace::lighten(base, amount = a), character(1)
      )
    } else {
      grDevices::colorRampPalette(c("#FFFFFF", base))(n)
    }
  }

  # --- build a palette only for PRESENT combos (per group) ---
  if (has_time) {
    # make sure time has an order we honor inside each group
    tlev_all <- levels(df[[time_col]])
    # list of present time levels per group (respect factor order)
    present_by_group <- lapply(groups, function(g) {
      lev <- levels(droplevels(df[df[[group_col]] == g, time_col, drop = TRUE]))
      if (is.null(lev)) character(0) else lev
    })
    names(present_by_group) <- groups

    shades_list <- lapply(groups, function(g) {
      lev <- present_by_group[[g]]
      cols <- .make_shades(base_cols[[g]], length(lev), from = lighten_from)
      stats::setNames(cols, paste(g, lev, sep = " | "))
    })
    shades <- unlist(shades_list, use.names = TRUE)

    # gt factor with levels exactly those present
    df$gt <- factor(paste(df[[group_col]], df[[time_col]], sep = " | "),
      levels = names(shades)
    )
  } else {
    # no time: color by group only
    shades <- base_cols
    names(shades) <- groups
    df$gt <- factor(df[[group_col]], levels = groups)
  }

  # --- centroids ---
  if (has_time) {
    cents <- df %>%
      dplyr::group_by(rank, !!gsym, !!tsym) %>%
      dplyr::summarise(pcoa1 = mean(.data$pcoa1), pcoa2 = mean(.data$pcoa2), .groups = "drop") %>%
      dplyr::mutate(gt = factor(paste(.data[[group_col]], .data[[time_col]], sep = " | "),
        levels = levels(df$gt)
      )) %>%
      dplyr::arrange(rank, !!gsym, !!tsym)
  } else {
    cents <- df %>%
      dplyr::group_by(rank, !!gsym) %>%
      dplyr::summarise(pcoa1 = mean(.data$pcoa1), pcoa2 = mean(.data$pcoa2), .groups = "drop") %>%
      dplyr::mutate(gt = factor(.data[[group_col]], levels = levels(df$gt)))
  }

  # --- plotting ---
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = .data$pcoa1, y = .data$pcoa2, color = .data$gt),
      size = 2.8, alpha = 0.85
    )

  if (isTRUE(ellipse)) {
    p <- p + ggplot2::stat_ellipse(
      data = df,
      ggplot2::aes(
        x = .data$pcoa1, y = .data$pcoa2,
        group = interaction(.data$rank, .data$gt, drop = TRUE),
        color = .data$gt
      ),
      type = "t", level = 0.95, linewidth = 0.9, linetype = "solid", show.legend = FALSE
    )
  }

  p <- p +
    ggplot2::geom_point(
      data = cents,
      ggplot2::aes(x = .data$pcoa1, y = .data$pcoa2, fill = .data$gt),
      shape = 21, size = centroid_size, stroke = centroid_stroke, color = "black"
    )

  if (isTRUE(draw_centroid_paths)) {
    # connect within each group, in time order (if present)
    if (has_time) {
      p <- p + ggplot2::geom_path(
        data = cents,
        ggplot2::aes(
          x = .data$pcoa1, y = .data$pcoa2,
          group = interaction(.data$rank, .data[[group_col]], drop = TRUE),
          color = .data$gt
        ),
        linewidth = 0.9, alpha = 0.9
      )
    } else {
      p <- p + ggplot2::geom_path(
        data = cents,
        ggplot2::aes(
          x = .data$pcoa1, y = .data$pcoa2,
          group = interaction(.data$rank, .data[[group_col]], drop = TRUE),
          color = .data$gt
        ),
        linewidth = 0.9, alpha = 0.9
      )
    }
  }

  # manual palette tied to gt levels (works in all layers)
  p <- p +
    ggplot2::scale_color_manual(
      values = shades, name = "Group | Time",
      drop = TRUE, guide = if (show_legend) "legend" else "none"
    ) +
    ggplot2::scale_fill_manual(values = shades, drop = TRUE, guide = "none") +
    ggplot2::labs(x = "PCoA 1", y = "PCoA 2") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = if (show_legend) "right" else "none"
    ) +
    ggplot2::facet_wrap(~strip, ncol = ncol, scales = facet_scales)

  p
}

#' Faceted beta-dispersion (distance-to-centroid) boxplots
#' Uses group base colors, shaded by time within each group (light → dark).
#' @export
plot_beta_dispersion <- function(beta,
                                 group_col = "group",
                                 custom_colors = NULL, # named by group
                                 facet_scales = "free_y",
                                 ncol = 2,
                                 time_col = NULL, # if NULL, tries 'time_fac'
                                 show_legend = TRUE,
                                 lighten_from = 0.85) { # stronger lightening for clarity
  stopifnot(inherits(beta, "phip_beta"))
  df <- beta$dispersion
  if (!nrow(df)) {
    return(ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "No dispersion available"))
  }

  # prefer the grouping actually used by betadisper if present (e.g., group×time)
  gcol <- if ("group_disp" %in% names(df)) "group_disp" else group_col

  # resolve time (categorical only)
  if (is.null(time_col) && "time_fac" %in% names(df)) time_col <- "time_fac"
  has_time <- !is.null(time_col) && time_col %in% names(df)
  if (has_time && !(is.factor(df[[time_col]]) || is.character(df[[time_col]]))) {
    stop("plot_beta_dispersion: `time_col` must be factor/character (continuous time disabled).")
  }
  if (has_time) df[[time_col]] <- factor(df[[time_col]])

  # base colors per group
  groups <- sort(unique(df[[group_col]]))
  if (is.null(custom_colors)) {
    base_cols <- scales::hue_pal()(length(groups))
    names(base_cols) <- groups
  } else {
    if (is.null(names(custom_colors))) {
      stopifnot(length(custom_colors) >= length(groups))
      base_cols <- custom_colors[seq_along(groups)]
      names(base_cols) <- groups
    } else {
      base_cols <- custom_colors[groups]
    }
  }

  # lightening helper
  .make_shades <- function(base, n, from = 0.85) {
    if (n <= 1) {
      return(rep(base, n))
    }
    if (rlang::is_installed("colorspace")) {
      vapply(
        seq(from, 0, length.out = n),
        function(a) colorspace::lighten(base, amount = a), character(1)
      )
    } else {
      grDevices::colorRampPalette(c("#FFFFFF", base))(n)
    }
  }

  # Build the palette for present combos and a gt factor for mapping
  if (has_time) {
    # present time levels per group (respect factor order)
    tlev_all <- levels(df[[time_col]])
    present_by_group <- lapply(groups, function(g) {
      lev <- levels(droplevels(df[df[[group_col]] == g, time_col, drop = TRUE]))
      if (is.null(lev)) character(0) else lev
    })
    names(present_by_group) <- groups

    shades_list <- lapply(groups, function(g) {
      lev <- present_by_group[[g]]
      cols <- .make_shades(base_cols[[g]], length(lev), from = lighten_from)
      stats::setNames(cols, paste(g, lev, sep = " | "))
    })
    shades <- unlist(shades_list, use.names = TRUE)

    df$gt <- factor(paste(df[[group_col]], df[[time_col]], sep = " | "),
      levels = names(shades)
    )
    x_lab <- "Group | Time"
  } else {
    shades <- base_cols
    names(shades) <- groups
    df$gt <- factor(df[[group_col]], levels = groups)
    x_lab <- "Group"
  }

  # Boxplot/jitter mapped to gt with manual palette (works across layers)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$gt, y = .data$distance, fill = .data$gt)) +
    ggplot2::geom_boxplot(outlier.shape = NA, show.legend = show_legend) +
    ggplot2::geom_jitter(color = "black", size = 1, width = 0.2, alpha = 0.3, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = shades, name = x_lab, drop = TRUE) +
    ggplot2::labs(x = x_lab, y = "Distance to centroid") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.6, hjust = 0.5)
    ) +
    ggplot2::facet_wrap(~rank, ncol = ncol, scales = facet_scales)

  p
}
