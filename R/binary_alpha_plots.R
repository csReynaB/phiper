# ------------------------------------------------------------------------------
# GENERICS
# ------------------------------------------------------------------------------
plot_enrichment_counts <- function(phip_data,
                                   group_cols = NULL,
                                   prevalence_threshold = 0.05,
                                   custom_colors = NULL,
                                   binwidth = 1,
                                   group_interaction = FALSE,
                                   interaction_sep = " * ",
                                   ...) {
  UseMethod("plot_enrichment_counts")
}

plot_alpha_diversity <- function(x, ...) {
  UseMethod("plot_alpha_diversity")
}

# ------------------------------------------------------------------------------
# METHODS (S3)
# ------------------------------------------------------------------------------

# -- plot_enrichment_counts: <phip_data> ---------------------------------------
plot_enrichment_counts.phip_data <- function(phip_data,
                                             group_cols = NULL,
                                             prevalence_threshold = 0.05,
                                             custom_colors = NULL,
                                             binwidth = 1,
                                             group_interaction = FALSE,
                                             interaction_sep = " * ",
                                             ...) {
  x <- phip_data
  stopifnot(inherits(x, "phip_data"))
  .data <- rlang::.data

  .ph_with_timing(
    headline = "Plotting enrichment counts (<phip_data>)",
    step = if (is.null(group_cols)) {
      "group_cols: <none>"
    } else {
      sprintf("group_cols: %s", paste(add_quotes(group_cols, 1L),
        collapse = ", "
      ))
    },
    expr = {
      tbl <- x$data_long

      if (isTRUE(x$meta$full_cross) && ("exist" %in% colnames(tbl))) {
        red_txt <- tryCatch(
          {
            ep <- as.numeric(x$meta$exist_prop)
            if (is.finite(ep) && ep > 0) sprintf("~%.1fx", 1 / ep) else "<unknown>"
          },
          error = function(e) "<unknown>"
        )
        .ph_log_info(
          "Full-cross detected; pruning non-existing rows before plotting",
          bullets = c(
            "rule: keep exist == 1",
            sprintf("estimated reduction: %s", red_txt)
          )
        )
        tbl <- dplyr::filter(tbl, .data$exist == 1L)
      }

      if (!is.null(group_cols) && length(group_cols)) {
        missing_gcs <- setdiff(group_cols, colnames(tbl))
        if (length(missing_gcs)) {
          .ph_abort(
            headline = "Missing grouping columns in data_long.",
            step = "input validation",
            bullets = sprintf(
              "missing: %s",
              paste(add_quotes(missing_gcs, 1L),
                collapse = ", "
              )
            )
          )
        }
      }

      if (is.null(group_cols) || !length(group_cols)) {
        .plot_enrichment_counts_one(
          tbl,
          group_col = NULL,
          prevalence_threshold = prevalence_threshold,
          custom_colors = custom_colors,
          binwidth = binwidth
        )
      } else {
        plots <- lapply(group_cols, function(gc) {
          .plot_enrichment_counts_one(
            tbl,
            group_col = gc,
            prevalence_threshold = prevalence_threshold,
            custom_colors = custom_colors,
            binwidth = binwidth,
            title_label = gc
          )
        })
        names(plots) <- group_cols

        if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
          combo_nm <- paste(group_cols, collapse = interaction_sep)
          inter_col <- "..phip_interaction.."
          tbl_inter <- dplyr::mutate(tbl, !!rlang::sym(inter_col) :=
            paste(!!!rlang::syms(group_cols),
              sep = interaction_sep
            ))

          plots[[combo_nm]] <- .plot_enrichment_counts_one(
            tbl_inter,
            group_col = inter_col,
            prevalence_threshold = prevalence_threshold,
            custom_colors = custom_colors,
            binwidth = binwidth,
            title_label = combo_nm
          )
        } else if (isTRUE(group_interaction) && length(group_cols) < 2L) {
          .ph_warn(
            headline = "group_interaction requested but fewer than
            2 group_cols supplied.",
            step     = "interaction plot",
            bullets  = "Interaction plot skipped."
          )
        }

        if (length(plots) == 1L) plots[[1L]] else plots
      }
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# -- plot_enrichment_counts: data.frame ----------------------------------------
plot_enrichment_counts.data.frame <- function(phip_data,
                                              group_cols = NULL,
                                              prevalence_threshold = 0.05,
                                              custom_colors = NULL,
                                              binwidth = 1,
                                              group_interaction = FALSE,
                                              interaction_sep = " * ",
                                              ...) {
  tbl <- phip_data
  .data <- rlang::.data

  .ph_with_timing(
    headline = "Plotting enrichment counts (data.frame)",
    step = if (is.null(group_cols)) {
      "group_cols: <none>"
    } else {
      sprintf("group_cols: %s", paste(add_quotes(group_cols, 1L),
        collapse = ", "
      ))
    },
    expr = {
      need_base <- c("sample_id", "peptide_id", "exist")
      miss_base <- setdiff(need_base, colnames(tbl))
      if (length(miss_base)) {
        .ph_abort(
          headline = "Missing required columns in data.frame.",
          step = "input validation",
          bullets = sprintf(
            "missing: %s",
            paste(add_quotes(miss_base, 1L),
              collapse = ", "
            )
          )
        )
      }

      if (is.null(group_cols) || !length(group_cols)) {
        .plot_enrichment_counts_one(
          tbl,
          group_col = NULL,
          prevalence_threshold = prevalence_threshold,
          custom_colors = custom_colors,
          binwidth = binwidth
        )
      } else {
        plots <- lapply(group_cols, function(gc) {
          if (!gc %in% colnames(tbl)) {
            .ph_abort(
              headline = "Grouping column not found in data.frame.",
              step     = "input validation",
              bullets  = sprintf("group_col: %s", add_quotes(gc, 1L))
            )
          }
          .plot_enrichment_counts_one(
            tbl,
            group_col = gc,
            prevalence_threshold = prevalence_threshold,
            custom_colors = custom_colors,
            binwidth = binwidth,
            title_label = gc
          )
        })
        names(plots) <- group_cols

        if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
          combo_nm <- paste(group_cols, collapse = interaction_sep)
          inter_col <- "..phip_interaction.."
          tbl_inter <- dplyr::mutate(tbl, !!rlang::sym(inter_col) :=
            paste(!!!rlang::syms(group_cols),
              sep = interaction_sep
            ))

          plots[[combo_nm]] <- .plot_enrichment_counts_one(
            tbl_inter,
            group_col = inter_col,
            prevalence_threshold = prevalence_threshold,
            custom_colors = custom_colors,
            binwidth = binwidth,
            title_label = combo_nm
          )
        }

        if (length(plots) == 1L) plots[[1L]] else plots
      }
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# -- plot_alpha_diversity: <phip_data> -----------------------------------------
plot_alpha_diversity.phip_data <- function(
  x,
  metric = c("richness", "shannon_diversity", "simpson_diversity"),
  group_col = "group",
  rank_col = "rank",
  filter_groups = NULL,
  filter_ranks = NULL,
  custom_colors = NULL,
  facet_by_rank = TRUE,
  ncol = 2,
  facet_scales = "fixed",
  # longitudinal
  time_col = NULL,
  continuous_mode = c("gam", "binned", "loess"),
  gam_k = 7,
  point_alpha = 0.25,
  # CIs
  ci_method = c("model", "bootstrap"),
  ci_level = 0.95,
  boot_R = 500,
  boot_seed = NULL,
  boot_progress = TRUE,
  ci_fill = "grey70",
  ci_alpha = 0.15,
  ...
) {
  stopifnot(inherits(x, "phip_data"))
  continuous_mode <- match.arg(continuous_mode)
  ci_method <- match.arg(ci_method)

  .ph_with_timing(
    headline = "Plotting alpha diversity (<phip_data>)",
    step = sprintf(
      "metrics: %s",
      paste(add_quotes(tolower(metric), 1L), collapse = ", ")
    ),
    expr = {
      if (!is.null(time_col) && isTRUE(!isTRUE(x$meta$longitudinal))) {
        .ph_abort(
          headline = "Longitudinal plotting not allowed for this object.",
          step     = "meta$longitudinal check",
          bullets  = "x$meta$longitudinal is FALSE"
        )
      }

      alpha_list <- compute_alpha_diversity(
        x,
        group_col  = group_col,
        ranks      = rank_col,
        carry_cols = time_col
      )

      plot_alpha_diversity(
        alpha_list,
        metric = metric,
        group_col = group_col,
        rank_col = rank_col,
        filter_groups = filter_groups,
        filter_ranks = filter_ranks,
        custom_colors = custom_colors,
        facet_by_rank = facet_by_rank,
        ncol = ncol,
        facet_scales = facet_scales,
        time_col = time_col,
        continuous_mode = continuous_mode,
        gam_k = gam_k,
        point_alpha = point_alpha,
        ci_method = ci_method,
        ci_level = ci_level,
        boot_R = boot_R,
        boot_seed = boot_seed,
        boot_progress = boot_progress,
        ci_fill = ci_fill,
        ci_alpha = ci_alpha
      )
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# -- plot_alpha_diversity: data.frame ------------------------------------------
plot_alpha_diversity.data.frame <- function(
  x,
  metric = c("richness", "shannon_diversity", "simpson_diversity"),
  group_col = "group",
  rank_col = "rank",
  filter_groups = NULL,
  filter_ranks = NULL,
  custom_colors = NULL,
  facet_by_rank = TRUE,
  ncol = 2,
  facet_scales = "fixed",
  # longitudinal
  time_col = NULL,
  continuous_mode = c("gam", "binned", "loess"),
  gam_k = 7,
  point_alpha = 0.25,
  # CIs
  ci_method = c("model", "bootstrap"),
  ci_level = 0.95,
  boot_R = 500,
  boot_seed = NULL,
  boot_progress = TRUE,
  ci_fill = "grey70",
  ci_alpha = 0.15,
  ...
) {
  continuous_mode <- match.arg(continuous_mode)
  ci_method <- match.arg(ci_method)

  .ph_with_timing(
    headline = "Plotting alpha diversity (data.frame)",
    step = sprintf(
      "metrics: %s",
      paste(add_quotes(tolower(metric), 1L), collapse = ", ")
    ),
    expr = {
      alpha_list <- compute_alpha_diversity(
        x,
        group_col  = group_col,
        ranks      = rank_col,
        carry_cols = time_col
      )

      plot_alpha_diversity(
        alpha_list,
        metric = metric,
        group_col = group_col,
        rank_col = rank_col,
        filter_groups = filter_groups,
        filter_ranks = filter_ranks,
        custom_colors = custom_colors,
        facet_by_rank = facet_by_rank,
        ncol = ncol,
        facet_scales = facet_scales,
        time_col = time_col,
        continuous_mode = continuous_mode,
        gam_k = gam_k,
        point_alpha = point_alpha,
        ci_method = ci_method,
        ci_level = ci_level,
        boot_R = boot_R,
        boot_seed = boot_seed,
        boot_progress = boot_progress,
        ci_fill = ci_fill,
        ci_alpha = ci_alpha
      )
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# -- plot_alpha_diversity: precomputed -----------------------------------------
plot_alpha_diversity.phip_alpha_diversity <- function(
  x,
  metric = c("richness", "shannon_diversity", "simpson_diversity"),
  group_col = "group",
  rank_col = "rank",
  filter_groups = NULL,
  filter_ranks = NULL,
  custom_colors = NULL,
  facet_by_rank = TRUE,
  ncol = 2,
  facet_scales = "fixed",
  # longitudinal
  time_col = NULL,
  continuous_mode = c("gam", "binned", "loess"),
  gam_k = 7,
  point_alpha = 0.25,
  # CIs
  ci_method = c("model", "bootstrap"),
  ci_level = 0.95,
  boot_R = 500,
  boot_seed = NULL,
  boot_progress = TRUE,
  ci_fill = "grey70",
  ci_alpha = 0.15,
  ...
) {
  continuous_mode <- match.arg(continuous_mode)
  ci_method <- match.arg(ci_method)

  .ph_with_timing(
    headline = "Plotting alpha diversity (precomputed)",
    step = sprintf(
      "metrics: %s",
      paste(add_quotes(tolower(metric), 1L), collapse = ", ")
    ),
    expr = {
      alpha_df <- if (is.list(x) && !is.null(group_col) &&
        group_col %in% names(x)) {
        tibble::as_tibble(x[[group_col]])
      } else if (is.list(x)) {
        inner <- x[vapply(
          x, function(el) inherits(el, "data.frame"),
          logical(1)
        )]
        if (length(inner) == 1L) {
          tibble::as_tibble(inner[[1L]])
        } else if (length(inner) >= 1L) {
          suppressWarnings(dplyr::bind_rows(inner))
        } else {
          if (length(x) && inherits(x[[1L]], "data.frame")) {
            tibble::as_tibble(x[[1L]])
          } else {
            tibble::as_tibble(x)
          }
        }
      } else {
        tibble::as_tibble(x)
      }

      if (!is.null(attr(x, "group_cols", exact = TRUE))) {
        attr(alpha_df, "group_cols") <- attr(x, "group_cols", exact = TRUE)
      }
      if (!is.null(attr(x, "ranks", exact = TRUE))) {
        attr(alpha_df, "ranks") <- attr(x, "ranks", exact = TRUE)
      }

      metrics_norm <- tolower(metric)
      for (m in metrics_norm) {
        if (!m %in% names(alpha_df)) {
          .ph_abort(
            headline = "Metric column not found in alpha data.",
            step     = "input validation",
            bullets  = sprintf("missing: %s", add_quotes(m, 1L))
          )
        }
      }

      if (!is.null(group_col) && !group_col %in% names(alpha_df)) {
        .ph_abort(
          headline = "Grouping column not found in precomputed alpha data.",
          step = "input validation",
          bullets = c(
            sprintf("requested group_col: %s", add_quotes(group_col, 1L)),
            sprintf(
              "available columns: %s",
              paste(add_quotes(colnames(alpha_df), 1L), collapse = ", ")
            )
          )
        )
      }

      if (!is.null(rank_col) && !rank_col %in% names(alpha_df)) {
        .ph_warn(
          headline = "Rank column not found in precomputed alpha data;
          disabling faceting.",
          step     = "input validation",
          bullets  = sprintf("rank_col: %s", add_quotes(rank_col, 1L))
        )
        rank_col <- NULL
      }

      if (!is.null(group_col) && !is.null(filter_groups)) {
        alpha_df <- dplyr::filter(
          alpha_df,
          .data[[group_col]] %in% !!filter_groups
        )
      }
      if (!is.null(rank_col) && !is.null(filter_ranks)) {
        alpha_df <- dplyr::filter(
          alpha_df,
          .data[[rank_col]] %in% !!filter_ranks
        )
      }

      plots <- lapply(metrics_norm, function(m) {
        .build_alpha_plot(
          alpha_df = alpha_df,
          metric = m,
          group_col = group_col,
          rank_col = rank_col,
          filter_groups = filter_groups,
          filter_ranks = filter_ranks,
          custom_colors = custom_colors,
          facet_by_rank = facet_by_rank,
          ncol = ncol,
          facet_scales = facet_scales,
          time_col = if (!is.null(time_col) &&
            time_col %in% names(alpha_df)) {
            time_col
          } else {
            NULL
          },
          continuous_mode = continuous_mode,
          gam_k = gam_k,
          point_alpha = point_alpha,
          ci_method = ci_method,
          ci_level = ci_level,
          boot_R = boot_R,
          boot_seed = boot_seed,
          boot_progress = boot_progress,
          ci_fill = ci_fill,
          ci_alpha = ci_alpha
        )
      })
      names(plots) <- metrics_norm
      if (length(plots) == 1L) plots[[1L]] else plots
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ------------------------------------------------------------------------------
# HELPERS (internal) – put at the end as requested
# ------------------------------------------------------------------------------

# -- enrichment: shared builder -------------------------------------------------
.plot_enrichment_counts_one <- function(tbl,
                                        group_col = NULL,
                                        prevalence_threshold = 0.05,
                                        custom_colors = NULL,
                                        binwidth = 1,
                                        title_label = if (is.null(group_col)) {
                                          "all samples"
                                        } else {
                                          group_col
                                        }) {
  .data <- rlang::.data
  grp_sym <- if (!is.null(group_col)) rlang::sym(group_col) else NULL

  .ph_with_timing(
    headline = "Building enrichment count plot",
    step = if (is.null(group_col)) {
      "no grouping (aggregate)"
    } else {
      sprintf("grouping variable: '%s'", group_col)
    },
    expr = {
      need <- c("sample_id", "peptide_id", "exist")
      if (!is.null(group_col)) need <- c(need, group_col)
      miss <- setdiff(need, colnames(tbl))
      if (length(miss)) {
        .ph_abort(
          headline = "Missing required columns for plotting.",
          step = "input validation",
          bullets = sprintf(
            "missing: %s",
            paste(add_quotes(miss, 1L), collapse = ", ")
          )
        )
      }

      tbl <- if (is.null(group_col)) {
        dplyr::mutate(tbl, Cohort = "All samples")
      } else {
        dplyr::mutate(tbl, Cohort = !!grp_sym)
      }

      tbl_present <- dplyr::filter(tbl, .data$exist > 0)

      cohort_sizes <- tbl |>
        dplyr::distinct(.data$sample_id, .data$Cohort) |>
        dplyr::count(.data$Cohort, name = "n_samples") |>
        dplyr::collect()

      if (nrow(cohort_sizes) == 0L || any(cohort_sizes$n_samples <= 0)) {
        .ph_warn(
          headline = "No samples for at least one cohort (facet may be empty).",
          step     = "cohort sizes"
        )
        cohort_sizes <- dplyr::filter(cohort_sizes, .data$n_samples > 0)
      }

      pep_counts <- tbl_present |>
        dplyr::group_by(.data$Cohort, .data$peptide_id) |>
        dplyr::summarise(
          n_present = dplyr::n_distinct(.data$sample_id),
          .groups = "drop"
        ) |>
        dplyr::filter(.data$n_present > 0) |>
        dplyr::collect()

      thresholds <- cohort_sizes |>
        dplyr::mutate(thresh = ceiling(.data$n_samples * prevalence_threshold))

      n_thresh_tbl <- pep_counts |>
        dplyr::inner_join(dplyr::select(thresholds, .data$Cohort, .data$thresh),
          by = "Cohort"
        ) |>
        dplyr::group_by(.data$Cohort) |>
        dplyr::summarise(
          n_peptides_thresh = sum(.data$n_present >= .data$thresh),
          .groups = "drop"
        )

      overall_tbl <- pep_counts |>
        dplyr::group_by(.data$Cohort) |>
        dplyr::summarise(
          n_overall = dplyr::n_distinct(.data$peptide_id),
          .groups = "drop"
        )

      thresholds_df <- thresholds |>
        dplyr::left_join(n_thresh_tbl, by = "Cohort") |>
        dplyr::left_join(overall_tbl, by = "Cohort") |>
        dplyr::mutate(
          n_peptides_thresh = dplyr::coalesce(.data$n_peptides_thresh, 0L),
          n_overall         = dplyr::coalesce(.data$n_overall, 0L),
          y_line            = pmax(.data$n_peptides_thresh, 1L),
          x_mid             = (.data$thresh + .data$n_samples) / 2
        )

      real_order <- thresholds_df$Cohort
      count_df <- pep_counts |>
        dplyr::mutate(Cohort = factor(.data$Cohort, levels = real_order))

      lvl <- levels(count_df$Cohort)
      pal_map <- if (is.null(custom_colors)) {
        .phip_palette_map(lvl)
      } else {
        if (is.null(names(custom_colors))) {
          stats::setNames(rep_len(custom_colors, length(lvl)), lvl)
        } else {
          stats::setNames(unname(custom_colors[lvl]), lvl)
        }
      }

      count_df <- dplyr::mutate(
        count_df,
        fill_col = unname(pal_map[as.character(.data$Cohort)])
      )

      label_size <- 4
      p <- ggplot2::ggplot(count_df, ggplot2::aes(
        x = .data$n_present,
        fill = .data$fill_col
      )) +
        ggplot2::geom_histogram(
          binwidth = binwidth,
          position = "identity",
          alpha = 1,
          colour = NA
        ) +
        ggplot2::scale_y_log10(
          breaks = 10^(0:6),
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          expand = ggplot2::expansion(mult = c(0, .15))
        ) +
        ggplot2::annotation_logticks(sides = "l", scaled = TRUE) +
        ggplot2::scale_fill_identity(guide = "none") +
        ggplot2::labs(
          x = "# of observations",
          y = expression("# of significantly bound peptides (" * log[10] * ")")
        ) +
        ggplot2::geom_segment(
          data = thresholds_df,
          ggplot2::aes(
            x = .data$thresh, xend = .data$n_samples,
            y = .data$y_line, yend = .data$y_line
          ),
          inherit.aes = FALSE,
          linetype = "dashed",
          color = "black",
          linewidth = 0.4,
          arrow = ggplot2::arrow(length = grid::unit(0.1, "cm"), ends = "both")
        ) +
        ggplot2::geom_text(
          data = thresholds_df,
          ggplot2::aes(
            x = .data$x_mid, y = .data$y_line,
            label = paste0(
              .data$n_peptides_thresh, " peptides in ≥",
              round(prevalence_threshold * 100), "%"
            )
          ),
          inherit.aes = FALSE, vjust = -0.6, size = label_size
        ) +
        ggplot2::geom_text(
          data = thresholds_df,
          ggplot2::aes(
            x = .data$x_mid, y = .data$y_line,
            label = paste0(.data$n_overall, " peptides overall")
          ),
          inherit.aes = FALSE, vjust = 1.6, size = label_size
        ) +
        theme_phip() +
        ggplot2::ggtitle(sprintf("Enrichment counts by %s", title_label))

      if (!is.null(group_col)) {
        p <- p + ggplot2::facet_wrap(~Cohort, ncol = 2, scales = "free_x")
      }

      .ph_log_ok("Plot built")
      p
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# -- alpha-diversity helpers ----------------------------------------------------
.norm_metric <- function(metric) {
  allowed <- c("richness", "shannon_diversity", "simpson_diversity")
  key <- tolower(metric[1])
  if (key %in% allowed) {
    return(key)
  }
  .ph_abort(
    headline = "Unknown metric.",
    step     = "argument validation",
    bullets  = sprintf("allowed: %s", paste(allowed, collapse = ", "))
  )
}

.safe_k <- function(x, k_req) {
  n_ux <- length(unique(stats::na.omit(x)))
  k_eff <- max(3L, min(as.integer(k_req), as.integer(n_ux) - 1L))
  if (!is.finite(k_eff) || k_eff < 3L) 3L else k_eff
}

.gam_band_one <- function(d, xvar, yvar,
                          k_req = 7, level = 0.95, nonneg = FALSE) {
  x <- d[[xvar]]
  rng <- range(x, na.rm = TRUE)
  xgrid <- if (rng[1] == rng[2]) sort(unique(x)) else seq(rng[1], rng[2], length.out = 200)
  k_eff <- .safe_k(x, k_req)
  form <- stats::as.formula(paste(yvar, "~ s(", xvar, ", k = ", k_eff, ")", sep = ""))
  fit <- mgcv::gam(form, data = d)
  nd <- data.frame(xgrid)
  names(nd) <- xvar
  pr <- mgcv::predict.gam(fit, newdata = nd, type = "response", se.fit = TRUE)
  z <- stats::qnorm(1 - (1 - level) / 2)
  mu <- as.numeric(pr$fit)
  se <- as.numeric(pr$se.fit)
  lwr <- mu - z * se
  upr <- mu + z * se
  if (isTRUE(nonneg)) {
    mu <- pmax(mu, 0)
    lwr <- pmax(lwr, 0)
  }
  tibble::tibble(.x = nd[[xvar]], .y = mu, lwr = lwr, upr = upr, .k = k_eff)
}

.bootstrap_gam_one <- function(d, xvar, yvar,
                               k_req = 7, R = 500, level = 0.95,
                               seed = NULL, progress = TRUE, nonneg = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  ctr <- .gam_band_one(d, xvar, yvar, k_req = k_req, level = 0.95, nonneg = nonneg)
  xgrid <- ctr$.x
  mu0 <- ctr$.y
  mat <- matrix(NA_real_, nrow = length(xgrid), ncol = R)
  pb <- if (isTRUE(progress)) utils::txtProgressBar(min = 0, max = R, style = 3) else NULL
  for (r in seq_len(R)) {
    dd <- d[sample.int(nrow(d), replace = TRUE), , drop = FALSE]
    k_eff <- .safe_k(dd[[xvar]], k_req)
    form <- stats::as.formula(paste(yvar, "~ s(", xvar, ", k = ", k_eff, ")", sep = ""))
    fit <- try(mgcv::gam(form, data = dd), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      nd <- data.frame(xgrid)
      names(nd) <- xvar
      mat[, r] <- as.numeric(mgcv::predict.gam(fit, newdata = nd, type = "response"))
    }
    if (!is.null(pb)) utils::setTxtProgressBar(pb, r)
  }
  if (!is.null(pb)) close(pb)
  alpha <- 1 - level
  lwr <- apply(mat, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE, type = 6)
  upr <- apply(mat, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, type = 6)
  if (isTRUE(nonneg)) lwr <- pmax(lwr, 0)
  tibble::tibble(.x = xgrid, .y = mu0, lwr = lwr, upr = upr, .k = unique(ctr$.k))
}

.add_phip_scales <- function(p, custom_colors = NULL) {
  if (is.null(custom_colors)) {
    p + scale_color_phip() + scale_fill_phip()
  } else {
    p + ggplot2::scale_color_manual(values = custom_colors) +
      ggplot2::scale_fill_manual(values = custom_colors)
  }
}

.maybe_facet_by_rank <- function(p, df, rank_col, facet_by_rank, ncol, facet_scales) {
  if (!isTRUE(facet_by_rank)) {
    return(p)
  }
  if (is.null(rank_col) || !rank_col %in% names(df)) {
    return(p)
  }
  if (length(unique(df[[rank_col]])) <= 1) {
    return(p)
  }
  p + ggplot2::facet_wrap(stats::as.formula(paste("~", rank_col)),
    ncol = ncol, scales = facet_scales
  )
}

.build_alpha_plot <- function(alpha_df,
                              metric,
                              group_col = "group",
                              rank_col = "rank",
                              filter_groups = NULL,
                              filter_ranks = NULL,
                              custom_colors = NULL,
                              facet_by_rank = TRUE,
                              ncol = 2,
                              facet_scales = "fixed",
                              time_col = NULL,
                              continuous_mode = c("gam", "binned", "loess"),
                              gam_k = 7,
                              point_alpha = 0.25,
                              ci_method = c("model", "bootstrap"),
                              ci_level = 0.95,
                              boot_R = 500,
                              boot_seed = NULL,
                              boot_progress = TRUE,
                              ci_fill = "grey70",
                              ci_alpha = 0.15) {
  metric_col <- .norm_metric(metric)
  continuous_mode <- match.arg(continuous_mode)
  ci_method <- match.arg(ci_method)
  df <- tibble::as_tibble(alpha_df)

  if (!metric_col %in% names(df)) {
    .ph_abort(
      headline = "Metric column not found in alpha data.",
      step     = "input validation",
      bullets  = sprintf("missing: %s", add_quotes(metric_col, 1L))
    )
  }
  if (!is.null(group_col) && !group_col %in% names(df)) {
    .ph_abort(
      headline = "Grouping column not found in alpha data.",
      step     = "input validation",
      bullets  = sprintf("group_col: %s", add_quotes(group_col, 1L))
    )
  }
  if (!is.null(rank_col) && !rank_col %in% names(df)) {
    .ph_warn(
      headline = "Rank column not found in alpha data; ignoring faceting.",
      step     = "input validation",
      bullets  = sprintf("rank_col: %s", add_quotes(rank_col, 1L))
    )
    rank_col <- NULL
  }

  if (!is.null(filter_groups) && !is.null(group_col)) {
    df <- df[df[[group_col]] %in% filter_groups, , drop = FALSE]
  }
  if (!is.null(filter_ranks) && !is.null(rank_col)) {
    df <- df[df[[rank_col]] %in% filter_ranks, , drop = FALSE]
  }

  keep_cols <- c(metric_col, group_col, rank_col, time_col)
  keep_cols <- keep_cols[!is.na(keep_cols)]
  if (length(keep_cols)) {
    df <- df[stats::complete.cases(df[, keep_cols, drop = FALSE]), , drop = FALSE]
  }

  ylab <- switch(metric_col,
    "richness"          = "Richness",
    "shannon_diversity" = "Shannon diversity",
    "simpson_diversity" = "Simpson diversity (1 - \u03A3 p^2)"
  )

  # CROSS-SECTIONAL ------------------------------------------------------------
  if (is.null(time_col) || !time_col %in% names(df)) {
    if (is.null(group_col)) {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = "", y = .data[[metric_col]])) +
        ggplot2::geom_boxplot(fill = "grey70", colour = "black", outlier.shape = NA) +
        ggplot2::geom_jitter(alpha = point_alpha, width = 0.1, size = 1) +
        ggplot2::labs(x = NULL, y = ylab) +
        theme_phip() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank())
      p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
      return(p)
    }

    gsym <- rlang::sym(group_col)
    msym <- rlang::sym(metric_col)
    df_counts <- df |>
      dplyr::group_by(!!gsym) |>
      dplyr::summarise(sample_count = dplyr::n(), .groups = "drop")
    xlab_map <- stats::setNames(
      paste0(df_counts[[group_col]], "\n(n = ", df_counts$sample_count, ")"),
      df_counts[[group_col]]
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = !!gsym, y = !!msym, fill = !!gsym)) +
      ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
      ggplot2::geom_jitter(
        color = "black", size = 1, width = 0.2,
        alpha = point_alpha, show.legend = FALSE
      ) +
      ggplot2::scale_x_discrete(labels = xlab_map) +
      ggplot2::labs(x = "Group", y = ylab, fill = group_col) +
      theme_phip()
    p <- .add_phip_scales(p, custom_colors)
    p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
    return(p)
  }

  # LONGITUDINAL ---------------------------------------------------------------
  if (continuous_mode == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      .ph_abort("mgcv is required for GAM mode.", step = "longitudinal")
    }

    split_keys <- c(group_col, if (!is.null(rank_col) && isTRUE(facet_by_rank)) rank_col else NULL)
    df_split <- if (length(split_keys)) {
      df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(split_keys))) |>
        dplyr::group_split()
    } else {
      list(df)
    }

    .ph_log_info("Fitting GAM smooths (auto-shrinking k per series)",
      bullets = c(
        sprintf("k requested: %d", gam_k),
        sprintf("series: %d", length(df_split))
      )
    )

    preds <- purrr::map_dfr(df_split, function(dsub) {
      dsub <- as.data.frame(dsub)
      out <- if (ci_method == "model") {
        .gam_band_one(dsub,
          xvar = time_col, yvar = metric_col,
          k_req = gam_k, level = ci_level,
          nonneg = (metric_col == "richness")
        )
      } else {
        .bootstrap_gam_one(dsub,
          xvar = time_col, yvar = metric_col,
          k_req = gam_k, R = boot_R, level = ci_level,
          seed = boot_seed, progress = boot_progress,
          nonneg = (metric_col == "richness")
        )
      }
      if (!is.null(group_col)) out[[group_col]] <- unique(dsub[[group_col]])
      if (!is.null(rank_col) && isTRUE(facet_by_rank)) out[[rank_col]] <- unique(dsub[[rank_col]])
      out
    })

    preds$.grp <- if (!is.null(group_col) && !is.null(rank_col) && isTRUE(facet_by_rank)) {
      interaction(preds[[group_col]], preds[[rank_col]], drop = TRUE)
    } else if (!is.null(group_col)) {
      preds[[group_col]]
    } else {
      factor("all")
    }

    p <- ggplot2::ggplot()
    if (!is.null(group_col)) {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]], color = .data[[group_col]]),
        alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]]),
        alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0, colour = "black"
      )
    }

    p <- p +
      ggplot2::geom_ribbon(
        data = preds,
        ggplot2::aes(x = .data$.x, ymin = .data$lwr, ymax = .data$upr, group = .data$.grp),
        inherit.aes = FALSE, fill = ci_fill, alpha = ci_alpha, colour = NA
      ) +
      {
        if (!is.null(group_col)) {
          ggplot2::geom_line(
            data = preds,
            ggplot2::aes(
              x = .data$.x, y = .data$.y,
              color = .data[[group_col]], group = .data$.grp
            ),
            linewidth = 1
          )
        } else {
          ggplot2::geom_line(
            data = preds,
            ggplot2::aes(x = .data$.x, y = .data$.y, group = .data$.grp),
            linewidth = 1, colour = "black"
          )
        }
      } +
      ggplot2::labs(x = time_col, y = ylab, color = group_col) +
      theme_phip()

    if (!is.null(group_col)) p <- .add_phip_scales(p, custom_colors)
    p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
    return(p)
  }

  if (continuous_mode == "binned") {
    rng <- range(df[[time_col]], na.rm = TRUE)
    nbins <- 20
    binwidth <- (rng[2] - rng[1]) / nbins
    make_bins <- function(x, width, xmin) floor((x - xmin) / width)
    df$.bin <- make_bins(df[[time_col]], binwidth, rng[1])
    df$.bin_mid <- rng[1] + (df$.bin + 0.5) * binwidth

    group_map <- c(
      if (!is.null(group_col)) group_col else NULL,
      if (!is.null(rank_col) && isTRUE(facet_by_rank)) rank_col else NULL,
      ".bin", ".bin_mid"
    )

    sdat <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_map))) |>
      dplyr::summarise(
        y = mean(.data[[metric_col]], na.rm = TRUE),
        sd = stats::sd(.data[[metric_col]], na.rm = TRUE),
        n = dplyr::n(), .groups = "drop"
      )

    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    sdat <- sdat |>
      dplyr::mutate(
        se = sd / sqrt(pmax(n, 1)),
        lwr = y - z * se, upr = y + z * se
      )

    sdat$.grp <- if (!is.null(group_col) && !is.null(rank_col) && isTRUE(facet_by_rank)) {
      interaction(sdat[[group_col]], sdat[[rank_col]], drop = TRUE)
    } else if (!is.null(group_col)) {
      sdat[[group_col]]
    } else {
      factor("all")
    }

    p <- ggplot2::ggplot()
    if (!is.null(group_col)) {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]], color = .data[[group_col]]),
        alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]]),
        alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0, colour = "black"
      )
    }

    p <- p +
      ggplot2::geom_ribbon(
        data = sdat,
        ggplot2::aes(x = .data$.bin_mid, ymin = .data$lwr, ymax = .data$upr, group = .data$.grp),
        fill = ci_fill, alpha = ci_alpha, colour = NA, inherit.aes = FALSE
      ) +
      {
        if (!is.null(group_col)) {
          ggplot2::geom_line(
            data = sdat,
            ggplot2::aes(
              x = .data$.bin_mid, y = .data$y,
              color = .data[[group_col]], group = .data$.grp
            ),
            linewidth = 1
          )
        } else {
          ggplot2::geom_line(
            data = sdat,
            ggplot2::aes(x = .data$.bin_mid, y = .data$y, group = .data$.grp),
            linewidth = 1, colour = "black"
          )
        }
      } +
      ggplot2::labs(x = time_col, y = ylab, color = group_col) +
      theme_phip()

    if (!is.null(group_col)) p <- .add_phip_scales(p, custom_colors)
    p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
    return(p)
  }

  # LOESS fallback
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]])) +
    {
      if (!is.null(group_col)) {
        ggplot2::aes(color = .data[[group_col]])
      } else {
        ggplot2::aes()
      }
    } +
    ggplot2::geom_point(alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, span = 0.75, level = ci_level) +
    ggplot2::labs(x = time_col, y = ylab, color = group_col) +
    theme_phip()

  if (!is.null(group_col)) p <- .add_phip_scales(p, custom_colors)
  p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
  p
}
