#' @export
# generic (so the name stays the same)
plot_enrichment_counts <- function(features_target,
                                   group_col,
                                   group_cols,
                                   prevalence_threshold = 0,
                                   custom_colors,
                                   binwidth = 1,
                                   ...) {
  UseMethod("plot_enrichment_counts")
}

# your existing data.frame version can remain as-is:
# plot_enrichment_counts.data.frame <- function(features_target, ...) { <your code> }

# --- phip_data method: derive counts from x$data_long$exist --------------------
#' @param sample_meta optional data.frame with columns c("sample_id", group_col)
#'   to be used if group_col is not present in x$data_long
#' @export
plot_enrichment_counts.phip_data <- function(features_target, # actually <phip_data> x
                                             group_col,
                                             group_cols = group_col,
                                             prevalence_threshold = 0,
                                             custom_colors,
                                             binwidth = 1,
                                             sample_meta = NULL,
                                             ...) {
  x <- features_target
  stopifnot(inherits(x, "phip_data"))
  .data <- rlang::.data
  grp <- rlang::sym(group_col)

  tbl <- x$data_long

  # attach external labels if needed
  if (!group_col %in% colnames(tbl)) {
    stopifnot(
      !is.null(sample_meta),
      all(c("sample_id", group_col) %in% colnames(sample_meta))
    )
    lab <- sample_meta |>
      dplyr::select(sample_id, dplyr::all_of(group_col)) |>
      dplyr::distinct()
    tbl <- tbl |>
      dplyr::left_join(lab, by = "sample_id")
  }

  # sanity checks
  need <- c("sample_id", "peptide_id", "exist", group_col)
  miss <- setdiff(need, colnames(tbl))
  if (length(miss)) stop("Missing columns in data_long: ", paste(miss, collapse = ", "))

  # keep only rows that count as "present" (exist > threshold; with 0/1 this means ==1)
  tbl_present <- tbl |>
    dplyr::filter(.data$exist > !!prevalence_threshold)

  # cohort sizes (how many samples per cohort)
  cohort_sizes <- tbl |>
    dplyr::distinct(.data$sample_id, !!grp) |>
    dplyr::count(Cohort = !!grp, name = "n_samples")

  # per-cohort per-peptide present counts (how many individuals where peptide is present)
  pep_counts <- tbl_present |>
    dplyr::group_by(Cohort = !!grp, .data$peptide_id) |>
    dplyr::summarise(n_present = dplyr::n_distinct(.data$sample_id), .groups = "drop") |>
    dplyr::filter(.data$n_present > 0)

  # thresholds per cohort (5%) and how many peptides exceed it
  thresholds <- cohort_sizes |>
    dplyr::mutate(thresh = ceiling(.data$n_samples * 0.05)) |>
    dplyr::left_join(
      pep_counts |>
        dplyr::group_by(.data$Cohort) |>
        dplyr::summarise(n_peptides5 = sum(.data$n_present >= dplyr::first(ceiling(NA_real_))), .groups = "drop"),
      by = "Cohort"
    )

  # The line above used a placeholder to align schema; recompute n_peptides5 properly after collect
  cohort_sizes_df <- cohort_sizes |> dplyr::collect()
  pep_counts_df <- pep_counts |> dplyr::collect()

  thresholds_df <- cohort_sizes_df |>
    dplyr::mutate(thresh = ceiling(.data$n_samples * 0.05)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      n_peptides5 = sum(pep_counts_df$n_present[pep_counts_df$Cohort == .data$Cohort] >= .data$thresh)
    ) |>
    dplyr::ungroup()

  # plotting data
  real_order <- cohort_sizes_df$Cohort
  count_df <- pep_counts_df |>
    dplyr::mutate(Cohort = factor(.data$Cohort, levels = real_order))

  # ---- plot (same look & feel as your original) ------------------------------
  ggplot2::ggplot(count_df, ggplot2::aes(x = .data$n_present, fill = .data$Cohort)) +
    ggplot2::geom_histogram(
      binwidth = binwidth,
      position = "identity",
      alpha    = 0.9
    ) +
    ggplot2::scale_y_log10(
      breaks = 10^(0:6),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      expand = ggplot2::expansion(mult = c(0, .15))
    ) +
    ggplot2::annotation_logticks(sides = "l", scaled = TRUE) +
    {
      if (missing(custom_colors) || is.null(custom_colors)) {
        ggplot2::scale_fill_discrete()
      } else {
        ggplot2::scale_fill_manual(values = custom_colors)
      }
    } +
    ggplot2::labs(
      x = "# of individuals",
      y = expression("# of significantly bound peptides (" * log[10] * ")")
    ) +
    ggplot2::geom_segment(
      data = thresholds_df,
      ggplot2::aes(
        x = .data$thresh, xend = .data$n_samples,
        y = .data$n_peptides5, yend = .data$n_peptides5
      ),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "black",
      size = 0.4,
      arrow = ggplot2::arrow(length = grid::unit(0.1, "cm"), ends = "both")
    ) +
    ggplot2::geom_text(
      data = thresholds_df,
      ggplot2::aes(
        x = (.data$thresh + .data$n_samples) / 2,
        y = .data$n_peptides5,
        label = paste0(.data$n_peptides5, " peptides in \u22655%")
      ),
      inherit.aes = FALSE,
      vjust = -0.5,
      size = 4
    ) +
    ggplot2::facet_wrap(~Cohort, ncol = 2, scales = "free_x") +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 14, colour = "black"),
      panel.grid.major = ggplot2::element_line(color = "grey90", linetype = "solid"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y.left = ggplot2::element_text(size = 13),
      axis.text.x.bottom = ggplot2::element_text(size = 13),
      axis.title.y = ggplot2::element_text(size = 14),
      axis.title.x = ggplot2::element_text(size = 14, margin = ggplot2::margin(t = 0)),
      plot.margin = ggplot2::margin(0, 3, 0, 0, unit = "pt")
    )
}

#' @export
plot_enrichment_counts_many <- function(x, group_cols, save_dir = NULL, ...) {
  stopifnot(inherits(x, "phip_data"))
  purrr::map(group_cols, function(gc) {
    p <- plot_enrichment_counts(
      x,
      group_col  = gc,
      group_cols = group_cols,
      ...
    )
    if (!is.null(save_dir)) {
      fs::dir_create(save_dir)
      ggplot2::ggsave(file.path(save_dir, paste0("enrichment_", gc, ".pdf")),
        p,
        width = 8, height = 6, limitsize = FALSE
      )
    }
    print(p)
    p
  }) |> rlang::set_names(group_cols)
}

# Ensure peptide_library is queryable from the SAME connection as data_long
.ensure_peplib_on_main <- function(x, schema_alias = "peplib") {
  # main counts connection
  main_con <- dbplyr::remote_con(x$data_long)

  # peptide library connection (prefer explicit meta$peptide_con)
  pep_con <- if (!is.null(x$meta$peptide_con)) x$meta$peptide_con else dbplyr::remote_con(x$peptide_library)

  # try zero-copy ATTACH when both are DuckDB connections
  if (inherits(main_con, "duckdb_connection") && inherits(pep_con, "duckdb_connection")) {
    pep_db_path <- try(pep_con@driver@dbdir, silent = TRUE)
    if (!inherits(pep_db_path, "try-error") && is.character(pep_db_path) && nzchar(pep_db_path)) {
      # ATTACH (ignore "already attached" errors)
      try(DBI::dbExecute(main_con, sprintf("ATTACH '%s' AS %s;", pep_db_path, schema_alias)), silent = TRUE)

      # get the base table name (e.g., "peptide_meta")
      base_name <- tryCatch(
        {
          nm <- dbplyr::remote_name(x$peptide_library)
          if (is.null(nm) || !nzchar(nm)) "peptide_meta" else sub("^.*\\.", "", nm)
        },
        error = function(e) "peptide_meta"
      )

      # IMPORTANT: use a subquery to avoid DBI::dbExistsTable(schema, table) checks
      return(
        dplyr::tbl(
          main_con,
          dbplyr::sql(sprintf("SELECT * FROM %s.%s", schema_alias, base_name))
        )
      )
    }
  }

  # Fallback: copy peptidelib into main_con as a TEMP table
  peplib_local <- dplyr::collect(x$peptide_library)
  tmp_name <- paste0("peptide_meta_tmp_", as.integer(Sys.time()))
  dplyr::copy_to(main_con, peplib_local, tmp_name, temporary = TRUE, overwrite = TRUE)
}

# generic
#' @export
compute_alpha_diversity <- function(x, ...) UseMethod("compute_alpha_diversity")

# phip_data method
#' @rdname compute_alpha_diversity
#' @export
compute_alpha_diversity.phip_data <- function(x,
                                              group_col,
                                              ranks = "peptide",
                                              present_by = c("exist", "fold_change"),
                                              fc_threshold = 0,
                                              sample_meta = NULL,
                                              shannon_log = c("ln", "log2", "log10"),
                                              carry_cols = NULL) {
  stopifnot(inherits(x, "phip_data"))
  present_by <- match.arg(present_by)
  shannon_log <- match.arg(shannon_log)

  # symbols
  sid <- rlang::sym("sample_id")
  pid <- rlang::sym("peptide_id")
  grp <- rlang::sym(group_col)
  rv <- rlang::sym("rank_val")

  # base table
  tbl <- x$data_long

  # attach external labels if needed (also bring carry_cols if provided there)
  if (!group_col %in% colnames(tbl)) {
    stopifnot(!is.null(sample_meta), all(c("sample_id", group_col) %in% colnames(sample_meta)))
    lab_cols <- unique(c("sample_id", group_col, carry_cols))
    lab <- sample_meta |>
      dplyr::select(dplyr::all_of(lab_cols)) |>
      dplyr::distinct()
    tbl <- tbl |>
      dplyr::left_join(lab, by = "sample_id")
  }

  # presence rule (no .data pronoun in lazy context)
  if (present_by == "exist") {
    stopifnot("exist" %in% colnames(tbl))
    tbl <- tbl |> dplyr::filter(!!rlang::sym("exist") == 1L)
  } else {
    stopifnot("fold_change" %in% colnames(tbl))
    tbl <- tbl |> dplyr::filter(!!rlang::sym("fold_change") > !!fc_threshold)
  }

  # ensure carry_cols are present
  if (!is.null(carry_cols)) {
    carry_cols <- intersect(carry_cols, colnames(tbl))
  }

  # distinct present (sample, peptide, cohort)
  tbl_present <- tbl |>
    dplyr::transmute(!!sid, !!pid, cohort = !!grp) |>
    dplyr::distinct()

  # ranks: accept synonyms & split
  ranks[ranks %in% c("peptide_id", "peptideID", "peptide_ids")] <- "peptide"
  ranks_nonpep <- setdiff(unique(ranks), "peptide")

  # prepare a compact rank map on the main con (collect + copy_to, robust across DBs)
  map_tbl <- NULL
  if (length(ranks_nonpep)) {
    peplib_cols <- colnames(x$peptide_library)
    ranks_ok <- intersect(ranks_nonpep, peplib_cols)
    ranks_miss <- setdiff(ranks_nonpep, ranks_ok)
    if (length(ranks_miss) && rlang::is_installed("cli")) {
      cli::cli_warn("Skipping ranks not found in peptide_library: {paste(ranks_miss, collapse=', ')}")
    }

    if (length(ranks_ok)) {
      main_con <- dbplyr::remote_con(x$data_long)
      map_df <- x$peptide_library |>
        dplyr::select(!!pid, dplyr::all_of(ranks_ok)) |>
        dplyr::collect()
      tmp_name <- paste0("peplib_map_tmp_", as.integer(Sys.time()))
      map_tbl <- dplyr::copy_to(main_con, map_df, tmp_name, temporary = TRUE, overwrite = TRUE)
    }
    ranks <- c("peptide", ranks_ok)
  } else {
    ranks <- "peptide"
  }
  ranks <- unique(ranks)

  # Shannon base change factor (H_b = H_ln / ln(b))
  ln_base <- switch(shannon_log,
    ln = 1.0,
    log2 = log(2),
    log10 = log(10)
  )

  # helper: compute one rank
  compute_one <- function(rank) {
    if (identical(rank, "peptide")) {
      ranked <- tbl_present |>
        dplyr::transmute(!!sid, cohort, rank_val = !!pid)
    } else {
      rank_sym <- rlang::sym(rank)
      ranked <- tbl_present |>
        dplyr::inner_join(
          map_tbl |> dplyr::select(!!pid, rank_val = !!rank_sym),
          by = "peptide_id"
        ) |>
        dplyr::filter(!is.na(!!rv))
    }

    # counts per (sample, rank_val) on DB
    per_cat <- ranked |>
      dplyr::group_by(!!sid, cohort, !!rv) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop")

    # bring compact counts to R and finish metrics per sample_id, cohort
    pc <- per_cat |> dplyr::collect()

    if (nrow(pc)) {
      by_sample <- pc |>
        dplyr::group_by(sample_id, cohort) |>
        dplyr::summarise(
          richness = dplyr::n_distinct(rank_val),
          H_ln = {
            p <- n / sum(n)
            -sum(p * log(p))
          },
          simpson = {
            p <- n / sum(n)
            1 - sum(p * p)
          },
          .groups = "drop"
        ) |>
        dplyr::mutate(shannon = H_ln / ln_base) |>
        dplyr::select(-H_ln)
    } else {
      by_sample <- tibble::tibble(
        sample_id = character(0), cohort = character(0),
        richness = integer(0), shannon = numeric(0), simpson = numeric(0)
      )
    }

    # all samples (carry meta columns if requested)
    keep_cols <- c("sample_id", group_col, carry_cols)
    keep_cols <- intersect(keep_cols, colnames(tbl))
    all_samples <- tbl |>
      dplyr::distinct(dplyr::across(dplyr::all_of(keep_cols))) |>
      dplyr::collect() |>
      dplyr::rename(group = !!grp)

    out <- all_samples |>
      dplyr::left_join(by_sample, by = c("sample_id", "group" = "cohort")) |>
      dplyr::mutate(
        richness = tidyr::replace_na(richness, 0L),
        shannon  = tidyr::replace_na(shannon, 0),
        simpson  = tidyr::replace_na(simpson, 0)
      ) |>
      dplyr::mutate(rank = rank, .before = 1)

    # normalize column names to snake_case/lowercase
    norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
    names(out) <- norm(names(out))

    # rename metrics to *_diversity and order columns
    out <- out |>
      dplyr::rename(shannon_diversity = shannon, simpson_diversity = simpson)

    ord <- c(
      "sample_id", "group", setdiff(norm(carry_cols), c("sample_id", "group")),
      "rank", "richness", "shannon_diversity", "simpson_diversity"
    )
    ord <- intersect(ord, names(out))
    out[, ord, drop = FALSE]
  }

  dplyr::bind_rows(lapply(ranks, compute_one))
}

#' Plot alpha diversity (cross-sectional or longitudinal) → returns a ggplot
#' Requires: ggplot2, dplyr, tidyr; and mgcv (+ mvtnorm) when ci_method = "posterior"
#' @export
plot_alpha_diversity <- function(alpha_df,
                                 metric = c(
                                   "richness", "shannon_diversity", "simpson_diversity",
                                   "Richness", "Shannon Diversity", "Simpson Diversity"
                                 ),
                                 group_col = "group",
                                 rank_col = "rank",
                                 custom_colors = NULL,
                                 sig_level = 0.05,
                                 label_format = "p.signif",
                                 facet_by_rank = TRUE,
                                 ncol = 2,
                                 facet_scales = "fixed",
                                 # longitudinal:
                                 time_col = NULL,
                                 continuous_mode = c("gam", "binned", "loess"),
                                 # gam options
                                 gam_k = 7,
                                 gam_se = TRUE,
                                 # binned options
                                 nbins = 20,
                                 binwidth = NULL,
                                 ci_level = 0.95,
                                 stat = c("mean", "median"),
                                 # CI engine
                                 ci_method = c("model", "bootstrap", "posterior"),
                                 # bootstrap knobs
                                 boot_R = 500,
                                 boot_seed = 42,
                                 boot_id_col = NULL, # cluster/block bootstrap if provided
                                 boot_ci = c("percentile", "basic"),
                                 boot_progress = TRUE,
                                 boot_sparse_fallback = c("loess", "skip"),
                                 boot_min_rep_frac = 0.6,
                                 boot_smooth_ci = TRUE,
                                 # posterior bands (GAM)
                                 posterior_grid_n = 200,
                                 # points on continuous plots
                                 point_alpha = 0.25,
                                 # visuals / constraints
                                 enforce_nonneg = NULL, # NULL → TRUE for richness
                                 ci_fill = "grey70", # neutral ribbons
                                 ci_alpha = 0.15) {
  # ---- arg handling ----
  continuous_mode <- match.arg(continuous_mode)
  ci_method <- match.arg(ci_method)
  boot_ci <- match.arg(boot_ci)
  boot_sparse_fallback <- match.arg(boot_sparse_fallback)
  stat <- match.arg(stat)

  metric_map <- c(
    "richness" = "richness",
    "shannon diversity" = "shannon_diversity", "shannon_diversity" = "shannon_diversity",
    "simpson diversity" = "simpson_diversity", "simpson_diversity" = "simpson_diversity"
  )
  metric_col <- metric_map[[tolower(metric[1])]]
  if (is.null(metric_col)) stop("Unknown metric: ", metric[1])

  # default: enforce non-negativity for richness only
  if (is.null(enforce_nonneg)) enforce_nonneg <- (metric_col == "richness")

  gsym <- rlang::sym(group_col)
  msym <- rlang::sym(metric_col)
  rsym <- if (!is.null(rank_col) && rank_col %in% names(alpha_df)) rlang::sym(rank_col) else NULL

  ylab <- switch(metric_col,
    "richness" = "Richness",
    "shannon_diversity" = "Shannon diversity",
    "simpson_diversity" = "Simpson diversity (1 - \u03A3 p^2)"
  )

  df <- alpha_df
  if (!is.null(rsym)) df[[rank_col]] <- factor(df[[rank_col]])

  add_color_scales <- function(p) {
    if (is.null(custom_colors)) {
      p
    } else {
      p + ggplot2::scale_color_manual(values = custom_colors) +
        ggplot2::scale_fill_manual(values = custom_colors)
    }
  }
  add_facets <- function(p) {
    if (!facet_by_rank || is.null(rsym) || length(unique(df[[rank_col]])) <= 1) {
      return(p)
    }
    p + ggplot2::facet_wrap(ggplot2::vars(!!rsym), ncol = ncol, scales = facet_scales)
  }

  # ---------- helpers ----------
  .resample_idx <- function(d) {
    if (!is.null(boot_id_col) && boot_id_col %in% names(d)) {
      ids <- unique(d[[boot_id_col]])
      samp_ids <- sample(ids, length(ids), replace = TRUE)
      unlist(lapply(samp_ids, function(id) which(d[[boot_id_col]] == id)), use.names = FALSE)
    } else {
      sample.int(nrow(d), replace = TRUE)
    }
  }
  .make_newdata <- function(xvar, xgrid) {
    nd <- data.frame(xgrid)
    names(nd) <- xvar
    nd
  }
  .pred_grid <- function(x, n = 200) {
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) sort(unique(x)) else seq(rng[1], rng[2], length.out = n)
  }

  # ---- bootstrap smoother (per-group) ----
  .bootstrap_smoother <- function(d, xvar, yvar,
                                  method = c("gam", "loess"),
                                  k = 7, R = 500, level = 0.95,
                                  seed = NULL, progress = TRUE,
                                  nonneg = FALSE,
                                  sparse_fallback = c("loess", "skip"),
                                  min_rep_frac = 0.6,
                                  smooth_ci = TRUE) {
    method <- match.arg(method)
    sparse_fallback <- match.arg(sparse_fallback)
    if (!is.null(seed)) set.seed(seed)

    x <- d[[xvar]]
    xgrid <- .pred_grid(x)

    # center line
    if (method == "gam") {
      n_ux <- length(unique(x))
      k_eff <- max(3L, min(k, n_ux - 1L))
      if (!is.finite(k_eff) || k_eff < 3L) {
        fit0 <- try(stats::loess(stats::as.formula(paste(yvar, "~", xvar)), data = d, span = 0.75), silent = TRUE)
        mu0 <- if (!inherits(fit0, "try-error")) {
          as.numeric(predict(fit0, newdata = .make_newdata(xvar, xgrid)))
        } else {
          rep(NA_real_, length(xgrid))
        }
      } else {
        form <- stats::as.formula(paste(yvar, "~ s(", xvar, ", k = ", k_eff, ")", sep = ""))
        fit0 <- mgcv::gam(form, data = d)
        mu0 <- as.numeric(predict(fit0, newdata = .make_newdata(xvar, xgrid), type = "response"))
      }
    } else {
      fit0 <- stats::loess(stats::as.formula(paste(yvar, "~", xvar)), data = d, span = 0.75)
      mu0 <- as.numeric(predict(fit0, newdata = .make_newdata(xvar, xgrid)))
    }
    if (isTRUE(nonneg)) mu0 <- pmax(mu0, 0)

    # bootstrap
    mat <- matrix(NA_real_, nrow = length(xgrid), ncol = R)
    pb <- if (isTRUE(progress)) utils::txtProgressBar(min = 0, max = R, style = 3) else NULL

    for (r in seq_len(R)) {
      idx <- .resample_idx(d)
      dd <- d[idx, , drop = FALSE]
      if (method == "gam") {
        n_ux <- length(unique(dd[[xvar]]))
        k_eff <- max(3L, min(k, n_ux - 1L))
        if (!is.finite(k_eff) || k_eff < 3L) {
          if (sparse_fallback == "skip") {
            if (!is.null(pb)) utils::setTxtProgressBar(pb, r)
            next
          } else {
            fit_r <- try(stats::loess(stats::as.formula(paste(yvar, "~", xvar)), data = dd, span = 0.75), silent = TRUE)
            if (!inherits(fit_r, "try-error")) {
              mat[, r] <- as.numeric(predict(fit_r, newdata = .make_newdata(xvar, xgrid)))
            }
          }
        } else {
          form <- stats::as.formula(paste(yvar, "~ s(", xvar, ", k = ", k_eff, ")", sep = ""))
          fit_r <- try(mgcv::gam(form, data = dd), silent = TRUE)
          if (!inherits(fit_r, "try-error")) {
            mat[, r] <- as.numeric(predict(fit_r, newdata = .make_newdata(xvar, xgrid), type = "response"))
          }
        }
      } else {
        fit_r <- try(stats::loess(stats::as.formula(paste(yvar, "~", xvar)), data = dd, span = 0.75), silent = TRUE)
        if (!inherits(fit_r, "try-error")) {
          mat[, r] <- as.numeric(predict(fit_r, newdata = .make_newdata(xvar, xgrid)))
        }
      }
      if (!is.null(pb)) utils::setTxtProgressBar(pb, r)
    }
    if (!is.null(pb)) close(pb)

    alpha <- 1 - level
    coverage <- rowSums(!is.na(mat))
    lwr <- apply(mat, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE, type = 6)
    upr <- apply(mat, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, type = 6)

    # coverage guard + optional smoothing
    min_reps <- ceiling(R * min_rep_frac)
    lwr[coverage < min_reps] <- NA_real_
    upr[coverage < min_reps] <- NA_real_
    if (isTRUE(smooth_ci)) {
      smooth_vec <- function(y) {
        pred <- try(stats::predict(stats::loess(y ~ seq_along(y),
          span = 0.25,
          na.action = stats::na.exclude
        )), silent = TRUE)
        if (inherits(pred, "try-error")) y else as.numeric(pred)
      }
      lwr <- smooth_vec(lwr)
      upr <- smooth_vec(upr)
    }
    if (isTRUE(nonneg)) lwr <- pmax(lwr, 0)

    tibble::tibble(.x = xgrid, .y = mu0, lwr = lwr, upr = upr)
  }

  # ---- posterior band for GAM (per-group) ----
  .posterior_band <- function(d, xvar, yvar, k = 7, level = 0.95,
                              nonneg = FALSE, grid_n = 200) {
    if (!requireNamespace("mgcv", quietly = TRUE) ||
      !requireNamespace("mvtnorm", quietly = TRUE)) {
      stop("ci_method='posterior' requires 'mgcv' and 'mvtnorm'.")
    }
    x <- d[[xvar]]
    rng <- range(x, na.rm = TRUE)
    xgrid <- if (rng[1] == rng[2]) sort(unique(x)) else seq(rng[1], rng[2], length.out = grid_n)

    n_ux <- length(unique(x))
    k_eff <- max(3L, min(k, n_ux - 1L))
    if (!is.finite(k_eff) || k_eff < 3L) {
      fit0 <- stats::loess(stats::as.formula(paste(yvar, "~", xvar)), data = d, span = 0.75)
      mu0 <- as.numeric(predict(fit0, newdata = .make_newdata(xvar, xgrid)))
      res <- d[[yvar]] - predict(fit0)
      se <- rep(stats::sd(res, na.rm = TRUE), length(xgrid))
      z <- stats::qnorm(1 - (1 - level) / 2)
      lwr <- mu0 - z * se
      upr <- mu0 + z * se
      if (isTRUE(nonneg)) {
        mu0 <- pmax(mu0, 0)
        lwr <- pmax(lwr, 0)
      }
      return(tibble::tibble(.x = xgrid, .y = mu0, lwr = lwr, upr = upr))
    }

    form <- stats::as.formula(paste(yvar, "~ s(", xvar, ", k = ", k_eff, ")", sep = ""))
    fit <- mgcv::gam(form, data = d)
    nd <- .make_newdata(xvar, xgrid)
    Xp <- mgcv::predict.gam(fit, newdata = nd, type = "lpmatrix")
    beta <- stats::coef(fit)
    Vb <- fit$Vp
    R <- 2000L
    B <- mvtnorm::rmvnorm(R, mean = beta, sigma = Vb)
    Y <- Xp %*% t(B)
    mu0 <- as.numeric(Xp %*% beta)
    alpha <- 1 - level
    lwr <- apply(Y, 1, stats::quantile, probs = alpha / 2, type = 6)
    upr <- apply(Y, 1, stats::quantile, probs = 1 - alpha / 2, type = 6)
    if (isTRUE(nonneg)) {
      mu0 <- pmax(mu0, 0)
      lwr <- pmax(lwr, 0)
    }
    tibble::tibble(.x = xgrid, .y = mu0, lwr = lwr, upr = upr)
  }

  # ---------- longitudinal mode ----------
  if (!is.null(time_col) && time_col %in% names(df)) {
    tsym <- rlang::sym(time_col)
    tvec <- df[[time_col]]
    is_cont <- is.numeric(tvec) || inherits(tvec, c("Date", "POSIXct", "POSIXt"))

    # clean rows
    df <- df[stats::complete.cases(df[, c(group_col, metric_col, time_col), drop = FALSE]), , drop = FALSE]

    if (!is_cont) {
      # categorical time → mean ± SD per timepoint & group
      group_vars <- list(gsym, tsym)
      if (!is.null(rsym) && facet_by_rank) group_vars <- c(group_vars, list(rsym))
      sum_df <- df |>
        dplyr::group_by(!!!group_vars) |>
        dplyr::summarise(
          mean = mean(!!msym, na.rm = TRUE),
          sd = stats::sd(!!msym, na.rm = TRUE),
          n = dplyr::n(), .groups = "drop"
        ) |>
        dplyr::mutate(sd = dplyr::coalesce(sd, 0))

      p <- ggplot2::ggplot(
        sum_df,
        ggplot2::aes(x = !!tsym, y = .data$mean, group = !!gsym, color = !!gsym)
      ) +
        ggplot2::geom_line(alpha = 0.9) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            ymin = pmax(.data$mean - .data$sd, 0),
            ymax = .data$mean + .data$sd
          ),
          width = 0.15, alpha = 0.7
        ) +
        ggplot2::labs(x = time_col, y = ylab, color = group_col, fill = group_col) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.6, hjust = 0.5)
        )
      p <- add_color_scales(p)
      p <- add_facets(p)
      return(p)
    }

    # ---- continuous time ----
    if (continuous_mode == "gam") {
      if (ci_method == "model") {
        # model-based (geom_smooth) — ggplot groups by color automatically
        group_keys <- list(gsym)
        if (!is.null(rsym) && facet_by_rank) group_keys <- c(group_keys, list(rsym))
        uniq_per_curve <- df |>
          dplyr::group_by(!!!group_keys) |>
          dplyr::summarise(n_u = dplyr::n_distinct(!!tsym), .groups = "drop")
        nmin <- if (nrow(uniq_per_curve)) min(uniq_per_curve$n_u, na.rm = TRUE) else 0L
        k_safe <- max(3L, min(gam_k, nmin - 1L))
        if (!is.finite(k_safe) || k_safe < 3L) {
          warning(sprintf("Not enough unique time points for GAM (min=%s); falling back to LOESS.", nmin))
          p <- ggplot2::ggplot(df, ggplot2::aes(x = !!tsym, y = !!msym, color = !!gsym)) +
            ggplot2::geom_point(alpha = point_alpha, size = 1) +
            ggplot2::geom_smooth(method = "loess", se = TRUE, span = 0.75, level = ci_level) +
            ggplot2::labs(x = time_col, y = ylab, color = group_col) +
            ggplot2::theme_bw(base_size = 13) +
            ggplot2::theme(panel.grid = ggplot2::element_blank())
          p <- add_color_scales(p)
          p <- add_facets(p)
          return(p)
        }
        p <- ggplot2::ggplot(df, ggplot2::aes(x = !!tsym, y = !!msym, color = !!gsym)) +
          ggplot2::geom_point(alpha = point_alpha, size = 1) +
          ggplot2::geom_smooth(
            method = "gam",
            formula = y ~ s(x, k = k_safe),
            se = gam_se,
            level = ci_level
          ) +
          ggplot2::labs(x = time_col, y = ylab, color = group_col) +
          ggplot2::theme_bw(base_size = 13) +
          ggplot2::theme(panel.grid = ggplot2::element_blank())
        p <- add_color_scales(p)
        p <- add_facets(p)
        return(p)
      }

      if (ci_method == "posterior") {
        if (!requireNamespace("mgcv", quietly = TRUE) ||
          !requireNamespace("mvtnorm", quietly = TRUE)) {
          stop("ci_method='posterior' requires 'mgcv' and 'mvtnorm'.")
        }
        split_keys <- c(group_col, if (!is.null(rsym) && facet_by_rank) rank_col else NULL)
        df_split <- df |>
          dplyr::group_by(dplyr::across(dplyr::all_of(split_keys))) |>
          dplyr::group_split()

        preds <- purrr::map_dfr(df_split, function(dsub) {
          dsub <- as.data.frame(dsub)
          gval <- unique(dsub[[group_col]])
          rval <- if (!is.null(rsym) && facet_by_rank) unique(dsub[[rank_col]]) else NA
          out <- .posterior_band(dsub,
            xvar = time_col, yvar = metric_col,
            k = gam_k, level = ci_level,
            nonneg = enforce_nonneg,
            grid_n = posterior_grid_n
          )
          out[[group_col]] <- gval
          if (!is.null(rsym) && facet_by_rank) out[[rank_col]] <- rval
          out
        }) |>
          dplyr::mutate(
            .grp = if (!is.null(rsym) && facet_by_rank) {
              interaction(.data[[group_col]], .data[[rank_col]], drop = TRUE)
            } else {
              .data[[group_col]]
            }
          ) |>
          dplyr::arrange(.grp, .data$.x) |>
          tidyr::drop_na(lwr, upr)

        p <- ggplot2::ggplot() +
          ggplot2::geom_point(
            data = df,
            ggplot2::aes(x = !!tsym, y = !!msym, color = !!gsym),
            alpha = point_alpha, size = 1
          ) +
          ggplot2::geom_ribbon(
            data = preds,
            ggplot2::aes(
              x = .data$.x, ymin = .data$lwr, ymax = .data$upr,
              group = .data$.grp
            ),
            fill = ci_fill, alpha = ci_alpha, colour = NA,
            inherit.aes = FALSE, na.rm = TRUE
          ) +
          ggplot2::geom_line(
            data = preds,
            ggplot2::aes(x = .data$.x, y = .data$.y, color = !!gsym, group = .data$.grp),
            linewidth = 1, na.rm = TRUE
          ) +
          ggplot2::labs(x = time_col, y = ylab, color = group_col) +
          ggplot2::theme_bw(base_size = 13) +
          ggplot2::theme(panel.grid = ggplot2::element_blank())
        p <- add_color_scales(p)
        p <- add_facets(p)
        return(p)
      }

      # ci_method == "bootstrap"
      if (!requireNamespace("mgcv", quietly = TRUE)) stop("mgcv is required for GAM bootstrap.")
      set.seed(boot_seed)
      split_keys <- c(group_col, if (!is.null(rsym) && facet_by_rank) rank_col else NULL)
      df_split <- df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(split_keys))) |>
        dplyr::group_split()

      preds <- purrr::map_dfr(df_split, function(dsub) {
        dsub <- as.data.frame(dsub)
        gval <- unique(dsub[[group_col]])
        rval <- if (!is.null(rsym) && facet_by_rank) unique(dsub[[rank_col]]) else NA
        out <- .bootstrap_smoother(
          dsub,
          xvar = time_col, yvar = metric_col,
          method = "gam", k = gam_k, R = boot_R,
          level = ci_level, seed = NULL,
          progress = boot_progress,
          nonneg = enforce_nonneg,
          sparse_fallback = boot_sparse_fallback,
          min_rep_frac = boot_min_rep_frac,
          smooth_ci = boot_smooth_ci
        )
        out[[group_col]] <- gval
        if (!is.null(rsym) && facet_by_rank) out[[rank_col]] <- rval
        out
      }) |>
        dplyr::mutate(
          .grp = if (!is.null(rsym) && facet_by_rank) {
            interaction(.data[[group_col]], .data[[rank_col]], drop = TRUE)
          } else {
            .data[[group_col]]
          }
        ) |>
        dplyr::arrange(.grp, .data$.x) |>
        tidyr::drop_na(lwr, upr)

      p <- ggplot2::ggplot() +
        ggplot2::geom_point(
          data = df,
          ggplot2::aes(x = !!tsym, y = !!msym, color = !!gsym),
          alpha = point_alpha, size = 1
        ) +
        ggplot2::geom_ribbon(
          data = preds,
          ggplot2::aes(
            x = .data$.x, ymin = .data$lwr, ymax = .data$upr,
            group = .data$.grp
          ),
          fill = ci_fill, alpha = ci_alpha, colour = NA,
          inherit.aes = FALSE, na.rm = TRUE
        ) +
        ggplot2::geom_line(
          data = preds,
          ggplot2::aes(x = .data$.x, y = .data$.y, color = !!gsym, group = .data$.grp),
          linewidth = 1, na.rm = TRUE
        ) +
        ggplot2::labs(x = time_col, y = ylab, color = group_col) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
      p <- add_color_scales(p)
      p <- add_facets(p)
      return(p)
    }

    if (continuous_mode == "binned") {
      rng <- range(df[[time_col]], na.rm = TRUE)
      if (is.null(binwidth)) binwidth <- (rng[2] - rng[1]) / nbins
      make_bins <- function(x, width, xmin) floor((x - xmin) / width)
      df$.bin <- make_bins(df[[time_col]], binwidth, rng[1])
      df$.bin_mid <- rng[1] + (df$.bin + 0.5) * binwidth

      group_vars <- list(gsym, rlang::sym(".bin"), rlang::sym(".bin_mid"))
      if (!is.null(rsym) && facet_by_rank) group_vars <- c(group_vars, list(rsym))

      sum_fun <- if (stat == "mean") mean else stats::median

      sdat <- df |>
        dplyr::group_by(!!!group_vars) |>
        dplyr::summarise(
          y = sum_fun(!!msym, na.rm = TRUE),
          sd = stats::sd(!!msym, na.rm = TRUE),
          n = dplyr::n(),
          .groups = "drop"
        )

      if (ci_method == "model") {
        z <- stats::qnorm(1 - (1 - ci_level) / 2)
        sdat <- sdat |>
          dplyr::mutate(
            se = sd / sqrt(pmax(n, 1)),
            lwr = y - z * se,
            upr = y + z * se
          )
      } else {
        set.seed(boot_seed)
        split_keys <- c(group_col, if (!is.null(rsym) && facet_by_rank) rank_col else NULL, ".bin", ".bin_mid")
        df_split <- df |>
          dplyr::group_by(dplyr::across(dplyr::all_of(split_keys))) |>
          dplyr::summarise(vals = list(!!msym), .groups = "drop")
        .boot_bin <- function(vals) {
          n <- length(vals)
          if (!n) {
            return(c(NA_real_, NA_real_))
          }
          fun <- if (stat == "mean") mean else stats::median
          boots <- replicate(boot_R, fun(sample(vals, n, replace = TRUE)))
          alpha <- 1 - ci_level
          stats::quantile(boots, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE, type = 6)
        }
        boots <- df_split |>
          dplyr::mutate(ci = purrr::map(vals, .boot_bin)) |>
          tidyr::unnest_wider(ci, names_sep = "_") |>
          dplyr::rename(lwr = ci_1, upr = ci_2) |>
          dplyr::select(-vals)
        sdat <- sdat |>
          dplyr::left_join(boots, by = c(group_col, if (!is.null(rsym) && facet_by_rank) rank_col, ".bin", ".bin_mid"))
      }

      if (isTRUE(enforce_nonneg)) {
        sdat <- sdat |>
          dplyr::mutate(y = pmax(y, 0), lwr = pmax(lwr, 0))
      }

      # per-cohort polygon grouping for ribbons
      sdat <- sdat |>
        dplyr::mutate(
          .grp = if (!is.null(rsym) && facet_by_rank) {
            interaction(.data[[group_col]], .data[[rank_col]], drop = TRUE)
          } else {
            .data[[group_col]]
          }
        ) |>
        dplyr::arrange(.grp, .data$.bin_mid) |>
        tidyr::drop_na(lwr, upr)

      p <- ggplot2::ggplot(sdat, ggplot2::aes(
        x = .data$.bin_mid, y = .data$y,
        color = !!gsym, group = .grp
      )) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lwr, ymax = .data$upr, group = .grp),
          inherit.aes = FALSE, fill = ci_fill, alpha = ci_alpha,
          colour = NA, na.rm = TRUE
        ) +
        ggplot2::geom_line(size = 1, na.rm = TRUE) +
        ggplot2::geom_point(size = 1.5, na.rm = TRUE) +
        ggplot2::labs(x = time_col, y = ylab, color = group_col) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
      p <- add_color_scales(p)
      p <- add_facets(p)
      return(p)
    }

    # LOESS (plain)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = !!tsym, y = !!msym, color = !!gsym)) +
      ggplot2::geom_point(alpha = point_alpha, size = 1) +
      ggplot2::geom_smooth(method = "loess", se = TRUE, span = 0.75, level = ci_level) +
      ggplot2::labs(x = time_col, y = ylab, color = group_col) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    p <- add_color_scales(p)
    p <- add_facets(p)
    return(p)
  }

  # ---------- cross-sectional mode ----------
  df <- df[stats::complete.cases(df[, c(group_col, metric_col), drop = FALSE]), , drop = FALSE]

  df_counts <- df |>
    dplyr::group_by(!!gsym) |>
    dplyr::summarise(sample_count = dplyr::n(), .groups = "drop")
  xlab_map <- stats::setNames(
    paste0(df_counts[[group_col]], "\n(n = ", df_counts$sample_count, ")"),
    df_counts[[group_col]]
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = !!gsym, y = !!msym, fill = !!gsym)) +
    ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
    ggplot2::geom_jitter(color = "black", size = 1, width = 0.2, alpha = 0.3, show.legend = FALSE) +
    ggplot2::scale_x_discrete(labels = xlab_map) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(x = "Group", y = ylab, fill = group_col) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x        = ggplot2::element_text(angle = 45, vjust = 0.6, hjust = 0.5),
      axis.text.x.bottom = ggplot2::element_text(size = 13),
      axis.text.y.left   = ggplot2::element_text(size = 13),
      axis.title.y       = ggplot2::element_text(size = 13),
      axis.title.x       = ggplot2::element_text(size = 13),
      plot.margin        = ggplot2::margin(0, 1, 0, 1, unit = "pt"),
      panel.grid         = ggplot2::element_blank()
    )

  if (!is.null(custom_colors)) {
    p <- p + ggplot2::scale_fill_manual(values = custom_colors)
  }

  # pairwise Wilcoxon (sig only)
  levs <- levels(factor(df[[group_col]]))
  if (length(levs) >= 2 && rlang::is_installed("ggpubr")) {
    pairwise <- combn(levs, 2, simplify = FALSE)
    sig_pairs <- purrr::keep(pairwise, function(pr) {
      x <- df[df[[group_col]] == pr[1], metric_col, drop = TRUE]
      y <- df[df[[group_col]] == pr[2], metric_col, drop = TRUE]
      stats::wilcox.test(x, y, exact = FALSE)$p.value < sig_level
    })
    if (length(sig_pairs)) {
      p <- p + ggpubr::stat_compare_means(
        method      = "wilcox.test",
        comparisons = sig_pairs,
        label       = label_format,
        hide.ns     = FALSE,
        size        = 4.5,
        tip.length  = 0.02
      )
    }
  }
  p <- add_facets(p)
  return(p)
}
