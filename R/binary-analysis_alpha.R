# -------------------------------------------------------------------
# PLOTTING - generic + methods with group_interaction support
# -------------------------------------------------------------------

#' @title Plot enrichment counts by cohort/group
#'
#' @description For each cohort (level of a grouping variable), draw a histogram
#'   of how many individuals a peptide is present in. A dashed horizontal line
#'   marks a **prevalence threshold** equal to `ceiling(n_samples *
#'   prevalence_threshold)`, with two labels:
#'
#' - “**N peptides in ≥threshold**” (above the dashed line)
#' - “**M peptides overall**” (below; peptides with ≥1 present sample)
#'
#' @details Works with a `<phip_data>` object or a long `data.frame`.
#' - If `group_cols` is a character vector, a **separate plot** is created for
#'   each grouping column.
#' - If `group_cols = NULL`, a **single non-facetted** plot over all samples is
#'   returned.
#' - If `group_interaction = TRUE` and `length(group_cols) >= 2`, an additional
#'   plot for the **interaction** of all `group_cols` is appended.
#'
#'   Styling defaults to `theme_phip()` and the **PHIP palette**. Colors are
#'   applied deterministically by attaching hex values to the data and using
#'   `scale_fill_identity()`, preventing ggplot from overriding hues.
#'
#' @param phip_data A `<phip_data>` object **or** a long `data.frame` with
#'   columns `sample_id`, `peptide_id`, `exist` (0/1 or logical), and the
#'   grouping column(s) when `group_cols` is not `NULL`.
#' @param group_cols Character vector of grouping column names, or `NULL`
#'   (no facets; one aggregate plot).
#' @param prevalence_threshold Numeric in `[0, 1]`. Fraction of cohort size used
#'   for the dashed horizontal **prevalence** line (default `0.05`). **Note**:
#'   presence per sample is defined as `exist > 0` and is not affected by this
#'   parameter.
#' @param custom_colors Optional named vector of hex colors to override the PHIP
#'   palette (used only if provided).
#' @param binwidth Histogram bin width for `# of individuals` (x-axis).
#' @param group_interaction Logical; also plot the **interaction** of all
#'   `group_cols` (default `FALSE`).
#' @param interaction_sep Separator string between interaction parts
#'   (default `" × "`).
#' @param ... Reserved.
#'
#' @return A single `ggplot` (when a single plot is produced) or a named list
#'   of `ggplot` objects (one per grouping column, plus optional interaction).
#'
#' @examples
#' \dontrun{
#' # phip_data input
#' plot_enrichment_counts(pd, group_cols = "Cohort",
#'                        prevalence_threshold = 0.05)
#'
#' # no facets (aggregate over all samples)
#' plot_enrichment_counts(pd, group_cols = NULL)
#'
#' # data.frame input + interaction
#' plot_enrichment_counts(df_long,
#'   group_cols = c("Cohort", "timepoint"),
#'   group_interaction = TRUE
#' )
#' }
#' @export
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

# ------------------------------------------------------------------------------
# Internal helper (shared): build one plot for one group var (or none)
# ------------------------------------------------------------------------------
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
  # clean eval
  .data <- rlang::.data
  grp_sym <- if (!is.null(group_col)) rlang::sym(group_col) else NULL

  # timing the plotting for logs
  .ph_with_timing(
    headline = "Building enrichment count plot",
    step = if (is.null(group_col)) {
      "no grouping (aggregate)"
    } else {
      sprintf("grouping variable: '%s'", group_col)
    },
    expr = {

      # ---- 0) sanity checks --------------------------------------------------
      need <- c("sample_id", "peptide_id", "exist")
      if (!is.null(group_col)) need <- c(need, group_col)
      miss <- setdiff(need, colnames(tbl))
      if (length(miss)) {
        .ph_abort(
          headline = "Missing required columns for plotting.",
          step     = "input validation",
          bullets  = sprintf("missing: %s",
                             paste(add_quotes(miss, 1L), collapse = ", "))
        )
      }

      # ---- ensure a 'Cohort' column is present (for unified downstream code) -
      tbl <- if (is.null(group_col)) {
        dplyr::mutate(tbl, Cohort = "All samples")
      } else {
        dplyr::mutate(tbl, Cohort = !!grp_sym)
      }

      # ---- 1) present rows ---------------------------------------------------
      # Presence is defined as exist > 0 (binary/logical compatible).
      tbl_present <- dplyr::filter(tbl, .data$exist > 0)

      # ---- 2) cohort sizes ---------------------------------------------------
      cohort_sizes <- tbl |>
        dplyr::distinct(.data$sample_id, .data$Cohort) |>
        dplyr::count(.data$Cohort, name = "n_samples") |>
        dplyr::collect()

      # rather a fallback warning
      if (nrow(cohort_sizes) == 0L || any(cohort_sizes$n_samples <= 0)) {
        .ph_warn(
          headline = "No samples for at least one cohort (facet may be empty).",
          step     = "cohort sizes"
        )
        cohort_sizes <- dplyr::filter(cohort_sizes, .data$n_samples > 0)
      }

      # ---- 3) per-cohort per-peptide counts ----------------------------------
      # aka fast and cheap peptide prevalence per group/Cohort
      # by selecting data only with distinct sample_ids we keep it
      # subject(+ timepoint)-specific
      # NOTE: it is on distinct sample_id!!! --> this means for cross-sectional
      # data basically each subject is treated separately, but for longitudinal
      # each subject x timepoint combo is treated separately; because of that in
      # some groupings the number of observations may seem to be higher than
      # the number of distinc subject_id --> basically because each sub x time
      # is treated separately, so the number of observations will be actually
      # the number of unique subject x timepoint combos !!!!
      pep_counts <- tbl_present |>
        dplyr::group_by(.data$Cohort, .data$peptide_id) |>
        dplyr::summarise(n_present = dplyr::n_distinct(.data$sample_id),
                         .groups = "drop") |>
        dplyr::filter(.data$n_present > 0) |>
        dplyr::collect()

      # ---- 4) thresholds & labels (use prevalence_threshold proportion) ------
      # filtering the peptide prevalence by the subject treshold
      thresholds <- cohort_sizes |>
        dplyr::mutate(thresh = ceiling(.data$n_samples * prevalence_threshold))

      # calculate how many peptides above/below treshold
      n_thresh_tbl <- pep_counts |>
        dplyr::inner_join(dplyr::select(thresholds, .data$Cohort, .data$thresh),
                          by = "Cohort") |>
        dplyr::group_by(.data$Cohort) |>
        dplyr::summarise(
          n_peptides_thresh = sum(.data$n_present >= .data$thresh),
          .groups = "drop")

      overall_tbl <- pep_counts |>
        dplyr::group_by(.data$Cohort) |>
        dplyr::summarise(n_overall = dplyr::n_distinct(.data$peptide_id),
                         .groups = "drop")

      # build the table for plotting
      thresholds_df <- thresholds |>
        dplyr::left_join(n_thresh_tbl, by = "Cohort") |>
        dplyr::left_join(overall_tbl, by = "Cohort") |>
        dplyr::mutate(
          n_peptides_thresh = dplyr::coalesce(.data$n_peptides_thresh, 0L),
          n_overall         = dplyr::coalesce(.data$n_overall, 0L),
          y_line            = pmax(.data$n_peptides_thresh, 1L), # for log scale
          x_mid             = (.data$thresh + .data$n_samples) / 2
        )

      # ---- 5) plotting data --------------------------------------------------
      # set facet order to the order in thresholds_df
      real_order <- thresholds_df$Cohort
      count_df <- pep_counts |>
        dplyr::mutate(Cohort = factor(.data$Cohort, levels = real_order))

      # fresh, named PHIP palette mapping (per plot)
      lvl <- levels(count_df$Cohort)
      pal_map <- if (is.null(custom_colors)) {
        .phip_palette_map(lvl)
      } else {
        # respect user colors; if not named, recycle in order
        if (is.null(names(custom_colors))) {
          stats::setNames(rep_len(custom_colors, length(lvl)), lvl)
        } else {
          stats::setNames(unname(custom_colors[lvl]), lvl)
        }
      }

      count_df <- dplyr::mutate(count_df,
                                fill_col = unname(pal_map[
                                  as.character(.data$Cohort)]))

      # ---- 6) plot -----------------------------------------------------------
      label_size <- 4

      p <- ggplot2::ggplot(count_df, ggplot2::aes(x = .data$n_present,
                                                  fill = .data$fill_col)) +
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
        ggplot2::annotation_logticks(sides = "l",
                                     scaled = TRUE) +
        ggplot2::scale_fill_identity(guide = "none") +
        ggplot2::labs(
          x = "# of observations",
          y = expression("# of significantly bound peptides (" * log[10] * ")")
        ) +
        ggplot2::geom_segment(
          data = thresholds_df,
          ggplot2::aes(x = .data$thresh,
                       xend = .data$n_samples,
                       y = .data$y_line,
                       yend = .data$y_line),
          inherit.aes = FALSE,
          linetype = "dashed",
          color = "black",
          linewidth = 0.4,
          arrow = ggplot2::arrow(length = grid::unit(0.1, "cm"), ends = "both")
        ) +
        ggplot2::geom_text(
          data = thresholds_df,
          ggplot2::aes(
            x = .data$x_mid,
            y = .data$y_line,
            label = paste0(.data$n_peptides_thresh,
                           " peptides in ≥",
                           round(prevalence_threshold * 100), "%")
          ),
          inherit.aes = FALSE,
          vjust = -0.6,
          size = label_size
        ) +
        ggplot2::geom_text(
          data = thresholds_df,
          ggplot2::aes(
            x = .data$x_mid,
            y = .data$y_line,
            label = paste0(.data$n_overall, " peptides overall")
          ),
          inherit.aes = FALSE,
          vjust = 1.6,
          size = label_size
        ) +
        theme_phip() +
        ggplot2::ggtitle(sprintf("Enrichment counts by %s", title_label))

      # facet only when we actually have a grouping variable
      if (!is.null(group_col)) {
        p <- p + ggplot2::facet_wrap(~Cohort, ncol = 2, scales = "free_x")
      }

      .ph_log_ok("Plot built")
      p
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ------------------------------------------------------------------------------
# phip_data method
# ------------------------------------------------------------------------------

#' @rdname plot_enrichment_counts
#'
#' @section Grouping columns:
#' For the `<phip_data>` method, all requested `group_cols` must be present in
#' `x$data_long`. If any grouping column is missing, the function aborts with a
#' phiper-style error explaining which columns are missing.
#'
#' @export
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

  # timing for logs
  .ph_with_timing(
    headline = "Plotting enrichment counts (<phip_data>)",
    step     = if (is.null(group_cols)) "group_cols: <none>" else
      sprintf("group_cols: %s", paste(add_quotes(group_cols, 1L),
                                      collapse = ", ")),
    expr = {
      # extract data_long
      tbl <- x$data_long

      ## prune introduced zeros if we know we're on a full-cross ---------------
      if (isTRUE(x$meta$full_cross) && ("exist" %in% colnames(tbl))) {
        red_txt <- tryCatch({
          ep <- as.numeric(x$meta$exist_prop)
          if (is.finite(ep) && ep > 0) sprintf("~%.1fx", 1/ep) else "<unknown>"
        }, error = function(e) "<unknown>")
        .ph_log_info(
          "Full-cross detected; pruning non-existing rows before plotting",
          bullets = c("rule: keep exist == 1",
                      sprintf("estimated reduction: %s", red_txt))
        )
        tbl <- dplyr::filter(tbl, .data$exist == 1L)
      }

      # if grouping is requested, ensure all columns exist in data_long
      if (!is.null(group_cols) && length(group_cols)) {
        missing_gcs <- setdiff(group_cols, colnames(tbl))
        if (length(missing_gcs)) {
          .ph_abort(
            headline = "Missing grouping columns in data_long.",
            step     = "input validation",
            bullets  = sprintf("missing: %s",
                               paste(add_quotes(missing_gcs, 1L),
                                     collapse = ", "))
          )
        }
      }

      # plotting dispatcher
      if (is.null(group_cols) || !length(group_cols)) {
        # single aggregate plot (no facets)
        .plot_enrichment_counts_one(
          tbl,
          group_col = NULL,
          prevalence_threshold = prevalence_threshold,
          custom_colors = custom_colors,
          binwidth = binwidth
        )
      } else {
        # one plot per group column
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

        # optional: interaction plot across ALL provided group_cols
        if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
          combo_nm  <- paste(group_cols, collapse = interaction_sep)
          inter_col <- "..phip_interaction.."
          tbl_inter <- dplyr::mutate(tbl, !!rlang::sym(inter_col) :=
                                       paste(!!!rlang::syms(group_cols),
                                             sep = interaction_sep))

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

# ------------------------------------------------------------------------------
# data.frame method
# ------------------------------------------------------------------------------

#' @rdname plot_enrichment_counts
#' @export
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

  # timing for logs
  .ph_with_timing(
    headline = "Plotting enrichment counts (data.frame)",
    step = if (is.null(group_cols)) {
      "group_cols: <none>"
    } else{
      sprintf("group_cols: %s", paste(add_quotes(group_cols, 1L),
                                      collapse = ", "))
    },
    expr = {
      ## minimal checks
      need_base <- c("sample_id", "peptide_id", "exist")
      miss_base <- setdiff(need_base, colnames(tbl))
      if (length(miss_base)) {
        .ph_abort(
          headline = "Missing required columns in data.frame.",
          step     = "input validation",
          bullets  = sprintf("missing: %s",
                             paste(add_quotes(miss_base, 1L),
                                   collapse = ", "))
        )
      }

      # plot dispatcher (bascially the same as up)
      if (is.null(group_cols) || !length(group_cols)) {
        .plot_enrichment_counts_one(
          tbl,
          group_col = NULL,
          prevalence_threshold = prevalence_threshold,
          custom_colors = custom_colors,
          binwidth = binwidth
        )
      } else {
        # one plot per group column
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

        # optional interaction across ALL group_cols
        if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
          combo_nm <- paste(group_cols, collapse = interaction_sep)
          inter_col <- "..phip_interaction.."
          tbl_inter <- dplyr::mutate(tbl, !!rlang::sym(inter_col) :=
                                       paste(!!!rlang::syms(group_cols),
                                             sep = interaction_sep))

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

# ==============================================================================
# Alpha diversity - generic + methods
# ==============================================================================

#' @title Compute alpha diversity per sample / group across ranks
#'
#' @description Computes **richness**, **Shannon**, and **Simpson** diversity
#'   per sample and grouping variable at one or more **ranks** (columns) that
#'   describe peptides.
#'
#' @details ## What are “ranks” here? Ranks are *the peptide identities or
#' characteristics you aggregate by*. They must be **exact column names**:
#' - For `<phip_data>`: columns in the PHIPER peptide library from Vogl Lab
#' (e.g., `peptide_id`, lineage/taxa fields available in your library).
#' - For `data.frame`: columns present on your long count table.
#'
#' ## Presence rule
#' - Default: `exist > 0`.
#' - If `fc_threshold` is numeric, presence is instead `fold_change > fc_threshold`.
#'
#' ## Grouping & interactions
#' - `group_cols` can be a character vector; the return value is a **named list**
#' of data frames, one per `group_col`.
#' - If `group_cols = NULL`, a single non-facetted table is returned under the
#' name `"all_samples"`.
#' - If `group_interaction = TRUE` (and you supplied ≥ 2 `group_cols`), an
#' additional element is computed for the interaction of all group columns, with
#' labels joined by `interaction_sep`.
#'
#' @param x A `<phip_data>` object or a long `data.frame`.
#' @param group_cols Character vector of grouping columns, or `NULL` for a
#'   single aggregate (non-facetted) table. **All columns must be present** on
#'   the input.
#' @param ranks Character vector of **exact column names** to aggregate by (see
#'   “What are ranks?”). Typical defaults include `"peptide_id"` or taxonomy /
#'   lineage columns present in your peptide library or data frame.
#' @param fc_threshold Numeric or `NULL`. If `NULL` (default), presence is
#'   `exist > 0`. If numeric, presence is `fold_change > fc_threshold`.
#' @param shannon_log One of `"ln"`, `"log2"`, `"log10"`; reporting base for the
#'   Shannon index (via base change from natural log).
#' @param carry_cols Optional character vector of extra columns to carry forward
#'   into the output if present (e.g., sample metadata already in your table).
#' @param group_interaction Logical; also compute the **interaction** of all
#'   `group_cols` (default `FALSE`).
#' @param interaction_sep Separator used for the interaction label (default `" *
#'   "`).
#'
#' @return A **named list** of data frames with S3 class
#'   `"phip_alpha_diversity"`. Each element (per `group_col`, plus optional
#'   interaction or `"all_samples"`) contains: `sample_id`, the grouping column
#'   (or `group` when `group_cols = NULL`), any `carry_cols`, and the metrics:
#'   `rank`, `richness`, `shannon_diversity`, `simpson_diversity`.
#'
#' @examples
#' \dontrun{
#' # phip_data input — peptide-level diversity by cohort
#' out <- compute_alpha_diversity(pd,
#'   group_cols = "Cohort",
#'   ranks = "peptide_id"
#' )
#'
#' # Use ranks present in the peptide library (e.g., lineage/taxa)
#' out_taxa <- compute_alpha_diversity(pd,
#'   group_cols = c("Cohort","timepoint"),
#'   ranks = c("peptide_id","family","genus"),
#'   group_interaction = TRUE
#' )
#'
#' # data.frame input — ranks must be columns in the data
#' out_df <- compute_alpha_diversity(df_long,
#'   group_cols = NULL,
#'   ranks = "peptide_id"
#' )
#'
#' # Presence via fold-change
#' out_fc <- compute_alpha_diversity(pd,
#'   group_cols = "Cohort",
#'   ranks = "peptide_id",
#'   fc_threshold = 1.5
#' )
#' }
#' @export
compute_alpha_diversity <- function(x, ...) UseMethod("compute_alpha_diversity")

# ------------------------------------------------------------------------------
# internals
# ------------------------------------------------------------------------------

# natural->base-b change factor for Shannon (H_b = H_ln / ln(b))
.phip_ln_base <- function(shannon_log = c("ln","log2","log10")) {
  shannon_log <- match.arg(shannon_log)
  switch(shannon_log,
         ln    = 1.0,
         log2  = log(2),
         log10 = log(10)
  )
}

# Engine used by both methods.
# - tbl: counts table (lazy or local) containing sample_id, peptide_id, exist,
#        optional fold_change, group_col(s).
# - map_provider(rank_name): returns a two-column table (peptide_id, rank_val),
#        lazy on same con (phip_data) or local tibble (data.frame). Must return
#        NULL to skip a rank (with a warning already emitted upstream).
.compute_alpha_for_group <- function(tbl,
                                     group_col = NULL,
                                     ranks,
                                     fc_threshold = NULL,
                                     shannon_log = c("ln","log2","log10"),
                                     carry_cols = NULL,
                                     map_provider) {
  # clean eval and helpers
  .data <- rlang::.data
  shannon_log <- match.arg(shannon_log)
  ln_base <- .phip_ln_base(shannon_log)

  # ---- validate required columns ---------------------------------------------
  need <- c("sample_id", "peptide_id", "exist")
  if (!is.null(fc_threshold)) need <- union(need, "fold_change")
  if (!is.null(group_col))    need <- union(need, group_col)
  miss <- setdiff(need, colnames(tbl))
  if (length(miss)) {
    .ph_abort(
      headline = "Missing required columns.",
      step     = "input validation",
      bullets  = sprintf("missing: %s", paste(add_quotes(miss, 1L),
                                              collapse = ", "))
    )
  }

  # ---- unified cohort column -------------------------------------------------
  tbl <- if (is.null(group_col)) {
    dplyr::mutate(tbl, cohort = "All samples")
  } else {
    dplyr::mutate(tbl, cohort = .data[[group_col]])
  }

  # ---- presence rule ---------------------------------------------------------
  pres_tbl <- if (is.null(fc_threshold)) {
    dplyr::filter(tbl, .data$exist > 0)
  } else {
    dplyr::filter(tbl, .data$fold_change > !!fc_threshold)
  }

  # distinct present tuples
  pres_min <- pres_tbl |>
    dplyr::transmute(sample_id = .data$sample_id,
                     peptide_id = .data$peptide_id,
                     cohort     = .data$cohort) |>
    dplyr::distinct()

  # ---- each rank must be an exact column name --------------------------------
  ranks <- unique(ranks)
  if (!length(ranks) || !all(vapply(ranks, is.character, logical(1)))) {
    .ph_abort("`ranks` must be a non-empty character vector
              of exact column names.")
  }

  # ---- per-rank computation --------------------------------------------------
  compute_one_rank <- function(rank_name) {
    # if the rank is exactly 'peptide_id', no mapping/join is needed:
    # use the left table's column directly to define rank_val.
    if (identical(rank_name, "peptide_id")) {
      ranked <- pres_min |>
        dplyr::transmute(
          sample_id,
          cohort,
          rank_val = peptide_id
        )
    } else {
      # obtain (peptide_id, rank_val) mapping lazily on the same connection
      map_tbl <- map_provider(rank_name)
      if (is.null(map_tbl)) return(NULL)  # already warned upstream

      # join & keep non-missing rank values
      ranked <- pres_min |>
        dplyr::inner_join(map_tbl, by = "peptide_id") |>
        dplyr::filter(!is.na(.data$rank_val))
    }

    # counts per (sample, cohort, rank_val) - stays lazy if possible
    per_cat <- ranked |>
      dplyr::group_by(sample_id, cohort, rank_val) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop")

    # compact compute
    pc <- dplyr::collect(per_cat)

    if (nrow(pc)) {
      by_sample <- pc |>
        dplyr::group_by(sample_id, cohort) |>
        dplyr::summarise(
          richness = dplyr::n_distinct(rank_val),
          H_ln     = { p <- n / sum(n); -sum(p * log(p)) },
          simpson  = { p <- n / sum(n); 1 - sum(p * p) },
          .groups  = "drop"
        ) |>
        dplyr::mutate(shannon = H_ln / ln_base) |>
        dplyr::select(-H_ln)
    } else {
      by_sample <- tibble::tibble(
        sample_id = character(0), cohort = character(0),
        richness = integer(0), shannon = numeric(0), simpson = numeric(0)
      )
    }

    # carry meta columns if requested
    keep_cols <- c("sample_id", "cohort", carry_cols)
    keep_cols <- intersect(keep_cols, colnames(tbl))
    all_samples <- tbl |>
      dplyr::distinct(dplyr::across(dplyr::all_of(keep_cols))) |>
      dplyr::collect()

    out <- all_samples |>
      dplyr::left_join(by_sample, by = c("sample_id", "cohort")) |>
      dplyr::mutate(
        richness = tidyr::replace_na(richness, 0L),
        shannon  = tidyr::replace_na(shannon, 0),
        simpson  = tidyr::replace_na(simpson, 0)
      ) |>
      dplyr::mutate(rank = rank_name, .before = 1)

    # normalise names; metrics to *_diversity
    norm <- function(x) gsub("[^a-z0-9]+", "_", tolower(x))
    names(out) <- norm(names(out))
    out |>
      dplyr::rename(shannon_diversity = shannon,
                    simpson_diversity = simpson)
  }


  res <- do.call(dplyr::bind_rows, lapply(ranks, compute_one_rank))

  # friendly grouping column naming
  if (is.null(group_col)) {
    res <- dplyr::rename(res, group = cohort)
  } else {
    res <- dplyr::rename(res, !!rlang::sym(group_col) := cohort)
  }

  res
}

# ------------------------------------------------------------------------------
# <phip_data> method
# ------------------------------------------------------------------------------

#' @rdname compute_alpha_diversity
#' @export
compute_alpha_diversity.phip_data <- function(x,
                                              group_cols = NULL,
                                              ranks = "peptide_id",
                                              fc_threshold = NULL,
                                              shannon_log = c("ln",
                                                              "log2",
                                                              "log10"),
                                              carry_cols = NULL,
                                              group_interaction = FALSE,
                                              interaction_sep = " * ") {
  # fast check and tidy eval
  stopifnot(inherits(x, "phip_data"))
  .data <- rlang::.data

  # log time
  .ph_with_timing(
    headline = "Computing alpha diversity (<phip_data>)",
    step     = paste0(
      "group_cols: ", if (is.null(group_cols)) {
        "<none>"
      } else {
        paste(add_quotes(group_cols, 1L), collapse = ", ")
      },
      "; ranks: ", paste(add_quotes(ranks, 1L), collapse = ", ")
    ),
    expr = {
      tbl <- x$data_long

      ## if we’re working on a full-cross, prune introduced zeros upfront
      if (isTRUE(x$meta$full_cross) && ("exist" %in% colnames(tbl))) {
        # estimate reduction factor if we know exist_prop
        red_txt <- tryCatch({
          ep <- as.numeric(x$meta$exist_prop)
          if (is.finite(ep) && ep > 0) sprintf("~%.1fx", 1/ep) else "<unknown>"
        }, error = function(e) "<unknown>")

        .ph_log_info(
          "Full-cross detected; pruning non-existing rows before alpha calc",
          bullets = c(
            "rule: keep exist == 1",
            sprintf("estimated reduction: %s", red_txt)
          )
        )

        tbl <- dplyr::filter(tbl, .data$exist == 1L)
      }

      # validate requested group columns
      if (!is.null(group_cols) && length(group_cols)) {
        miss_gc <- setdiff(group_cols, colnames(tbl))
        if (length(miss_gc)) {
          .ph_abort(
            headline = "Grouping columns not found in data_long.",
            step     = "input validation",
            bullets  = sprintf("missing: %s", paste(add_quotes(miss_gc, 1L),
                                                    collapse = ", "))
          )
        }
      }

      # peptide-library rank availability (must include ranks) --> extract it
      # from the Vogl's lab peptide library by default
      peplib_main <- .ensure_peplib_on_main(x)  # <- (lazy tbl on main con)
      peplib_cols <- colnames(peplib_main)

      .ph_log_info(
        "Peptide library attached on main connection",
        bullets = c(
          sprintf("available columns: %s%s",
                  paste(utils::head(peplib_cols, 8), collapse = ", "),
                  if (length(peplib_cols) > 8) {
                    sprintf(" …(+%d)", length(peplib_cols) - 8)
                  } else {
                    ""
                  })
        )
      )

      # map provider: return a lazy (same-connection) mapping (peptide_id,
      # rank_val)
      map_provider <- function(rank_name) {
        if (!(rank_name %in% peplib_cols)) {
          .ph_warn(
            headline = "Rank not found in peptide_library (skipping).",
            step     = "rank mapping",
            bullets  = sprintf("rank: %s", add_quotes(rank_name, 1L))
          )
          return(NULL)
        }
        peplib_main |>
          dplyr::select(peptide_id, rank_val = .data[[rank_name]]) |>
          dplyr::distinct()
      }

      out_list <- list()

      # simple dispatcher based on the structure of the provided args
      if (is.null(group_cols) || !length(group_cols)) {
        out_list[["all_samples"]] <-
          .compute_alpha_for_group(
            tbl, group_col = NULL, ranks = ranks,
            fc_threshold = fc_threshold, shannon_log = shannon_log,
            carry_cols = carry_cols, map_provider = map_provider
          )
      } else {
        for (gc in group_cols) {
          out_list[[gc]] <-
            .compute_alpha_for_group(
              tbl, group_col = gc, ranks = ranks,
              fc_threshold = fc_threshold, shannon_log = shannon_log,
              carry_cols = carry_cols, map_provider = map_provider
            )
        }
        if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
          inter_col <- "..phip_interaction.."
          combo_nm  <- paste(group_cols, collapse = interaction_sep)
          tbl_inter <- dplyr::mutate(tbl, !!rlang::sym(inter_col) :=
                                       paste(!!!rlang::syms(group_cols),
                                             sep = interaction_sep))
          out_list[[combo_nm]] <-
            .compute_alpha_for_group(
              tbl_inter, group_col = inter_col, ranks = ranks,
              fc_threshold = fc_threshold, shannon_log = shannon_log,
              carry_cols = carry_cols, map_provider = map_provider
            )
        }
      }

      # Tag output with S3 class + attributes for downstream plotting
      out_list <- lapply(out_list, tibble::as_tibble)
      class(out_list) <- c("phip_alpha_diversity", class(out_list))
      attr(out_list, "group_cols")   <- group_cols
      attr(out_list, "ranks")        <- ranks
      attr(out_list, "fc_threshold") <- fc_threshold
      attr(out_list, "shannon_log")  <- match.arg(shannon_log)

      out_list
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ------------------------------------------------------------------------------
# data.frame method
# ------------------------------------------------------------------------------

#' @rdname compute_alpha_diversity
#' @export
compute_alpha_diversity.data.frame <- function(x,
                                               group_cols = NULL,
                                               ranks = "peptide_id",
                                               fc_threshold = NULL,
                                               shannon_log = c("ln",
                                                               "log2",
                                                               "log10"),
                                               carry_cols = NULL,
                                               group_interaction = FALSE,
                                               interaction_sep = " * ") {
  # clean eval
  tbl <- x
  .data <- rlang::.data

  # timing as always
  .ph_with_timing(
    headline = "Computing alpha diversity (data.frame)",
    step     = paste0(
      "group_cols: ", if (is.null(group_cols)) {
        "<none>"
      } else {
        paste(add_quotes(group_cols, 1L), collapse = ", ")
      },
      "; ranks: ", paste(add_quotes(ranks, 1L), collapse = ", ")
    ),
    expr = {
      # validate group columns
      if (!is.null(group_cols) && length(group_cols)) {
        miss_gc <- setdiff(group_cols, colnames(tbl))
        if (length(miss_gc)) {
          .ph_abort(
            headline = "Grouping columns not found in data.frame.",
            step     = "input validation",
            bullets  = sprintf("missing: %s", paste(add_quotes(miss_gc, 1L),
                                                    collapse = ", "))
          )
        }
      }

      # provider for data.frame: ranks must be exact columns in tbl; no peptide
      # library default
      map_provider <- function(rank_name) {
        if (!(rank_name %in% colnames(tbl))) {
          .ph_warn(
            headline = "Rank not found in data.frame (skipping).",
            step     = "rank mapping",
            bullets  = sprintf("rank: %s", add_quotes(rank_name, 1L))
          )
          return(NULL)
        }
        tibble::tibble(
          peptide_id = tbl$peptide_id,
          rank_val   = tbl[[rank_name]]
        ) |>
          dplyr::distinct()
      }

      out_list <- list()

      # dispatcher --> the same as for phip_data; uses different calls based on
      # the input structure
      if (is.null(group_cols) || !length(group_cols)) {
        out_list[["all_samples"]] <-
          .compute_alpha_for_group(
            tbl, group_col = NULL, ranks = ranks,
            fc_threshold = fc_threshold, shannon_log = shannon_log,
            carry_cols = carry_cols, map_provider = map_provider
          )
      } else {
        for (gc in group_cols) {
          out_list[[gc]] <-
            .compute_alpha_for_group(
              tbl, group_col = gc, ranks = ranks,
              fc_threshold = fc_threshold, shannon_log = shannon_log,
              carry_cols = carry_cols, map_provider = map_provider
            )
        }
        if (isTRUE(group_interaction) && length(group_cols) >= 2L) {
          inter_col <- "..phip_interaction.."
          combo_nm  <- paste(group_cols, collapse = interaction_sep)
          tbl_inter <- dplyr::mutate(tbl, !!rlang::sym(inter_col) :=
                                       paste(!!!rlang::syms(group_cols),
                                             sep = interaction_sep))
          out_list[[combo_nm]] <-
            .compute_alpha_for_group(
              tbl_inter, group_col = inter_col, ranks = ranks,
              fc_threshold = fc_threshold, shannon_log = shannon_log,
              carry_cols = carry_cols, map_provider = map_provider
            )
        }
      }

      out_list <- lapply(out_list, tibble::as_tibble)
      class(out_list) <- c("phip_alpha_diversity", class(out_list))
      attr(out_list, "group_cols")     <- group_cols
      attr(out_list, "ranks")          <- ranks
      attr(out_list, "fc_threshold")   <- fc_threshold
      attr(out_list, "shannon_log")    <- match.arg(shannon_log)

      out_list
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ==============================================================================
# Alpha-diversity plotting (generic + methods) with PHIP styling
# ==============================================================================

#' @title Plot alpha diversity (cross-sectional or longitudinal)
#'
#' @description Visualize per-sample **alpha diversity** (e.g., *richness*,
#' *Shannon*, *Simpson*) across groups and peptide ranks.
#'
#' - For `<phip_data>` and long `data.frame` inputs, the function first calls
#' [compute_alpha_diversity()] to produce a tidy summary and then plots it.
#' - For a precomputed object of class **`phip_alpha_diversity`**, it plots
#' directly (no recomputation).
#'
#' Styling defaults to `theme_phip()` and PHIP colour/fill scales unless you
#' provide `custom_colors`.
#'
#' @section Metrics & ranks:
#' Allowed metrics are **"richness"**, **"shannon_diversity"**,
#' **"simpson_diversity"** (case-insensitive). If you pass several metrics,
#' a **separate plot** is returned for each.
#'
#' Ranks are the peptide features you aggregated over when computing alpha
#' diversity (e.g., `"peptide_id"`, `"genus"`, `"species"`, …). The plotting
#' interface expects the result column to be named `rank` (this is the default
#' of [compute_alpha_diversity()]).
#'
#' @section Longitudinal mode:
#' Provide `time_col` to enable longitudinal plots. For the **`phip_data`**
#' method, longitudinal mode is allowed **only** if `x$meta$longitudinal` is
#' `TRUE`. For `data.frame` and `phip_alpha_diversity`, no such guard is
#' enforced.
#'
#' Continuous curves use GAM smoothing with an **auto-shrunk per-series k**
#' (`k_eff ≤ n_unique(x)-1`, min 3) to avoid over-parameterization. Confidence
#' bands are either **model-based** (`ci_method = "model"`) or **bootstrap**
#' (`"bootstrap"` with progress bar via `boot_progress = TRUE`).
#'
#' @param x Input object: a `<phip_data>`, a long `data.frame`, or a
#'   precomputed **`phip_alpha_diversity`**.
#' @param metric Character vector of metrics to plot. Allowed values:
#'   `"richness"`, `"shannon_diversity"`, `"simpson_diversity"`.
#'   Case-insensitive. If multiple are given, one plot per metric is returned.
#' @param group_col Name of the grouping column in the **alpha** data. For
#'   objects computed by [compute_alpha_diversity()] the default is `"group"`.
#'   Use `NULL` to plot without grouping (single colour/curve).
#' @param rank_col Name of the rank column in the **alpha** data (default
#'   `"rank"`). Used only for faceting and optional filtering.
#' @param filter_groups Optional character vector. If provided, only these
#'   groups are retained **for plotting**.
#' @param filter_ranks Optional character vector. If provided, only these ranks
#'   are retained **for plotting**.
#' @param custom_colors Optional named vector of colours for groups. If omitted,
#'   PHIP palette/scales are used.
#' @param facet_by_rank Logical; facet by `rank` if there is more than one level
#'   (default `TRUE`).
#' @param ncol Facet column count when `facet_by_rank = TRUE` (default `2`).
#' @param facet_scales Facet scales, passed to `ggplot2::facet_wrap()`; one of
#'   `"fixed"`, `"free_x"`, `"free_y"`, `"free"`. Default `"fixed"`.
#'
#' @section Longitudinal options:
#' @param time_col Optional name of a time/age variable (numeric, Date, or
#'   POSIXt). If provided, the plot switches to longitudinal mode.
#' @param continuous_mode One of `"gam"`, `"binned"`, `"loess"` (default
#'   `"gam"`). `"gam"` uses `mgcv::gam()` with per-series safe k. `"binned"`
#'   summarizes within time bins. `"loess"` draws a LOESS smoother.
#' @param gam_k Requested GAM basis size k (per series). The effective k is
#'   auto-shrunk to `k_eff ≤ n_unique(x)-1` (min 3). Default `7`.
#' @param point_alpha Alpha for raw points overlay (default `0.25`; set `0` to
#'   hide points).
#'
#' @section Confidence bands:
#' @param ci_method `"model"` (default) or `"bootstrap"`.
#' @param ci_level Confidence level for bands (default `0.95`).
#' @param boot_R Number of bootstrap replicates when `ci_method="bootstrap"`
#'   (default `500`).
#' @param boot_seed Optional integer seed for bootstrap reproducibility.
#' @param boot_progress Logical; show a textual progress bar during bootstrap
#'   (default `TRUE`).
#' @param ci_fill Fill colour for ribbons (default `"grey70"`).
#' @param ci_alpha Alpha for ribbons (default `0.15`).
#'
#' @return A single `ggplot` object or a **named list** of plots (when multiple
#'   `metric` values are provided).
#'
#' @family phip-plotting
#' @export
plot_alpha_diversity <- function(x, ...) {
  UseMethod("plot_alpha_diversity")
}

# ------------------------------------------------------------------------------
# Helpers (internal)
# ------------------------------------------------------------------------------

# make metric case-insensitive
.norm_metric <- function(metric) {
  allowed <- c("richness", "shannon_diversity", "simpson_diversity")
  key <- tolower(metric[1])
  if (key %in% allowed) return(key)
  .ph_abort(
    headline = "Unknown metric.",
    step     = "argument validation",
    bullets  = sprintf("allowed: %s", paste(allowed, collapse = ", "))
  )
}

# fallback for the gam_k argument -- sometimes the gam_k excceds the maximum
# degrees of freedom available for the longitudinal data (eg when you set gam_k
# to 9 and your data has only 3 unique timepoints) --> take the maximal allowed
# value then
.safe_k <- function(x, k_req) {
  n_ux  <- length(unique(stats::na.omit(x)))
  k_eff <- max(3L, min(as.integer(k_req), as.integer(n_ux) - 1L))
  if (!is.finite(k_eff) || k_eff < 3L) 3L else k_eff
}

# actually a wrapper for the mgcv::gam() smoother; calculates the gam curve over
# a grid of xvar for the given values of yvar in the d
.gam_band_one <- function(d,
                          xvar,
                          yvar,
                          k_req = 7,
                          level = 0.95,
                          nonneg = FALSE) {
  # extract the variable and get the x-grid to smooth/predict on
  x <- d[[xvar]]
  rng   <- range(x, na.rm = TRUE)
  xgrid <- if (rng[1] == rng[2]) {
    sort(unique(x))
  } else {
    seq(rng[1], rng[2], length.out = 200)
  }

  k_eff <- .safe_k(x, k_req) # previous wrapper to get the number of k
  form  <- stats::as.formula(paste(yvar,
                                   "~ s(",
                                   xvar,
                                   ", k = ",
                                   k_eff,
                                   ")", sep = ""))

  fit  <- mgcv::gam(form, data = d)  # workhorse --> fits gam
  nd   <- data.frame(xgrid); names(nd) <- xvar # for the model predictions

  ## get the smoothed gam curve
  pr   <- mgcv::predict.gam(fit, newdata = nd, type = "response", se.fit = TRUE)

  # extract upper/lower bands and the mean (actual prediction --> BLUP)
  z   <- stats::qnorm(1 - (1 - level) / 2)
  mu  <- as.numeric(pr$fit)
  se  <- as.numeric(pr$se.fit)
  lwr <- mu - z * se
  upr <- mu + z * se

  # when nonneg=TRUE, do not allow mu or lower band to go below 0; it doesnt
  # actually make sense for the predictions or bands for count data to be below
  # 0; counts are measured on an absolute scale
  if (isTRUE(nonneg)) { mu <- pmax(mu, 0); lwr <- pmax(lwr, 0) }

  tibble::tibble(.x = nd[[xvar]], .y = mu, lwr = lwr, upr = upr, .k = k_eff)
}

# kinda the same clean logic as up, but the lower/upper bands are calculated
# with bootstrap; it is better for the final runs of plotting, cause the model
# extracted bands tend to be too narrow
.bootstrap_gam_one <- function(d,
                               xvar,
                               yvar,
                               k_req = 7,
                               R = 500,
                               level = 0.95,
                               seed = NULL,
                               progress = TRUE,
                               nonneg = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  # use the previous gam-helper to get the BLUPs, grids and meta
  ctr  <- .gam_band_one(d, xvar, yvar, k_req = k_req, level = 0.95,
                        nonneg = nonneg)
  xgrid <- ctr$.x
  mu0   <- ctr$.y

  # setting up a template matrix for the bootstrap
  mat <- matrix(NA_real_, nrow = length(xgrid), ncol = R)

  # some bp can be long, so a progress bar is nice to have; maybe rewrite this
  # in the future to match the phiper logging style?
  pb  <- if (isTRUE(progress)) {
    utils::txtProgressBar(min = 0, max = R, style = 3)
  } else {
    NULL
  }

  # begin bootstrap for R iters
  for (r in seq_len(R)) {
    dd <- d[sample.int(nrow(d), replace = TRUE), , drop = FALSE] # sampling
    k_eff <- .safe_k(dd[[xvar]], k_req)
    form  <- stats::as.formula(paste(yvar, "~ s(", xvar, ", k = ", k_eff, ")",
                                     sep = ""))
    fit   <- try(mgcv::gam(form, data = dd), silent = TRUE) # fitting
    if (!inherits(fit, "try-error")) {
      nd <- data.frame(xgrid); names(nd) <- xvar # extracting BLUPs
      mat[, r] <- as.numeric(mgcv::predict.gam(fit,
                                               newdata = nd,
                                               type = "response"))
    }
    if (!is.null(pb)) utils::setTxtProgressBar(pb, r) # update progress
  }
  if (!is.null(pb)) close(pb) # close bar

  # get lwr/upr bands based on the quantiles defined by the user; default is
  # standard 0.95
  alpha <- 1 - level
  lwr <- apply(mat, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE,
               type = 6)
  upr <- apply(mat, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE,
               type = 6)
  if (isTRUE(nonneg)) lwr <- pmax(lwr, 0) # dont allow negative lowers

  tibble::tibble(.x = xgrid, .y = mu0, lwr = lwr, upr = upr,
                 .k = unique(ctr$.k))
}

# small helper to deal with the colors (default + customs)
.add_phip_scales <- function(p,
                             custom_colors = NULL) {
  if (is.null(custom_colors)) {
    p + scale_color_phip() + scale_fill_phip()
  } else {
    p + ggplot2::scale_color_manual(values = custom_colors) +
      ggplot2::scale_fill_manual(values = custom_colors)
  }
}

# facet by rank if possible and explicitly defined
.maybe_facet_by_rank <- function(p,
                                 df,
                                 rank_col,
                                 facet_by_rank,
                                 ncol,
                                 facet_scales) {
  if (!isTRUE(facet_by_rank)) return(p)
  if (is.null(rank_col) || !rank_col %in% names(df)) return(p)
  if (length(unique(df[[rank_col]])) <= 1) return(p)
  p + ggplot2::facet_wrap(stats::as.formula(paste("~", rank_col)),
                          ncol = ncol, scales = facet_scales)
}

# ------------------------------------------------------------------------------
# central builder once alpha data is available
# ------------------------------------------------------------------------------
.build_alpha_plot <- function(alpha_df,
                              metric,
                              group_col = "group",
                              rank_col  = "rank",
                              filter_groups = NULL,
                              filter_ranks  = NULL,
                              custom_colors = NULL,
                              facet_by_rank = TRUE,
                              ncol = 2,
                              facet_scales = "fixed",
                              time_col = NULL,
                              continuous_mode = c("gam","binned","loess"),
                              gam_k = 7,
                              point_alpha = 0.25,
                              ci_method = c("model","bootstrap"),
                              ci_level  = 0.95,
                              boot_R = 500,
                              boot_seed = NULL,
                              boot_progress = TRUE,
                              ci_fill = "grey70",
                              ci_alpha = 0.15) {

  # normalize inputs
  metric_col <- .norm_metric(metric)
  continuous_mode <- match.arg(continuous_mode)
  ci_method       <- match.arg(ci_method)

  df <- tibble::as_tibble(alpha_df)

  # basic checks + filtering ---------------------------------------------------
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

  # filtering
  if (!is.null(filter_groups) && !is.null(group_col)) {
    df <- df[df[[group_col]] %in% filter_groups, , drop = FALSE]
  }
  if (!is.null(filter_ranks) && !is.null(rank_col)) {
    df <- df[df[[rank_col]] %in% filter_ranks, , drop = FALSE]
  }

  # clean NA rows in key columns -----------------------------------------------
  keep_cols <- c(metric_col, group_col, rank_col, time_col)
  keep_cols <- keep_cols[!is.na(keep_cols)]
  if (length(keep_cols)) {
    df <- df[stats::complete.cases(df[, keep_cols, drop = FALSE]), ,
             drop = FALSE]
  }

  # aesthetics/labels ----------------------------------------------------------
  ylab <- switch(metric_col,
                 "richness"          = "Richness",
                 "shannon_diversity" = "Shannon diversity",
                 "simpson_diversity" = "Simpson diversity (1 - \u03A3 p^2)"
  )

  # cross-sectional ------------------------------------------------------------
  if (is.null(time_col) || !time_col %in% names(df)) {
    if (is.null(group_col)) {
      ## base plot
      p <- ggplot2::ggplot(df, ggplot2::aes(x = "", y = .data[[metric_col]])) +
        ggplot2::geom_boxplot(fill = "grey70", colour = "black",
                              outlier.shape = NA) +
        ggplot2::geom_jitter(alpha = point_alpha, width = 0.1, size = 1) +
        ggplot2::labs(x = NULL, y = ylab) +
        theme_phip() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank())

      ## apply facets
      p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol,
                                facet_scales)
      return(p)
    }

    # tidy eval of group and metric
    gsym <- rlang::sym(group_col)
    msym <- rlang::sym(metric_col)

    # annotate n per group
    df_counts <- df |>
      dplyr::group_by(!!gsym) |>
      dplyr::summarise(sample_count = dplyr::n(), .groups = "drop")

    # labs mapping
    xlab_map <- stats::setNames(
      paste0(df_counts[[group_col]], "\n(n = ", df_counts$sample_count, ")"),
      df_counts[[group_col]]
    )

    # plotting with group_col defined
    p <- ggplot2::ggplot(df, ggplot2::aes(x = !!gsym, y = !!msym,
                                          fill = !!gsym)) +
      ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE) +
      ggplot2::geom_jitter(color = "black", size = 1, width = 0.2,
                           alpha = point_alpha, show.legend = FALSE) +
      ggplot2::scale_x_discrete(labels = xlab_map) +
      ggplot2::labs(x = "Group", y = ylab, fill = group_col) +
      theme_phip()

    # color scales and facetting
    p <- .add_phip_scales(p, custom_colors)
    p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol,
                              facet_scales)
    return(p)
  }

  # longitudinal ---------------------------------------------------------------
  if (continuous_mode == "gam") { # gam should be in suggests
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      .ph_abort("mgcv is required for GAM mode.", step = "longitudinal")
    }

    # cheracter vector of columns to split by
    split_keys <- c(group_col,
                    if (!is.null(rank_col) && isTRUE(facet_by_rank)) {
                      rank_col
                    } else {
                      NULL
                    })

    # if there's at least one split key :
    # group the data by every unique combination
    # if no, return NULL
    df_split <- if (length(split_keys)) {
      df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(split_keys))) |>
        dplyr::group_split()
    } else {
      list(df)
    }

    # the whole splitting logic is made to fit the model separately to x-groups.
    # eg, when you define the group_col, you want to fit the model separately
    # by group to get the smooths by group --> if you additionally define
    # multiple ranks, you want to split by rank x group
    .ph_log_info("Fitting GAM smooths (auto-shrinking k per series)",
                 bullets = c(sprintf("k requested: %d", gam_k),
                             sprintf("series: %d", length(df_split))))

    # workhorse of the longitudinal workflow; fits the model separately to each
    # group defined in the df_split using the previously defined helper
    preds <- purrr::map_dfr(df_split, function(dsub) {
      dsub <- as.data.frame(dsub)

      # standard, narrow model CI
      out <- if (ci_method == "model") {
        .gam_band_one(dsub, xvar = time_col, yvar = metric_col,
                      k_req = gam_k, level = ci_level,
                      nonneg = (metric_col == "richness"))
      } else {
        # bootstrapped CI
        .bootstrap_gam_one(dsub, xvar = time_col, yvar = metric_col,
                           k_req = gam_k, R = boot_R, level = ci_level,
                           seed = boot_seed, progress = boot_progress,
                           nonneg = (metric_col == "richness"))
      }
      if (!is.null(group_col)) out[[group_col]] <- unique(dsub[[group_col]])
      if (!is.null(rank_col)  && isTRUE(facet_by_rank)) {
        out[[rank_col]] <- unique(dsub[[rank_col]])
      }
      out # purrr rbinds rows
    })

    preds$.grp <- if (!is.null(group_col) &&
                      !is.null(rank_col) &&
                      isTRUE(facet_by_rank)) {
      # 1) if we have both a group (e.g., big_group) AND a rank facet (e.g.,
      #    genus),  combine them into one id. This prevents lines/ribbons from
      #    crossing between ranks or groups.
      interaction(preds[[group_col]], preds[[rank_col]], drop = TRUE)
    } else if (!is.null(group_col)) {
      # 2) if we only have a group (no rank faceting), group by that.
      preds[[group_col]]
    } else {
      # 3) no grouping column at all --> a single series called "all".
      factor("all")
    }

    p <- ggplot2::ggplot()

    if (!is.null(group_col)) {
      # baseline case with groups
      p <- p +
        ggplot2::geom_point(
          data = df,
          ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]],
                       color = .data[[group_col]]),
          alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0
        )
    } else {
      # no grouping at all (leaves colors aes out)
      p <- p +
        ggplot2::geom_point(
          data = df,
          ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]]),
          alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0,
          colour = "black"
        )
    }

    # use the previously calculated smooths/confidence intervals
    p <- p +
      ggplot2::geom_ribbon(
        data = preds,
        ggplot2::aes(x = .data$.x, ymin = .data$lwr, ymax = .data$upr,
                     group = .data$.grp),
        inherit.aes = FALSE, fill = ci_fill, alpha = ci_alpha, colour = NA
      ) +
      {
        if (!is.null(group_col)) {
          # apply different line colors by group
          ggplot2::geom_line(
            data = preds,
            ggplot2::aes(x = .data$.x, y = .data$.y,
                         color = .data[[group_col]], group = .data$.grp),
            linewidth = 1
          )
        } else {
          # single black line when no grouping
          ggplot2::geom_line(
            data = preds,
            ggplot2::aes(x = .data$.x, y = .data$.y, group = .data$.grp),
            linewidth = 1, colour = "black"
          )
        }
      } +
      ggplot2::labs(x = time_col, y = ylab, color = group_col) +
      theme_phip()

    # apply custom_colors when user provided it
    if (!is.null(group_col)) p <- .add_phip_scales(p, custom_colors)

    # facet by rank where possible/wanted + pass args to facetter
    p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol,
                              facet_scales)
    return(p)
  }

  if (continuous_mode == "binned") {
    # by default 20 bins --> works actually only for trully continuous data,
    # maybe make the number of bins an arg in the future?

    # discretize the range
    rng <- range(df[[time_col]], na.rm = TRUE)
    nbins <- 20
    binwidth <- (rng[2] - rng[1]) / nbins
    make_bins <- function(x, width, xmin) floor((x - xmin) / width)
    df$.bin <- make_bins(df[[time_col]], binwidth, rng[1])
    df$.bin_mid <- rng[1] + (df$.bin + 0.5) * binwidth # calcaulate the midpoint

    # grouping map --> the same as up in the gam section
    group_map <- c(if (!is.null(group_col)) group_col else NULL,
                   if (!is.null(rank_col) && isTRUE(facet_by_rank)) {
                     rank_col
                   } else {
                     NULL
                   },
                   ".bin", ".bin_mid")

    # calculate the mean/sd per bin
    sdat <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_map))) |>
      dplyr::summarise(y = mean(.data[[metric_col]], na.rm = TRUE),
                       sd = stats::sd(.data[[metric_col]], na.rm = TRUE),
                       n = dplyr::n(), .groups = "drop")

    # calculate the se +- lwr/upr CI bands based on the normal distribution
    # (simplified scenario)
    z <- stats::qnorm(1 - (1 - ci_level) / 2)
    sdat <- sdat |>
      dplyr::mutate(se = sd / sqrt(pmax(n, 1)),
                    lwr = y - z * se, upr = y + z * se)

    # the same as in the gam section (see comments up)
    sdat$.grp <- if (!is.null(group_col) &&
                     !is.null(rank_col) &&
                     isTRUE(facet_by_rank)) {
      interaction(sdat[[group_col]], sdat[[rank_col]], drop = TRUE)
    } else if (!is.null(group_col)) {
      sdat[[group_col]]
    } else {
      factor("all")
    }

    p <- ggplot2::ggplot()

    # same logic as up
    if (!is.null(group_col)) {
      p <- p +
        ggplot2::geom_point(
          data = df,
          ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]],
                       color = .data[[group_col]]),
          alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0
        )
    } else {
      p <- p +
        ggplot2::geom_point(
          data = df,
          ggplot2::aes(x = .data[[time_col]], y = .data[[metric_col]]),
          alpha = point_alpha, size = if (point_alpha > 0) 1.6 else 0,
          colour = "black"
        )
    }

    p <- p +
      ggplot2::geom_ribbon(
        data = sdat,
        ggplot2::aes(x = .data$.bin_mid, ymin = .data$lwr, ymax = .data$upr,
                     group = .data$.grp),
        fill = ci_fill, alpha = ci_alpha, colour = NA, inherit.aes = FALSE
      ) +
      {
        if (!is.null(group_col)) {
          ggplot2::geom_line(
            data = sdat,
            ggplot2::aes(x = .data$.bin_mid, y = .data$y,
                         color = .data[[group_col]], group = .data$.grp),
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
    p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol,
                              facet_scales)
    return(p)
  }

  # LOESS fallback for regression/CI; actually not recommended (the best is gam)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[time_col]],
                                        y = .data[[metric_col]])) +
    {
      if (!is.null(group_col)) {
        ggplot2::aes(color = .data[[group_col]])
      } else {
        ggplot2::aes()
      }
    } +
    ggplot2::geom_point(alpha = point_alpha,
                        size = if (point_alpha > 0) 1.6 else 0) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, span = 0.75,
                         level = ci_level) +
    ggplot2::labs(x = time_col, y = ylab, color = group_col) +
    theme_phip()

  if (!is.null(group_col)) p <- .add_phip_scales(p, custom_colors)
  p <- .maybe_facet_by_rank(p, df, rank_col, facet_by_rank, ncol, facet_scales)
  p
}

# ------------------------------------------------------------------------------
# phip_data method
# ------------------------------------------------------------------------------
# i defined a generic earlier and now define 3 methods: for phip_data,
# data.frame or for precomputed phip_alpha_diversity --> different paths for the
# same result
#' @rdname plot_alpha_diversity
#' @export
plot_alpha_diversity.phip_data <- function(
    x,
    metric = c("richness","shannon_diversity","simpson_diversity"),
    group_col = "group",
    rank_col  = "rank",
    filter_groups = NULL,
    filter_ranks  = NULL,
    custom_colors = NULL,
    facet_by_rank = TRUE,
    ncol = 2,
    facet_scales = "fixed",
    # longitudinal
    time_col = NULL,
    continuous_mode = c("gam","binned","loess"),
    gam_k = 7,
    point_alpha = 0.25,
    # CIs
    ci_method = c("model","bootstrap"),
    ci_level  = 0.95,
    boot_R = 500,
    boot_seed = NULL,
    boot_progress = TRUE,
    ci_fill = "grey70",
    ci_alpha = 0.15,
    ...
) {
  stopifnot(inherits(x, "phip_data"))
  continuous_mode <- match.arg(continuous_mode)
  ci_method       <- match.arg(ci_method)

  # logging time
  .ph_with_timing(
    headline = "Plotting alpha diversity (<phip_data>)",
    step     = sprintf("metrics: %s",
                       paste(add_quotes(tolower(metric), 1L), collapse = ", ")),
    expr = {
      # longitudinal mode reserved for the phip_alpha_diversity
      if (!is.null(time_col) && isTRUE(!isTRUE(x$meta$longitudinal))) {
        .ph_abort(
          headline = "Longitudinal plotting not allowed for this object.",
          step     = "meta$longitudinal check",
          bullets  = "x$meta$longitudinal is FALSE"
        )
      }

      # compute, then defer to the precomputed method to plot
      alpha_list <- compute_alpha_diversity(
        x,
        group_col  = group_col,
        ranks      = rank_col,
        carry_cols = time_col
      )

      # generic call
      plot_alpha_diversity(
        alpha_list,
        metric         = metric,
        group_col      = group_col,
        rank_col       = rank_col,
        filter_groups  = filter_groups,
        filter_ranks   = filter_ranks,
        custom_colors  = custom_colors,
        facet_by_rank  = facet_by_rank,
        ncol           = ncol,
        facet_scales   = facet_scales,
        time_col       = time_col,
        continuous_mode = continuous_mode,
        gam_k          = gam_k,
        point_alpha    = point_alpha,
        ci_method      = ci_method,
        ci_level       = ci_level,
        boot_R         = boot_R,
        boot_seed      = boot_seed,
        boot_progress  = boot_progress,
        ci_fill        = ci_fill,
        ci_alpha       = ci_alpha
      )
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ------------------------------------------------------------------------------
# data.frame method
# ------------------------------------------------------------------------------

#' @rdname plot_alpha_diversity
#' @export
plot_alpha_diversity.data.frame <- function(
    x,
    metric = c("richness","shannon_diversity","simpson_diversity"),
    group_col = "group",
    rank_col  = "rank",
    filter_groups = NULL,
    filter_ranks  = NULL,
    custom_colors = NULL,
    facet_by_rank = TRUE,
    ncol = 2,
    facet_scales = "fixed",
    # longitudinal
    time_col = NULL,
    continuous_mode = c("gam","binned","loess"),
    gam_k = 7,
    point_alpha = 0.25,
    # CIs
    ci_method = c("model","bootstrap"),
    ci_level  = 0.95,
    boot_R = 500,
    boot_seed = NULL,
    boot_progress = TRUE,
    ci_fill = "grey70",
    ci_alpha = 0.15,
    ...
) {
  continuous_mode <- match.arg(continuous_mode)
  ci_method       <- match.arg(ci_method)

  # timelog
  .ph_with_timing(
    headline = "Plotting alpha diversity (data.frame)",
    step     = sprintf("metrics: %s",
                       paste(add_quotes(tolower(metric), 1L), collapse = ", ")),
    expr = {
      # firstly compute the diversity, no additional logs needed
      alpha_list <- compute_alpha_diversity(
        x,
        group_col  = group_col,
        ranks      = rank_col,
        carry_cols = time_col
      )

      # generic call --> dispatcher
      plot_alpha_diversity(
        alpha_list,
        metric         = metric,
        group_col      = group_col,
        rank_col       = rank_col,
        filter_groups  = filter_groups,
        filter_ranks   = filter_ranks,
        custom_colors  = custom_colors,
        facet_by_rank  = facet_by_rank,
        ncol           = ncol,
        facet_scales   = facet_scales,
        time_col       = time_col,
        continuous_mode = continuous_mode,
        gam_k          = gam_k,
        point_alpha    = point_alpha,
        ci_method      = ci_method,
        ci_level       = ci_level,
        boot_R         = boot_R,
        boot_seed      = boot_seed,
        boot_progress  = boot_progress,
        ci_fill        = ci_fill,
        ci_alpha       = ci_alpha
      )
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}

# ------------------------------------------------------------------------------
# phip_alpha_diversity method (precomputed)
# ------------------------------------------------------------------------------

#' @rdname plot_alpha_diversity
#' @export
plot_alpha_diversity.phip_alpha_diversity <- function(
    x,
    metric = c("richness","shannon_diversity","simpson_diversity"),
    group_col = "group",
    rank_col  = "rank",
    filter_groups = NULL,
    filter_ranks  = NULL,
    custom_colors = NULL,
    facet_by_rank = TRUE,
    ncol = 2,
    facet_scales = "fixed",
    # longitudinal
    time_col = NULL,
    continuous_mode = c("gam","binned","loess"),
    gam_k = 7,
    point_alpha = 0.25,
    # CIs
    ci_method = c("model","bootstrap"),
    ci_level  = 0.95,
    boot_R = 500,
    boot_seed = NULL,
    boot_progress = TRUE,
    ci_fill = "grey70",
    ci_alpha = 0.15,
    ...
) {
  # all modes possible with this type of x arg
  continuous_mode <- match.arg(continuous_mode)
  ci_method       <- match.arg(ci_method)

  # log time
  .ph_with_timing(
    headline = "Plotting alpha diversity (precomputed)",
    step     = sprintf("metrics: %s",
                       paste(add_quotes(tolower(metric), 1L), collapse = ", ")),
    expr = {
      # choose the correct tibble from the list
      alpha_df <- if (is.list(x) && !is.null(group_col) &&
                       group_col %in% names(x)) {
        tibble::as_tibble(x[[group_col]])
      } else if (is.list(x)) {
        # ---- robust handling for 1-element lists (e.g., group_col=NULL -->
        # $all_samples)
        # keep only elements that are data.frames/tibbles
        inner <- x[vapply(x, function(el) inherits(el, "data.frame"),
                          logical(1))]
        if (length(inner) == 1L) {
          tibble::as_tibble(inner[[1L]])
        } else if (length(inner) >= 1L) {
          suppressWarnings(dplyr::bind_rows(inner))
        } else {
          # last resort: if x[[1]] is a data.frame use it, otherwise coerce x
          if (length(x) && inherits(x[[1L]], "data.frame")) {
            tibble::as_tibble(x[[1L]])
          } else {
            tibble::as_tibble(x)
          }
        }
      } else {
        tibble::as_tibble(x)
      }

      # soft-attrs passthrough
      if (!is.null(attr(x, "group_cols", exact = TRUE)))
        attr(alpha_df, "group_cols") <- attr(x, "group_cols", exact = TRUE)
      if (!is.null(attr(x, "ranks", exact = TRUE)))
        attr(alpha_df, "ranks") <- attr(x, "ranks", exact = TRUE)

      print(alpha_df)
      # validate essentials
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

      # validate grouping column
      if (!is.null(group_col) && !group_col %in% names(alpha_df)) {
        .ph_abort(
          headline = "Grouping column not found in precomputed alpha data.",
          step     = "input validation",
          bullets  = c(
            sprintf("requested group_col: %s", add_quotes(group_col, 1L)),
            sprintf("available columns: %s",
                    paste(add_quotes(colnames(alpha_df), 1L), collapse = ", "))
          )
        )
      }

      # validating the rank
      if (!is.null(rank_col) && !rank_col %in% names(alpha_df)) {
        .ph_warn(
          headline = "Rank column not found in precomputed alpha data;
          disabling faceting.",
          step     = "input validation",
          bullets  = sprintf("rank_col: %s", add_quotes(rank_col, 1L))
        )
        rank_col <- NULL
      }

      # optional filters
      if (!is.null(group_col) && !is.null(filter_groups)) {
        alpha_df <- dplyr::filter(alpha_df,
                                  .data[[group_col]] %in% !!filter_groups)
      }
      if (!is.null(rank_col) && !is.null(filter_ranks)) {
        alpha_df <- dplyr::filter(alpha_df,
                                  .data[[rank_col]] %in% !!filter_ranks)
      }

      # make one plot per metric
      plots <- lapply(metrics_norm, function(m) {
        .build_alpha_plot(
          alpha_df      = alpha_df,
          metric        = m,
          group_col     = group_col,
          rank_col      = rank_col,
          filter_groups = filter_groups,
          filter_ranks  = filter_ranks,
          custom_colors = custom_colors,
          facet_by_rank = facet_by_rank,
          ncol          = ncol,
          facet_scales  = facet_scales,
          time_col      = if (!is.null(time_col) &&
                              time_col %in% names(alpha_df)) time_col else NULL,
          continuous_mode = continuous_mode,
          gam_k         = gam_k,
          point_alpha   = point_alpha,
          ci_method     = ci_method,
          ci_level      = ci_level,
          boot_R        = boot_R,
          boot_seed     = boot_seed,
          boot_progress = boot_progress,
          ci_fill       = ci_fill,
          ci_alpha      = ci_alpha
        )
      })
      names(plots) <- metrics_norm
      if (length(plots) == 1L) plots[[1L]] else plots
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}
