#' Scatterplots (percent1 vs percent2) from ph_prevalence_compare()
#'
#' @param comparison_tbl tibble or DuckDB lazy tbl from ph_prevalence_compare()
#' @param ranks character vector of ranks to plot (default: all present)
#' @param facet_pairs logical; TRUE (default) facets all (group1,group2) pairs per rank,
#'   FALSE returns multiple plots per pair
#' @param facet_nrow,facet_ncol integers controlling facet_wrap layout (optional)
#' @param sparse_rank character(1) rank name whose "not significant" peptides you want to downsample
#' @param sparse_drop_pct numeric in [0,1]; fraction of "not significant" rows to drop
#'   within sparse_rank (default 0 = keep all). Always keeps all significant rows.
#' @param sparse_seed integer or NULL; set to make downsampling reproducible
#' @param significant_colors named vector for category colors
#' @param interactive TRUE (default) for plotly, FALSE for static ggplot
#'
#' @return A named list of plots. If facet_pairs=TRUE: one plot per rank.
#'         If facet_pairs=FALSE: for each rank a sub-list of plots per (group1 vs group2).
#' @export
ph_scatterplot <- function(
    comparison_tbl,
    ranks         = NULL,
    facet_pairs   = TRUE,
    facet_nrow    = NULL,
    facet_ncol    = NULL,
    sparse_rank   = NULL,
    sparse_drop_pct = 0,
    sparse_seed   = NULL,
    significant_colors = c(
      "not significant"                 = "#386cb0",
      "significant prior correction"    = "#1b9e77",
      "significant post FDR correction" = "#e31a1c"
    ),
    interactive   = TRUE
) {
  .ph_with_timing(
    headline = "make_interactive_scatterplot",
    expr = {
      need <- c("rank","feature","group_col","group1","group2","percent1","percent2","category")
      have <- tryCatch(colnames(comparison_tbl), error = function(...) character())
      if (length(have) == 0 && inherits(comparison_tbl, "tbl_sql")) have <- colnames(comparison_tbl)
      miss <- setdiff(need, have)
      if (length(miss)) .ph_abort("Comparison table is missing required columns", bullets = paste("-", miss))

      keep_cols <- c("rank","feature","group_col","group1","group2","percent1","percent2","category")
      if (!is.null(ranks) && length(ranks)) comparison_tbl <- dplyr::filter(comparison_tbl, rank %in% ranks)
      df <- comparison_tbl |> dplyr::select(dplyr::all_of(keep_cols)) |> dplyr::collect()
      if (!nrow(df)) { .ph_warn("No rows after filtering; returning empty list."); return(invisible(list())) }

      df <- df |>
        dplyr::mutate(
          rank      = as.character(rank),
          group_col = as.character(group_col),
          group1    = as.character(group1),
          group2    = as.character(group2),
          percent1  = as.numeric(percent1),
          percent2  = as.numeric(percent2),
          category  = as.character(category),
          pair_label = paste0(group1, " vs ", group2)
        )

      .downsample_uniform_2d <- function(dat, pool_idx, xcol, ycol, keep_n, seed = NULL) {
        if (keep_n >= length(pool_idx) || keep_n <= 0L) return(pool_idx)
        if (!is.null(seed)) set.seed(as.integer(seed))
        x <- dat[[xcol]][pool_idx]; y <- dat[[ycol]][pool_idx]
        G <- max(1L, ceiling(sqrt(keep_n)))
        xr <- range(x, finite = TRUE); if (diff(xr) == 0) xr <- xr + c(-0.5, 0.5)
        yr <- range(y, finite = TRUE); if (diff(yr) == 0) yr <- yr + c(-0.5, 0.5)
        bx <- pmin(G, pmax(1L, as.integer(floor((x - xr[1]) / diff(xr) * G) + 1L)))
        by <- pmin(G, pmax(1L, as.integer(floor((y - yr[1]) / diff(yr) * G) + 1L)))
        bins <- split(seq_along(pool_idx), paste(bx, by, sep = "_"))
        order_bins <- sample(names(bins), length(bins))
        keep <- integer(0)
        while (length(keep) < keep_n && length(order_bins) > 0) {
          for (b in order_bins) {
            idxs <- bins[[b]]
            if (length(idxs)) {
              pick <- sample(idxs, 1L); keep <- c(keep, pick)
              bins[[b]] <- setdiff(bins[[b]], pick)
              if (length(keep) >= keep_n) break
            }
          }
          order_bins <- Filter(function(b) length(bins[[b]]) > 0, order_bins)
        }
        pool_idx[keep]
      }

      if (!is.null(sparse_rank) && sparse_drop_pct > 0) {
        if (!sparse_rank %in% unique(df$rank)) {
          .ph_warn(paste0("sparse_rank '", sparse_rank, "' not found in data; skipping downsampling"))
        } else {
          pool_idx <- which(df$rank == sparse_rank & df$category == "not significant")
          if (length(pool_idx) > 0) {
            n_drop <- floor(sparse_drop_pct * length(pool_idx))
            keep_n <- length(pool_idx) - n_drop
            if (keep_n > 0L) {
              keep_pool <- .downsample_uniform_2d(df, pool_idx, "percent1", "percent2", keep_n, sparse_seed)
              keep_mask <- rep(TRUE, nrow(df)); keep_mask[setdiff(pool_idx, keep_pool)] <- FALSE
              df <- df[keep_mask, , drop = FALSE]
              .ph_log_info("Downsampled 'not significant' rows (uniform 2D coverage)",
                           bullets = c(
                             paste0("rank: ", sparse_rank),
                             paste0("pool: ", length(pool_idx)),
                             paste0("kept (ns): ", keep_n),
                             paste0("dropped (ns): ", n_drop)
                           ))
            }
          }
        }
      }

      df <- df |>
        dplyr::mutate(
          tooltip_txt = paste0(
            "Feature: ", feature, "<br>",
            "Rank: ",   rank,    "<br>",
            group1, " %: ", signif(percent1, 3), "<br>",
            group2, " %: ", signif(percent2, 3), "<br>",
            "Pair: ", pair_label, "<br>",
            "Category: ", category
          )
        )

      add_color_scale <- function() {
        if (is.null(significant_colors)) {
          scale_color_phip(name = NULL)
        } else {
          ggplot2::scale_color_manual(values = significant_colors, name = NULL)
        }
      }

      # fixed limits & ticks
      lims   <- c(-10, 110)
      breaks <- c(0, 25, 50, 75, 100)

      build_rank_plot <- function(dat_rank) {
        if (facet_pairs) {
          xlab <- "% with presence (group1 varies by panel)"
          ylab <- "% with presence (group2 varies by panel)"
        } else {
          xlab <- paste0("% ", unique(dat_rank$group1)[1], " with presence")
          ylab <- paste0("% ", unique(dat_rank$group2)[1], " with presence")
        }

        base <- ggplot2::ggplot(dat_rank, ggplot2::aes(x = percent1, y = percent2)) +
          {
            if (interactive) {
              ggplot2::geom_point(ggplot2::aes(color = category, text = tooltip_txt), alpha = 0.7)
            } else {
              ggplot2::geom_point(ggplot2::aes(color = category), alpha = 0.7)
            }
          } +
          add_color_scale() +
          # NOTE: ticks without hard scale limits
          ggplot2::scale_x_continuous(breaks = breaks, labels = breaks, expand = c(0, 0)) +
          ggplot2::scale_y_continuous(breaks = breaks, labels = breaks, expand = c(0, 0)) +
          # visible window set here (prevents shaving on outer facets)
          ggplot2::coord_cartesian(xlim = lims, ylim = lims) +
          ggplot2::labs(
            x = xlab, y = ylab,
            title = paste0("Rank: ", unique(dat_rank$rank)[1])
          ) +
          theme_phip(base_size = 12) +
          ggplot2::theme(
            legend.position = "none",
            # a touch more breathing room
            plot.margin     = grid::unit(c(10, 8, 10, 8), "pt"),
            panel.spacing   = grid::unit(7, "pt")
          )

        if (facet_pairs) {
          base <- base + ggplot2::facet_wrap(
            ~ pair_label,
            scales = "fixed",
            nrow = facet_nrow, ncol = facet_ncol
          )

          if (interactive) {
            p <- plotly::ggplotly(base, tooltip = "text")

            # make sure plotly respects the same view; prevent trace clipping
            p <- p |>
              plotly::style(cliponaxis = FALSE) |>
              plotly::layout(
                showlegend = FALSE,
                xaxis = list(range = lims, scaleanchor = "y", scaleratio = 1, automargin = FALSE),
                yaxis = list(range = lims, scaleanchor = "x", scaleratio = 1, automargin = FALSE),
                margin = list(l = 60, r = 60, b = 60, t = 60)
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
              plt <- plotly::ggplotly(p, tooltip = "text") |>
                plotly::style(cliponaxis = FALSE) |>
                plotly::layout(
                  showlegend = FALSE,
                  xaxis = list(range = lims, scaleanchor = "y", scaleratio = 1, automargin = TRUE),
                  yaxis = list(range = lims, scaleanchor = "x", scaleratio = 1, automargin = TRUE),
                  margin = list(l = 86, r = 60, b = 86, t = 80)
                )
              plt
            } else p
          })
          names(out) <- names(pieces)
          out
        }
      }

      .ph_log_info("Building plots per rank")
      by_rank <- split(df, df$rank, drop = TRUE)
      plots <- lapply(names(by_rank), function(rk) {
        dat <- by_rank[[rk]]
        if (!nrow(dat)) { .ph_warn(paste0("No rows for rank '", rk, "'")); return(NULL) }
        build_rank_plot(dat)
      })
      names(plots) <- names(by_rank)

      .ph_log_ok("Plots built", bullets = c(
        paste0("ranks: ", paste(names(plots), collapse = ", ")),
        paste0("facet_pairs: ", facet_pairs),
        if (!is.null(sparse_rank)) paste0("sparse_rank: ", sparse_rank, " (drop=", sparse_drop_pct, ")") else NULL
      )[!vapply(c(
        paste0("ranks: ", paste(names(plots), collapse = ", ")),
        paste0("facet_pairs: ", facet_pairs),
        if (!is.null(sparse_rank)) paste0("sparse_rank: ", sparse_rank, " (drop=", sparse_drop_pct, ")") else NA_character_
      ), is.na, logical(1))])

      plots
    }
  )
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
    ranks           = NULL,
    facet_pairs     = TRUE,
    facet_nrow      = NULL,
    facet_ncol      = NULL,
    sparse_rank     = NULL,
    sparse_drop_pct = 0,
    sparse_seed     = NULL,
    fc_cut          = 1,
    p_cut           = 0.05,
    significant_colors = c(
      "not significant"                 = "#386cb0",
      "significant prior correction"    = "#1b9e77",
      "significant post FDR correction" = "#e31a1c"
    ),
    interactive     = TRUE
) {
  .ph_with_timing(
    headline = "make_interactive_volcano",
    expr = {
      # required + optional
      required_cols <- c("rank","feature","group_col","group1","group2","category")
      optional_cols <- c("ratio","prop1","prop2","p_raw","p_adj","pvals_not_adj")

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
          rank      = as.character(rank),
          group_col = as.character(group_col),
          group1    = as.character(group1),
          group2    = as.character(group2),
          category  = as.character(category),
          pair_label = paste0(group1, " vs ", group2)
        )

      # ratio fallback
      if (!"ratio" %in% colnames(df)) {
        if (!all(c("prop1","prop2") %in% colnames(df))) {
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
      p_col <- intersect(c("p_raw","pvals_not_adj","p_adj"), colnames(df))[1]
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
                           ))
            }
          }
        }
      }

      # tooltips
      df <- df |>
        dplyr::mutate(
          tooltip_txt = paste0(
            "Feature: ", feature, "<br>",
            "Rank: ",   rank,    "<br>",
            "log2(ratio): ", signif(log2ratio, 3), "<br>",
            "-log10(p): ",     signif(nlog10p, 3), "<br>",
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
          base <- base + ggplot2::facet_wrap(~ pair_label, scales = "fixed",
                                             nrow = facet_nrow, ncol = facet_ncol)
          if (interactive) {
            p <- plotly::ggplotly(base, tooltip = c("text","x","y"))
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
              plt <- plotly::ggplotly(p, tooltip = c("text","x","y"))
              plt <- plt |>
                plotly::layout(
                  showlegend = FALSE,
                  xaxis = list(automargin = TRUE),
                  yaxis = list(automargin = TRUE, title = list(standoff = 10))
                )
              plt
            } else p
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
