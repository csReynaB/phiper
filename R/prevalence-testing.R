# ==============================================================================
# Prevalence (counts/percent) + pairwise tests across multiple group vars & ranks
# DuckDB-first; optional materialized TABLE return; auto-parallel via future plan
# ==============================================================================

#' @title Prevalence by group with pairwise tests (multi-rank; phiper logging)
#' @description
#' Compute per-group counts/percentages for one or more ranks and build all
#' pairwise comparisons within each grouping variable (or a combined interaction).
#' All joins/counts/ratios are executed in DuckDB. P-values are computed in R
#' (chunked; auto-parallel if a `future` plan is active). BH/FDR is computed
#' inside DuckDB per `group_col`.
#'
#' Testing is **adaptive**:
#' - If all four expected counts in the 2x2 are ≥ 5, use chi-square with Yates
#'   continuity correction (`stats::prop.test(..., correct = TRUE)`).
#' - Otherwise use Fisher’s exact test.
#'
#' @param x phip_data or long data.frame with sample_id, peptide_id, `exist_col`, `group_cols`
#' @param rank_cols character vector of ranks (e.g., c("peptide_id","species"))
#' @param group_cols character vector of grouping columns
#' @param exist_col presence indicator (default "exist")
#' @param prevalence_threshold percent threshold for pairwise keep (default 0)
#' @param p_adjust kept for future extensibility (current: BH in DuckDB)
#' @param parallel logical or NULL. NULL = auto (uses future plan if available)
#' @param compute_ratios_db logical; compute ratio/delta_ratio in DB (default TRUE)
#' @param interaction logical; if TRUE combine two grouping cols and test only that
#' @param combine_cols length-2 character vector from `group_cols` when interaction=TRUE
#' @param interaction_sep separator between combined levels (default "::")
#' @param collect logical; if TRUE return a collected tibble; else a DuckDB TABLE (lazy tbl)
#' @param register_name optional TABLE name when collect=FALSE; if NULL timestamped.
#'
#' @return collected tibble or materialized DuckDB TABLE (lazy `tbl`).
#' @export
ph_prevalence_compare <- function(x,
                                  rank_cols,
                                  group_cols,
                                  exist_col = "exist",
                                  prevalence_threshold = 0,
                                  p_adjust = "BH",
                                  parallel = NULL,
                                  compute_ratios_db = TRUE,
                                  interaction = FALSE,
                                  combine_cols = NULL,
                                  interaction_sep = "::",
                                  collect = TRUE,
                                  register_name = NULL) {
  .ph_with_timing(
    headline = "prevalence_compare",
    step = NULL,
    bullets = NULL,
    expr = {
      .q <- function(con, nm) as.character(DBI::dbQuoteIdentifier(con, nm))
      .nowstamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
      .sym <- rlang::sym
      .syms <- rlang::syms
      chunk_n <- getOption("phiper.prev.chunk", 1e6L) # was 250k

      # ---- arg checks ----------------------------------------------------------
      tryCatch(
        {
          chk::chk_character(rank_cols)
          chk::chk_true(length(rank_cols) >= 1)
          chk::chk_character(group_cols)
          chk::chk_true(length(group_cols) >= 1)
          chk::chk_string(exist_col)
          chk::chk_number(prevalence_threshold)
          chk::chk_true(prevalence_threshold >= 0)
          if (isTRUE(interaction)) {
            if (is.null(combine_cols)) combine_cols <- utils::head(group_cols, 2L)
            if (length(combine_cols) != 2L) {
              .ph_abort("When `interaction = TRUE`, `combine_cols` must be length 2.")
            }
            miss_int <- setdiff(combine_cols, group_cols)
            if (length(miss_int)) {
              .ph_abort("`combine_cols` must be subset of `group_cols`.", bullets = paste("-", miss_int))
            }
          }
        },
        error = function(e) .ph_abort("Invalid arguments", bullets = e$message)
      )

      .ph_log_info("Preparing input data", bullets = c(
        paste0("ranks: ", paste(rank_cols, collapse = ", ")),
        paste0("group_cols: ", paste(group_cols, collapse = ", ")),
        if (interaction) paste0("interaction on: ", paste(combine_cols, collapse = " + ")) else NULL,
        paste0("exist_col: ", exist_col),
        paste0("threshold(%): ", prevalence_threshold),
        paste0("collect: ", collect)
      )[!vapply(c(
        paste0("ranks: ", paste(rank_cols, collapse = ", ")),
        paste0("group_cols: ", paste(group_cols, collapse = ", ")),
        if (interaction) paste0("interaction on: ", paste(combine_cols, collapse = " + ")) else NA_character_,
        paste0("exist_col: ", exist_col),
        paste0("threshold(%): ", prevalence_threshold),
        paste0("collect: ", collect)
      ), is.na, logical(1))])

      # ---- normalize to long df -----------------------------------------------
      df_long <- try(
        {
          if (inherits(x, "phip_data")) {
            x$data_long |>
              dplyr::select(tidyselect::any_of(c("sample_id", "peptide_id", exist_col, group_cols)))
          } else {
            chk::chk_data(x)
            need <- c("sample_id", "peptide_id", exist_col, group_cols)
            miss <- setdiff(need, colnames(x))
            if (length(miss)) .ph_abort("Missing required columns", bullets = paste("-", miss))
            tibble::as_tibble(x) |> dplyr::select(tidyselect::any_of(need))
          }
        },
        silent = TRUE
      )
      if (inherits(df_long, "try-error")) .ph_abort("Could not prepare input data.")

      df_long <- df_long |>
        dplyr::mutate(!!.sym(exist_col) := dplyr::coalesce(as.integer(!!.sym(exist_col)), 0L))

      # connection
      con <- NULL
      if (inherits(x, "phip_data")) {
        con <- tryCatch(x$meta$con, error = function(...) NULL)
      }
      if (is.null(con) && inherits(df_long, "tbl_sql")) con <- dbplyr::remote_con(df_long)
      if (is.null(con)) .ph_abort("No DuckDB connection found. Use <phip_data>$meta$con or a DuckDB-backed tbl_sql.")

      view_const <- if (inherits(x, "phip_data")) {
        attr(x, "view") %||% (x$meta$view %||% NA_character_)
      } else {
        NA_character_
      }

      # ---- ranks & library -----------------------------------------------------
      ranks_needing_lib <- setdiff(rank_cols, "peptide_id")
      if (length(ranks_needing_lib)) {
        if (!inherits(x, "phip_data") || is.null(x$peptide_library)) {
          .ph_abort("Peptide library required for non-peptide ranks.")
        }
        miss_tax <- setdiff(ranks_needing_lib, colnames(x$peptide_library))
        if (length(miss_tax)) .ph_abort("Requested taxonomy columns not in peptide_library.", bullets = paste("-", miss_tax))
      }

      df_ranked <- df_long
      if (length(ranks_needing_lib)) {
        lib_cols <- c("peptide_id", ranks_needing_lib)
        lib_min <- x$peptide_library |>
          dplyr::select(tidyselect::all_of(lib_cols)) |>
          dplyr::distinct()
        df_ranked <- df_ranked |>
          dplyr::left_join(lib_min, by = "peptide_id", copy = TRUE)
      }

      available_ranks <- intersect(rank_cols, colnames(df_ranked))
      if (!length(available_ranks)) .ph_abort("None of the requested rank_cols are available.")
      .ph_log_info("Ranks resolved", bullets = paste("- available:", paste(available_ranks, collapse = ", ")))

      df_ranked_long <- df_ranked |>
        tidyr::pivot_longer(
          cols = tidyselect::all_of(available_ranks),
          names_to = "rank", values_to = "feature"
        )

      # ---- grouping mode -------------------------------------------------------
      .ph_log_info(if (interaction) "Using interaction-only grouping" else "Using per-column grouping")

      build_group_view <- function(tbl) {
        if (!interaction) {
          tbl |>
            tidyr::pivot_longer(
              cols = tidyselect::all_of(group_cols),
              names_to = "group_col", values_to = "group_value"
            ) |>
            dplyr::filter(!is.na(group_value))
        } else {
          c1 <- combine_cols[1]
          c2 <- combine_cols[2]
          comb_name <- paste(combine_cols, collapse = " + ")
          tbl |>
            dplyr::filter(!is.na(!!.sym(c1)) & !is.na(!!.sym(c2))) |>
            dplyr::mutate(
              group_col   = comb_name,
              group_value = paste0(!!.sym(c1), interaction_sep, !!.sym(c2))
            )
        }
      }

      # ---- denominators & counts (lazy/DB) ------------------------------------
      .ph_log_info("Computing cohort sizes (N) per grouping variable")
      gs_view <- build_group_view(df_ranked_long)

      group_sizes <- gs_view |>
        dplyr::distinct(group_col, group_value, sample_id) |>
        dplyr::count(group_col, group_value, name = "N")

      gs_n <- group_sizes |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::collect()
      if (gs_n$n == 0L) .ph_abort("No non-missing group values; cannot compute denominators.")

      .ph_log_info("Counting present samples per feature (vectorised)")
      present_counts <- gs_view |>
        dplyr::filter(!!.sym(exist_col) > 0) |>
        dplyr::distinct(group_col, group_value, rank, feature, sample_id) |>
        dplyr::count(group_col, group_value, rank, feature, name = "n_present")

      features_per_rank <- present_counts |>
        dplyr::distinct(rank, feature)

      fpr_n <- features_per_rank |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::collect()
      if (fpr_n$n == 0L) {
        .ph_warn("No features present in any cohort; returning empty")
        return(if (collect) tibble::tibble()[0, ] else dplyr::tbl(con, dplyr::sql("SELECT 1 WHERE 0=1")))
      }

      base_grid <- group_sizes |>
        dplyr::select(group_col, group_value) |>
        dplyr::distinct() |>
        dplyr::mutate(.dummy = 1L) |>
        dplyr::inner_join(
          features_per_rank |> dplyr::mutate(.dummy = 1L),
          by = ".dummy"
        ) |>
        dplyr::select(group_col, group_value, rank, feature)

      stats_long <- base_grid |>
        dplyr::left_join(present_counts,
          by = c("group_col", "group_value", "rank", "feature")
        ) |>
        dplyr::left_join(group_sizes,
          by = c("group_col", "group_value")
        ) |>
        dplyr::mutate(
          n_present = dplyr::coalesce(n_present, 0L),
          prop      = dplyr::if_else(N > 0, n_present / N, NA_real_),
          percent   = 100 * prop
        )

      if (!is.null(view_const) && !is.na(view_const)) {
        stats_long <- stats_long |> dplyr::mutate(view = !!view_const, .before = 1)
      }

      # ---- pairwise (lazy/DB) -------------------------------------------------
      .ph_log_info("Building all pairwise comparisons within each group_col")
      by_cols <- intersect(c("view", "rank", "feature", "group_col"), colnames(stats_long))
      has_view <- "view" %in% colnames(stats_long)

      pairs_joined <- stats_long |>
        dplyr::inner_join(
          stats_long,
          by = by_cols,
          suffix = c("_1", "_2")
        ) |>
        dplyr::filter(group_value_1 != group_value_2) |>
        dplyr::mutate(
          group1_chr = dplyr::if_else(group_value_1 <= group_value_2, group_value_1, group_value_2),
          group2_chr = dplyr::if_else(group_value_1 <= group_value_2, group_value_2, group_value_1),
          n1_val = dplyr::if_else(group_value_1 <= group_value_2, n_present_1, n_present_2),
          n1_tot = dplyr::if_else(group_value_1 <= group_value_2, N_1, N_2),
          p1_val = dplyr::if_else(group_value_1 <= group_value_2, prop_1, prop_2),
          pct1val = dplyr::if_else(group_value_1 <= group_value_2, percent_1, percent_2),
          n2_val = dplyr::if_else(group_value_1 <= group_value_2, n_present_2, n_present_1),
          n2_tot = dplyr::if_else(group_value_1 <= group_value_2, N_2, N_1),
          p2_val = dplyr::if_else(group_value_1 <= group_value_2, prop_2, prop_1),
          pct2val = dplyr::if_else(group_value_1 <= group_value_2, percent_2, percent_1)
        ) |>
        dplyr::filter(group_value_1 == group1_chr)

      pairs_lazy <- (if (has_view) {
        pairs_joined |>
          dplyr::transmute(
            view, rank, feature, group_col,
            group1 = group1_chr,
            n1 = n1_val,
            n1_tot = n1_tot,
            prop1 = p1_val,
            percent1 = pct1val,
            group2 = group2_chr,
            n2 = n2_val,
            n2_tot = n2_tot,
            prop2 = p2_val,
            percent2 = pct2val
          )
      } else {
        pairs_joined |>
          dplyr::transmute(
            rank, feature, group_col,
            group1 = group1_chr,
            n1 = n1_val,
            n1_tot = n1_tot,
            prop1 = p1_val,
            percent1 = pct1val,
            group2 = group2_chr,
            n2 = n2_val,
            n2_tot = n2_tot,
            prop2 = p2_val,
            percent2 = pct2val
          )
      })

      if (!isTRUE(prevalence_threshold <= 0)) {
        pairs_lazy <- pairs_lazy |>
          dplyr::filter(percent1 >= prevalence_threshold | percent2 >= prevalence_threshold)
      }

      if (isTRUE(compute_ratios_db)) {
        pairs_lazy <- pairs_lazy |>
          dplyr::mutate(
            prop1_eps = dplyr::if_else(n1 == 0L, (n1 + 0.5) / n1_tot, prop1),
            prop2_eps = dplyr::if_else(n2 == 0L, (n2 + 0.5) / n2_tot, prop2),
            ratio = prop1_eps / prop2_eps,
            d1 = dplyr::if_else(n1 == 0L, (n1 + 1.0) / n1_tot, prop1),
            d2 = dplyr::if_else(n2 == 0L, (n2 + 1.0) / n2_tot, prop2),
            delta_ratio = dplyr::if_else(d1 >= d2, d1 / d2 - 1, -(d2 / d1 - 1))
          )
      }

      pairs_lazy <- pairs_lazy |>
        dbplyr::window_order(group_col, rank, feature, group1, group2) |>
        dplyr::mutate(row_id = dplyr::row_number()) |>
        dplyr::filter(n1_tot > 0, n2_tot > 0) # keep zero counts; only denominators must be >0

      # ---- materialize TABLE ---------------------------------------------------
      if (is.null(register_name)) register_name <- paste0("ph_prev_", .nowstamp())
      tbl_q <- .q(con, register_name)

      sql_body <- dbplyr::sql_render(pairs_lazy, con)
      DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", tbl_q))
      DBI::dbExecute(con, paste0("CREATE TABLE ", tbl_q, " AS ", sql_body))

      .ph_log_ok("Materialized DuckDB TABLE",
        bullets = c(
          paste0("name: ", register_name),
          "computing p-values (auto: Fisher vs chisq+Yates) + BH in-table"
        )
      )

      # ---- P-values (chunked; auto-parallel if future plan) --------------------
      p_core <- function(n1, N1, n2, N2) {
        n1 <- as.double(n1)
        N1 <- as.double(N1)
        n2 <- as.double(n2)
        N2 <- as.double(N2)
        if (any(!is.finite(c(n1, N1, n2, N2)))) {
          return(NA_real_)
        }
        if (N1 <= 0 || N2 <= 0) {
          return(NA_real_)
        }

        a <- max(0, min(n1, N1))
        b <- max(0, min(n2, N2))
        c1 <- N1 - a
        c2 <- N2 - b

        tot <- N1 + N2
        pos <- a + b
        if (tot <= 0) {
          return(NA_real_)
        }
        exp_a <- N1 * pos / tot
        exp_b <- N2 * pos / tot
        exp_c1 <- N1 - exp_a
        exp_c2 <- N2 - exp_b

        safe_chisq <- isTRUE(all(c(exp_a, exp_b, exp_c1, exp_c2) >= 5))

        out <- if (safe_chisq) {
          try(stats::prop.test(x = c(a, b), n = c(N1, N2), correct = TRUE)$p.value, silent = TRUE)
        } else {
          try(stats::fisher.test(matrix(c(a, c1, b, c2), 2, byrow = TRUE))$p.value, silent = TRUE)
        }
        if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
      }

      DBI::dbExecute(con, paste0("ALTER TABLE ", tbl_q, " ADD COLUMN IF NOT EXISTS p_raw DOUBLE"))
      DBI::dbExecute(con, paste0("ALTER TABLE ", tbl_q, " ADD COLUMN IF NOT EXISTS p_adj DOUBLE"))
      DBI::dbExecute(con, paste0("ALTER TABLE ", tbl_q, " ADD COLUMN IF NOT EXISTS passed_raw BOOLEAN"))
      DBI::dbExecute(con, paste0("ALTER TABLE ", tbl_q, " ADD COLUMN IF NOT EXISTS passed_adj BOOLEAN"))
      DBI::dbExecute(con, paste0("ALTER TABLE ", tbl_q, " ADD COLUMN IF NOT EXISTS category VARCHAR"))

      tmp_q <- .q(con, paste0(register_name, "_p_tmp_", format(Sys.time(), "%H%M%S")))
      DBI::dbExecute(con, paste0("CREATE TEMP TABLE ", tmp_q, " (row_id BIGINT, p_raw DOUBLE)"))

      rs <- DBI::dbSendQuery(con, paste0(
        "SELECT row_id, n1, n1_tot, n2, n2_tot FROM ", tbl_q, " ORDER BY row_id"
      ))
      on.exit(try(DBI::dbClearResult(rs), silent = TRUE), add = TRUE)

      want_parallel <- if (is.null(parallel)) {
        requireNamespace("future", quietly = TRUE) &&
          requireNamespace("future.apply", quietly = TRUE) &&
          tryCatch(future::nbrOfWorkers() > 1L, error = function(...) FALSE)
      } else {
        isTRUE(parallel)
      }
      if (want_parallel &&
        (!requireNamespace("future", quietly = TRUE) ||
          !requireNamespace("future.apply", quietly = TRUE))) {
        .ph_abort("Parallel requested but {future}/{future.apply} not available.")
      }

      repeat {
        chunk <- DBI::dbFetch(rs, n = chunk_n)
        if (!nrow(chunk)) break

        colnames(chunk) <- tolower(colnames(chunk))
        req <- c("row_id", "n1", "n1_tot", "n2", "n2_tot")
        miss <- setdiff(req, colnames(chunk))
        if (length(miss)) .ph_abort("Fetched chunk is missing required columns.", bullets = paste("-", miss))

        chunk$row_id <- as.numeric(chunk$row_id)
        chunk$n1 <- as.double(chunk$n1)
        chunk$n1_tot <- as.double(chunk$n1_tot)
        chunk$n2 <- as.double(chunk$n2)
        chunk$n2_tot <- as.double(chunk$n2_tot)

        n <- nrow(chunk)
        workers <- if (want_parallel) tryCatch(future::nbrOfWorkers(), error = function(e) 1L) else 1L

        # extract plain vectors (small, easy to serialize once)
        n1v <- chunk$n1
        N1v <- chunk$n1_tot
        n2v <- chunk$n2
        N2v <- chunk$n2_tot

        if (want_parallel && workers > 1L) {
          # ~one block per worker; adjust breaks if you want more/less granularity
          splits <- split(seq_len(n), cut(seq_len(n), breaks = workers, labels = FALSE))

          parts <- future.apply::future_lapply(
            splits,
            FUN = function(ii, n1v, N1v, n2v, N2v) {
              # define a tiny, self-contained helper INSIDE the worker
              p_core <- function(n1, N1, n2, N2) {
                n1 <- as.double(n1)
                N1 <- as.double(N1)
                n2 <- as.double(n2)
                N2 <- as.double(N2)
                if (any(!is.finite(c(n1, N1, n2, N2)))) {
                  return(NA_real_)
                }
                if (N1 <= 0 || N2 <= 0) {
                  return(NA_real_)
                }
                a <- max(0, min(n1, N1))
                b <- max(0, min(n2, N2))
                c1 <- N1 - a
                c2 <- N2 - b
                tot <- N1 + N2
                pos <- a + b
                if (tot <= 0) {
                  return(NA_real_)
                }
                exp_a <- N1 * pos / tot
                exp_b <- N2 * pos / tot
                exp_c1 <- N1 - exp_a
                exp_c2 <- N2 - exp_b
                safe_chisq <- isTRUE(all(c(exp_a, exp_b, exp_c1, exp_c2) >= 5))
                out <- if (safe_chisq) {
                  try(stats::prop.test(x = c(a, b), n = c(N1, N2), correct = TRUE)$p.value, silent = TRUE)
                } else {
                  try(stats::fisher.test(matrix(c(a, c1, b, c2), 2, byrow = TRUE))$p.value, silent = TRUE)
                }
                if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
              }
              vapply(ii, function(i) p_core(n1v[i], N1v[i], n2v[i], N2v[i]), numeric(1))
            },
            n1v = n1v, N1v = N1v, n2v = n2v, N2v = N2v,
            future.seed = TRUE,
            future.scheduling = 1,
            future.globals = FALSE # <- do NOT ship parent envs
          )
          pvals <- unlist(parts, use.names = FALSE)
        } else {
          # sequential fallback (same helper, local)
          p_core <- function(n1, N1, n2, N2) {
            n1 <- as.double(n1)
            N1 <- as.double(N1)
            n2 <- as.double(n2)
            N2 <- as.double(N2)
            if (any(!is.finite(c(n1, N1, n2, N2)))) {
              return(NA_real_)
            }
            if (N1 <= 0 || N2 <= 0) {
              return(NA_real_)
            }
            a <- max(0, min(n1, N1))
            b <- max(0, min(n2, N2))
            c1 <- N1 - a
            c2 <- N2 - b
            tot <- N1 + N2
            pos <- a + b
            if (tot <= 0) {
              return(NA_real_)
            }
            exp_a <- N1 * pos / tot
            exp_b <- N2 * pos / tot
            exp_c1 <- N1 - exp_a
            exp_c2 <- N2 - exp_b
            safe_chisq <- isTRUE(all(c(exp_a, exp_b, exp_c1, exp_c2) >= 5))
            out <- if (safe_chisq) {
              try(stats::prop.test(x = c(a, b), n = c(N1, N2), correct = TRUE)$p.value, silent = TRUE)
            } else {
              try(stats::fisher.test(matrix(c(a, c1, b, c2), 2, byrow = TRUE))$p.value, silent = TRUE)
            }
            if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
          }
          pvals <- vapply(seq_len(n), function(i) p_core(n1v[i], N1v[i], n2v[i], N2v[i]), numeric(1))
        }

        DBI::dbAppendTable(
          con,
          name  = gsub('^\"|\"$', "", tmp_q),
          value = data.frame(row_id = chunk$row_id, p_raw = as.numeric(pvals))
        )
      }

      DBI::dbExecute(con, paste0(
        "UPDATE ", tbl_q, " t SET p_raw = s.p_raw FROM ", tmp_q, " s WHERE t.row_id = s.row_id"
      ))
      DBI::dbExecute(con, paste0("DROP TABLE ", tmp_q))

      # ---- BH/FDR in DuckDB (monotone, per group_col) --------------------------
      DBI::dbExecute(con, paste0("
      UPDATE ", tbl_q, " t
         SET p_adj = a.p_adj,
             passed_raw = (t.p_raw < 0.05),
             passed_adj = (a.p_adj < 0.05),
             category = CASE
               WHEN a.p_adj < 0.05 THEN 'significant post FDR correction'
               WHEN t.p_raw < 0.05 THEN 'significant prior correction'
               ELSE 'not significant' END
        FROM (
          WITH ranked AS (
            SELECT row_id, group_col, p_raw,
                   COUNT(*) OVER (PARTITION BY group_col) AS m,
                   ROW_NUMBER() OVER (PARTITION BY group_col ORDER BY p_raw ASC) AS i
              FROM ", tbl_q, "
             WHERE p_raw IS NOT NULL
          ),
          adj AS (
            SELECT row_id, group_col, p_raw, m, i, (p_raw * m / i) AS bh_raw
              FROM ranked
          ),
          mono AS (
            SELECT row_id, group_col,
                   MIN(bh_raw) OVER (
                     PARTITION BY group_col
                     ORDER BY p_raw DESC
                     ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
                   ) AS p_adj_raw
              FROM adj
          )
          SELECT row_id, LEAST(1.0, p_adj_raw) AS p_adj
            FROM mono
        ) a
       WHERE t.row_id = a.row_id
    "))

      # ---- return --------------------------------------------------------------
      if (isTRUE(collect)) {
        # 1) pull the materialized TABLE lazily, sort for readability
        out_lazy <- dplyr::tbl(con, register_name) |>
          dplyr::arrange(group_col, rank, feature, group1, group2)

        # 2) collect FIRST (avoids translating tidyselect helpers to SQL)
        out <- out_lazy |> dplyr::collect()

        # 3) select + rename in-memory (tidyselect works here)
        out <- out |>
          dplyr::select(
            tidyselect::any_of("view"),
            rank, feature, group_col, group1, group2,
            n1,
            N1 = n1_tot, prop1, percent1,
            n2, N2 = n2_tot, prop2, percent2,
            tidyselect::any_of(c("ratio", "delta_ratio")),
            p_raw, p_adj, passed_raw, passed_adj, category
          )

        .ph_log_ok("Done", bullets = c(
          paste0("rows: ", nrow(out)),
          paste0("ranks: ", paste(unique(out$rank), collapse = ", ")),
          paste0("group_cols: ", paste(unique(out$group_col), collapse = ", ")),
          if (interaction) "mode: interaction-only" else "mode: per-column"
        ))

        return(out)
      } else {
        .ph_log_ok("Materialization complete (TABLE; no lazy recompute)", bullets = register_name)
        return(dplyr::tbl(con, register_name))
      }
    }
  )
}
