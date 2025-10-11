# ==============================================================================
# Prevalence (counts/percent) + pairwise Fisher tests across multiple group vars
# POP framework (Peptide-OR-Presence), DuckDB-first, chunked p-values in R
# Single FDR family PER RANK (across all requested universes); BH and weighted BH
# ==============================================================================

#' @title Prevalence by group with pairwise tests (multi-rank; POP; FDR per-rank)
#' @description
#' - POP: a (rank,feature) is present in a sample if ≥1 peptide of that feature is present (exist > 0).
#' - Pairwise group comparisons: Fisher's exact test only.
#' - FDR is computed PER RANK across ALL hypotheses of that rank:
#'   * all features (POOL) × all pairwise group contrasts across universes
#'   * universes = each group_col separately, plus interaction if interaction=TRUE,
#'     or only the combined interaction if combine_cols is provided.
#' - Weighted BH: weights ∝ #peptides per (rank,feature), assigned to all tests from that feature,
#'   then scaled so that sum(weights) = m_rank (number of tests in that rank).
#' - Logs:
#'   * POOL per rank (#unique features),
#'   * levels and pair counts per universe,
#'   * total m per rank.
#'
#' @param x             phip_data or DuckDB-backed long df with sample_id, peptide_id, exist_col, group_cols
#' @param rank_cols     character vector (e.g., c("species","genus","family"))
#' @param group_cols    character vector of grouping columns
#' @param exist_col     presence indicator (default "exist")
#' @param weight_mode   c("peptide_count","none") for weighted BH
#' @param parallel      NULL=auto via future, or TRUE/FALSE
#' @param compute_ratios_db logical; if TRUE add ratio/delta_ratio (epsilon-stabilized) in DB
#' @param interaction   if TRUE, include BOTH: per-column universes + interaction of two cols
#' @param combine_cols  optional length-2 subset of group_cols; if provided, include ONLY the combined interaction
#' @param interaction_sep separator for interaction levels
#' @param collect       if TRUE return a tibble; else a DuckDB TABLE
#' @param register_name optional TABLE name if collect=FALSE
#' @return tibble or lazy tbl with raw p, BH per-rank, weighted BH per-rank
#' @export
ph_prevalence_compare <- function(x,
                                  rank_cols,
                                  group_cols,
                                  exist_col         = "exist",
                                  weight_mode       = c("peptide_count","none"),
                                  parallel          = NULL,
                                  compute_ratios_db = TRUE,
                                  interaction       = FALSE,
                                  combine_cols      = NULL,
                                  interaction_sep   = "::",
                                  collect           = TRUE,
                                  register_name     = NULL) {

  weight_mode <- match.arg(weight_mode)

  .ph_with_timing(
    headline = "prevalence_compare (per-rank FDR)",
    step     = NULL,
    bullets  = NULL,
    expr     = {

      .q   <- function(con, nm) as.character(DBI::dbQuoteIdentifier(con, nm))
      .sym <- rlang::sym
      chunk_n <- getOption("phiper.prev.chunk", 1e6L)

      # ---- args / logging ----------------------------------------------------
      tryCatch({
        chk::chk_character(rank_cols);  chk::chk_true(length(rank_cols) >= 1)
        chk::chk_character(group_cols); chk::chk_true(length(group_cols) >= 1)
        chk::chk_string(exist_col)

        if (!is.null(combine_cols)) {
          if (length(combine_cols) != 2L)
            .ph_abort("combine_cols must be length 2 when provided.")
          miss_int <- setdiff(combine_cols, group_cols)
          if (length(miss_int))
            .ph_abort("combine_cols must be subset of group_cols.", bullets = paste("-", miss_int))
        }
      }, error = function(e) .ph_abort("Invalid arguments", bullets = e$message))

      .ph_log_info(
        "Preparing input data",
        bullets = c(
          paste0("ranks: ", paste(rank_cols, collapse = ", ")),
          paste0("group_cols: ", paste(group_cols, collapse = ", ")),
          if (!is.null(combine_cols)) paste0("combine_cols: ", paste(combine_cols, collapse = " + "))
          else if (interaction) "interaction: TRUE (additive to per-column)",
          paste0("exist_col: ", exist_col),
          paste0("weight_mode: ", weight_mode),
          paste0("collect: ", collect)
        )[!vapply(c(
          paste0("ranks: ", paste(rank_cols, collapse = ", ")),
          paste0("group_cols: ", paste(group_cols, collapse = ", ")),
          if (!is.null(combine_cols)) paste0("combine_cols: ", paste(combine_cols, collapse = " + "))
          else if (interaction) "interaction: TRUE (additive to per-column)" else NA_character_,
          paste0("exist_col: ", exist_col),
          paste0("weight_mode: ", weight_mode),
          paste0("collect: ", collect)
        ), is.na, logical(1))]
      )

      # ---- normalize to long df ----------------------------------------------
      df_long <- try({
        if (inherits(x, "phip_data")) {
          x$data_long |>
            dplyr::select(tidyselect::any_of(c("sample_id","peptide_id", exist_col, group_cols)))
        } else {
          chk::chk_data(x)
          need <- c("sample_id","peptide_id", exist_col, group_cols)
          miss <- setdiff(need, colnames(x))
          if (length(miss)) .ph_abort("Missing required columns", bullets = paste("-", miss))
          tibble::as_tibble(x) |>
            dplyr::select(tidyselect::any_of(need))
        }
      }, silent = TRUE)
      if (inherits(df_long, "try-error")) .ph_abort("Could not prepare input data.")

      df_long <- df_long |>
        dplyr::mutate(!! .sym(exist_col) := dplyr::coalesce(as.integer(!! .sym(exist_col)), 0L))

      # connection
      con <- NULL
      if (inherits(x, "phip_data")) con <- tryCatch(x$meta$con, error = function(...) NULL)
      if (is.null(con) && inherits(df_long, "tbl_sql")) con <- dbplyr::remote_con(df_long)
      if (is.null(con)) .ph_abort("No DuckDB connection found. Use <phip_data>$meta$con or a DuckDB-backed tbl_sql.")

      view_const <- if (inherits(x, "phip_data")) attr(x, "view") %||% (x$meta$view %||% NA_character_) else NA_character_

      # ---- ranks & library ---------------------------------------------------
      ranks_needing_lib <- setdiff(rank_cols, "peptide_id")
      if (length(ranks_needing_lib)) {
        if (!inherits(x, "phip_data") || is.null(x$peptide_library))
          .ph_abort("Peptide library required for non-peptide ranks.")
        miss_tax <- setdiff(ranks_needing_lib, colnames(x$peptide_library))
        if (length(miss_tax)) .ph_abort("Requested taxonomy columns not in peptide_library.", bullets = paste("-", miss_tax))
      }

      df_ranked <- df_long
      if (length(ranks_needing_lib)) {
        lib_cols <- c("peptide_id", ranks_needing_lib)
        lib_min  <- x$peptide_library |>
          dplyr::select(tidyselect::all_of(lib_cols)) |>
          dplyr::distinct()
        df_ranked <- df_ranked |>
          dplyr::left_join(lib_min, by = "peptide_id", copy = TRUE)
      }

      available_ranks <- intersect(rank_cols, colnames(df_ranked))
      if (!length(available_ranks)) .ph_abort("None of the requested rank_cols are available.")
      .ph_log_info("Ranks resolved", bullets = paste("- available:", paste(available_ranks, collapse = ", ")))

      df_ranked_long <- df_ranked %>%
        tidyr::pivot_longer(
          cols      = tidyselect::all_of(available_ranks),
          names_to  = "rank",
          values_to = "feature"
        ) %>%
        dplyr::mutate(feature = as.character(feature))  # ensure TEXT/VARCHAR in DB

      # ---- construct universes (per-column and/or interaction) ---------------
      make_interaction <- function(tbl, c1, c2, sep) {
        comb_name <- paste(c1, c2, sep = " + ")
        tbl |>
          dplyr::filter(!is.na(!! .sym(c1)) & !is.na(!! .sym(c2))) |>
          dplyr::mutate(
            group_col   = comb_name,
            group_value = paste0(!! .sym(c1), sep, !! .sym(c2))
          )
      }

      per_column <- df_ranked_long |>
        tidyr::pivot_longer(
          cols      = tidyselect::all_of(group_cols),
          names_to  = "group_col",
          values_to = "group_value"
        ) |>
        dplyr::filter(!is.na(group_value))

      if (!is.null(combine_cols)) {
        gs_view <- make_interaction(df_ranked_long, combine_cols[1], combine_cols[2], interaction_sep)
        .ph_log_info("Grouping universes", bullets = c(
          paste0("- ONLY interaction of: ", paste(combine_cols, collapse = " + "))
        ))
      } else if (isTRUE(interaction)) {
        # include BOTH per-column and the interaction of the first two group_cols
        if (length(group_cols) < 2L) .ph_abort("interaction=TRUE needs at least two group_cols.")
        inter_view <- make_interaction(df_ranked_long, group_cols[1], group_cols[2], interaction_sep)
        gs_view <- dplyr::bind_rows(per_column, inter_view)
        .ph_log_info("Grouping universes", bullets = c(
          paste0("- per-column: ", paste(group_cols, collapse = ", ")),
          paste0("- PLUS interaction: ", paste(group_cols[1:2], collapse = " + "))
        ))
      } else {
        gs_view <- per_column
        .ph_log_info("Grouping universes", bullets = c(
          paste0("- per-column only: ", paste(group_cols, collapse = ", "))
        ))
      }

      # ---- denominators & POP counts -----------------------------------------
      .ph_log_info("Computing cohort sizes (N) per universe")
      group_sizes <- gs_view |>
        dplyr::distinct(group_col, group_value, sample_id) |>
        dplyr::count(group_col, group_value, name = "N")

      gs_n <- group_sizes |> dplyr::summarise(n = dplyr::n()) |> dplyr::collect()
      if (gs_n$n == 0L) .ph_abort("No non-missing group values; cannot compute denominators.")

      .ph_log_info("Counting present samples per feature (POP)")
      present_counts <- gs_view |>
        dplyr::filter(!! .sym(exist_col) > 0) |>
        dplyr::distinct(group_col, group_value, rank, feature, sample_id) |>
        dplyr::count(group_col, group_value, rank, feature, name = "n_present")

      features_per_rank <- present_counts |> dplyr::distinct(rank, feature)

      # ---- POOL + universe sizes logging -------------------------------------
      pool_tbl <- features_per_rank |>
        dplyr::count(rank, name = "POOL") |>
        dplyr::arrange(rank) |>
        dplyr::collect()

      lev_tbl <- group_sizes |>
        dplyr::count(group_col, name = "k_levels") |>
        dplyr::mutate(n_pairs = dplyr::if_else(k_levels >= 2, (k_levels * (k_levels - 1)) / 2, 0)) |>
        dplyr::collect()

      sum_pairs <- sum(lev_tbl$n_pairs, na.rm = TRUE)
      log_bullets <- c(
        paste0("POOL per rank: ", paste(paste0(pool_tbl$rank, "=", pool_tbl$POOL), collapse = "; ")),
        paste0("Universes: ", paste(paste0(lev_tbl$group_col, " (k=", lev_tbl$k_levels, ", pairs=", lev_tbl$n_pairs, ")"), collapse = "; ")),
        paste0("Pairs per feature (sum across universes): ", sum_pairs),
        paste0("Total tests m per rank = POOL * pairs: ",
               paste(paste0(pool_tbl$rank, "=", pool_tbl$POOL * sum_pairs), collapse = "; "))
      )
      .ph_log_info("FDR accounting", bullets = log_bullets)

      # ---- assemble grid & pairwise summaries --------------------------------
      base_grid <- group_sizes |>
        dplyr::select(group_col, group_value) |>
        dplyr::distinct() |>
        dplyr::mutate(.dummy = 1L) |>
        dplyr::inner_join(features_per_rank |> dplyr::mutate(.dummy = 1L), by = ".dummy") |>
        dplyr::select(group_col, group_value, rank, feature)

      stats_long <- base_grid |>
        dplyr::left_join(present_counts, by = c("group_col","group_value","rank","feature")) |>
        dplyr::left_join(group_sizes,   by = c("group_col","group_value")) |>
        dplyr::mutate(
          n_present = dplyr::coalesce(n_present, 0L),
          prop      = dplyr::if_else(N > 0, n_present / N, NA_real_),
          percent   = 100 * prop
        )

      if (!is.na(view_const)) {
        stats_long <- stats_long |> dplyr::mutate(view = !! view_const, .before = 1)
      }

      .ph_log_info("Building pairwise comparisons")
      by_cols  <- intersect(c("view","rank","feature","group_col"), colnames(stats_long))
      has_view <- "view" %in% colnames(stats_long)

      pairs_joined <- stats_long |>
        dplyr::inner_join(stats_long, by = by_cols, suffix = c("_1", "_2")) |>
        dplyr::filter(group_value_1 != group_value_2) |>
        dplyr::mutate(
          group1_chr = dplyr::if_else(group_value_1 <= group_value_2, group_value_1, group_value_2),
          group2_chr = dplyr::if_else(group_value_1 <= group_value_2, group_value_2, group_value_1),

          n1_val  = dplyr::if_else(group_value_1 <= group_value_2, n_present_1, n_present_2),
          n1_tot  = dplyr::if_else(group_value_1 <= group_value_2, N_1,        N_2),
          p1_val  = dplyr::if_else(group_value_1 <= group_value_2, prop_1,     prop_2),
          pct1val = dplyr::if_else(group_value_1 <= group_value_2, percent_1,  percent_2),

          n2_val  = dplyr::if_else(group_value_1 <= group_value_2, n_present_2, n_present_1),
          n2_tot  = dplyr::if_else(group_value_1 <= group_value_2, N_2,         N_1),
          p2_val  = dplyr::if_else(group_value_1 <= group_value_2, prop_2,      prop_1),
          pct2val = dplyr::if_else(group_value_1 <= group_value_2, percent_2,   percent_1)
        ) |>
        dplyr::filter(group_value_1 == group1_chr)

      pairs_lazy <- (if (has_view) {
        pairs_joined |>
          dplyr::transmute(
            view, rank, feature, group_col,
            group1 = group1_chr, n1 = n1_val, n1_tot = n1_tot, prop1 = p1_val, percent1 = pct1val,
            group2 = group2_chr, n2 = n2_val, n2_tot = n2_tot, prop2 = p2_val, percent2 = pct2val
          )
      } else {
        pairs_joined |>
          dplyr::transmute(
            rank, feature, group_col,
            group1 = group1_chr, n1 = n1_val, n1_tot = n1_tot, prop1 = p1_val, percent1 = pct1val,
            group2 = group2_chr, n2 = n2_val, n2_tot = n2_tot, prop2 = p2_val, percent2 = pct2val
          )
      }) |>
        dbplyr::window_order(group_col, rank, feature, group1, group2) |>
        dplyr::mutate(row_id = dplyr::row_number()) |>
        dplyr::filter(n1_tot > 0, n2_tot > 0)

      if (isTRUE(compute_ratios_db)) {
        pairs_lazy <- pairs_lazy |>
          dplyr::mutate(
            prop1_eps = dplyr::if_else(n1 == 0L, (n1 + 0.5)/n1_tot, prop1),
            prop2_eps = dplyr::if_else(n2 == 0L, (n2 + 0.5)/n2_tot, prop2),
            ratio     = prop1_eps / prop2_eps,
            d1        = dplyr::if_else(n1 == 0L, (n1 + 1.0)/n1_tot, prop1),
            d2        = dplyr::if_else(n2 == 0L, (n2 + 1.0)/n2_tot, prop2),
            delta_ratio = dplyr::if_else(d1 >= d2, d1/d2 - 1, -(d2/d1 - 1))
          )
      }

      # ---- materialize TABLE -------------------------------------------------
      if (is.null(register_name)) register_name <- paste0("ph_prev_", format(Sys.time(), "%Y%m%d_%H%M%S"))
      tbl_q <- .q(con, register_name)

      sql_body <- dbplyr::sql_render(pairs_lazy, con)
      DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", tbl_q))
      DBI::dbExecute(con, paste0("CREATE TABLE ", tbl_q, " AS ", sql_body))

      .ph_log_ok(
        "Materialized DuckDB TABLE",
        bullets = c(paste0("name: ", register_name),
                    "computing p-values (Fisher-only); then FDR per rank (BH / wBH)")
      )

      # ---- P-values (chunked; Fisher-only) ----------------------------------
      p_core_fisher <- function(n1, N1, n2, N2) {
        n1 <- as.double(n1); N1 <- as.double(N1)
        n2 <- as.double(n2); N2 <- as.double(N2)
        if (any(!is.finite(c(n1,N1,n2,N2)))) return(NA_real_)
        if (N1 <= 0 || N2 <= 0) return(NA_real_)
        a  <- max(0, min(n1, N1)); b <- max(0, min(n2, N2))
        c1 <- N1 - a;              c2 <- N2 - b
        out <- try(stats::fisher.test(matrix(c(a, c1, b, c2), 2, byrow = TRUE))$p.value, silent = TRUE)
        if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
      }

      for (col in c("p_raw","p_adj_rank","p_adj_rank_wbh",
                    "passed_rank_bh","passed_rank_wbh",
                    "category_rank_bh","category_rank_wbh")) {
        DBI::dbExecute(
          con,
          paste0(
            "ALTER TABLE ", tbl_q,
            " ADD COLUMN IF NOT EXISTS ", col, " ",
            if (grepl("^p_", col)) "DOUBLE" else if (grepl("^passed", col)) "BOOLEAN" else "VARCHAR"
          )
        )
      }

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
      } else isTRUE(parallel)
      if (want_parallel &&
          (!requireNamespace("future", quietly = TRUE) ||
           !requireNamespace("future.apply", quietly = TRUE))) {
        .ph_abort("Parallel requested but {future}/{future.apply} not available.")
      }

      repeat {
        chunk <- DBI::dbFetch(rs, n = chunk_n)
        if (!nrow(chunk)) break

        colnames(chunk) <- tolower(colnames(chunk))
        req  <- c("row_id","n1","n1_tot","n2","n2_tot")
        miss <- setdiff(req, colnames(chunk))
        if (length(miss)) .ph_abort("Fetched chunk is missing required columns.", bullets = paste("-", miss))

        chunk$row_id <- as.numeric(chunk$row_id)
        n <- nrow(chunk)

        if (want_parallel && tryCatch(future::nbrOfWorkers() > 1L, error = function(...) FALSE)) {
          splits <- split(seq_len(n), cut(seq_len(n), breaks = future::nbrOfWorkers(), labels = FALSE))
          parts <- future.apply::future_lapply(
            splits,
            FUN = function(ii, dat) {
              vapply(ii, function(i) p_core_fisher(dat$n1[i], dat$n1_tot[i], dat$n2[i], dat$n2_tot[i]), numeric(1))
            },
            dat = chunk, future.seed = TRUE, future.scheduling = 1, future.globals = FALSE
          )
          pvals <- unlist(parts, use.names = FALSE)
        } else {
          pvals <- vapply(
            seq_len(n),
            function(i) p_core_fisher(chunk$n1[i], chunk$n1_tot[i], chunk$n2[i], chunk$n2_tot[i]),
            numeric(1)
          )
        }

        DBI::dbAppendTable(
          con,
          name  = gsub('^\"|\"$', '', tmp_q),
          value = data.frame(row_id = chunk$row_id, p_raw = as.numeric(pvals))
        )
      }

      DBI::dbExecute(con, paste0(
        "UPDATE ", tbl_q, " t SET p_raw = s.p_raw FROM ", tmp_q, " s WHERE t.row_id = s.row_id"
      ))
      DBI::dbExecute(con, paste0("DROP TABLE ", tmp_q))

      # ---- Build weights per (rank,feature) ---------------------------------
      weight_tbl_name <- paste0(register_name, "_weights")
      wq <- .q(con, weight_tbl_name)
      DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", wq))

      if (weight_mode == "peptide_count" && length(ranks_needing_lib)) {

        lib_cols_needed <- c("peptide_id", ranks_needing_lib)
        lib_src <- x$peptide_library %>%
          dplyr::select(tidyselect::all_of(lib_cols_needed)) %>%
          dplyr::distinct()

        tmp_lib <- .q(con, paste0(register_name, "_libtmp"))
        DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", tmp_lib))

        if (inherits(lib_src, "tbl_sql")) {
          lib_df <- lib_src %>% dplyr::collect()
          DBI::dbWriteTable(con, name = gsub('^\"|\"$', '', tmp_lib),
                            value = tibble::as_tibble(lib_df),
                            temporary = TRUE)
        } else {
          DBI::dbWriteTable(con, name = gsub('^\"|\"$', '', tmp_lib),
                            value = tibble::as_tibble(lib_src),
                            temporary = TRUE)
        }

        bad <- ranks_needing_lib[!grepl("^[A-Za-z][A-Za-z0-9_]*$", ranks_needing_lib)]
        if (length(bad)) {
          .ph_abort("Taxa column names must be alphanumeric/underscore (start with a letter).",
                    bullets = paste("-", bad))
        }

        selects <- vapply(
          ranks_needing_lib,
          function(rc) {
            rc_q <- .q(con, rc)
            paste0(
              "SELECT '", rc, "' AS rank, ",
              rc_q, " AS feature, peptide_id FROM ", tmp_lib,
              " WHERE ", rc_q, " IS NOT NULL"
            )
          },
          character(1)
        )
        sql_union <- paste(selects, collapse = " UNION ALL ")

        DBI::dbExecute(con, paste0(
          "CREATE TABLE ", wq, " AS ",
          "SELECT rank, feature, COUNT(DISTINCT peptide_id) AS n_peptides ",
          "FROM (", sql_union, ") ",
          "GROUP BY rank, feature"
        ))

        if ("peptide_id" %in% available_ranks) {
          DBI::dbExecute(con, paste0("
            INSERT INTO ", wq, " (rank, feature, n_peptides)
            SELECT 'peptide_id' AS rank, feature, 1 AS n_peptides
            FROM (
              SELECT DISTINCT peptide_id AS feature FROM ", tbl_q, " WHERE rank = 'peptide_id'
            )
            ON CONFLICT DO NOTHING
          "))
        }

        DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS ", tmp_lib))

      } else {
        DBI::dbExecute(con, paste0("
          CREATE TABLE ", wq, " AS
          SELECT DISTINCT rank, feature, 1::INTEGER AS n_peptides
          FROM ", tbl_q, "
        "))
      }

      # ---- BH per rank -------------------------------------------------------
      DBI::dbExecute(con, paste0("
        UPDATE ", tbl_q, " t
           SET p_adj_rank = a.p_adj
          FROM (
            WITH ranked AS (
              SELECT row_id, ", if (has_view) "view," else "", " rank, p_raw,
                     COUNT(*) OVER (PARTITION BY ", if (has_view) "view," else "", " rank) AS m,
                     ROW_NUMBER() OVER (PARTITION BY ", if (has_view) "view," else "", " rank ORDER BY p_raw ASC) AS i
                FROM ", tbl_q, "
               WHERE p_raw IS NOT NULL
            ),
            adj AS (
              SELECT row_id, (p_raw * m / i) AS bh_raw
                FROM ranked
            ),
            mono AS (
              SELECT r.row_id,
                     LEAST(1.0,
                       MIN(a.bh_raw) OVER (
                         PARTITION BY ", if (has_view) "r.view," else "", " r.rank
                         ORDER BY r.p_raw DESC
                         ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
                       )
                     ) AS p_adj
                FROM ranked r
                JOIN adj a USING (row_id)
            )
            SELECT row_id, p_adj FROM mono
          ) a
         WHERE t.row_id = a.row_id
      "))

      DBI::dbExecute(con, paste0("
        UPDATE ", tbl_q, "
           SET passed_rank_bh    = (p_adj_rank < 0.05),
               category_rank_bh  = CASE
                   WHEN p_adj_rank < 0.05 THEN 'significant (BH, per rank)'
                   WHEN p_raw      < 0.05 THEN 'nominal only'
                   ELSE 'not significant' END
      "))

      # ---- Weighted BH per rank ---------------------------------------------
      DBI::dbExecute(con, paste0("
        UPDATE ", tbl_q, " t
           SET p_adj_rank_wbh = z.p_adj_wbh,
               passed_rank_wbh = (z.p_adj_wbh < 0.05),
               category_rank_wbh = CASE
                 WHEN z.p_adj_wbh < 0.05 THEN 'significant (wBH, per rank)'
                 WHEN t.p_raw < 0.05 THEN 'nominal only'
                 ELSE 'not significant' END
          FROM (
            WITH base AS (
              SELECT t.row_id, ", if (has_view) "t.view," else "", " t.rank, t.feature, t.p_raw,
                     COALESCE(w.n_peptides, 1)::DOUBLE AS w_base
                FROM ", tbl_q, " t
                LEFT JOIN ", wq, " w
                  ON t.rank = w.rank AND t.feature = w.feature
            ),
            w_sums AS (
              SELECT ", if (has_view) "view," else "", " rank,
                     COUNT(*)::DOUBLE AS m,
                     SUM(w_base)::DOUBLE AS sum_w
                FROM base
               GROUP BY ", if (has_view) "view," else "", " rank
            ),
            w_scaled AS (
              SELECT b.row_id, ", if (has_view) "b.view," else "", " b.rank, b.p_raw,
                     (b.w_base * s.m / NULLIF(s.sum_w,0.0)) AS w,
                     (b.p_raw / NULLIF((b.w_base * s.m / NULLIF(s.sum_w,0.0)),0.0)) AS p_over_w,
                     ROW_NUMBER() OVER (
                       PARTITION BY ", if (has_view) "b.view," else "", " b.rank
                       ORDER BY (b.p_raw / NULLIF((b.w_base * s.m / NULLIF(s.sum_w,0.0)),0.0)) ASC
                     ) AS i,
                     s.m
                FROM base b
                JOIN w_sums s
                  ON ", if (has_view) "b.view = s.view AND ", "", " b.rank = s.rank
            ),
            adj AS (
              SELECT row_id, ", if (has_view) "view," else "", " rank, m, p_over_w,
                     (m * p_over_w / i) AS wbh_raw
                FROM w_scaled
            ),
            final AS (
              SELECT row_id,
                     LEAST(1.0,
                       MIN(wbh_raw) OVER (
                         PARTITION BY ", if (has_view) "view," else "", " rank
                         ORDER BY p_over_w DESC
                         ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
                       )
                     ) AS p_adj_wbh
                FROM adj
            )
            SELECT row_id, p_adj_wbh FROM final
          ) z
         WHERE t.row_id = z.row_id
      "))

      # ---- return ------------------------------------------------------------
      if (isTRUE(collect)) {
        # ONLY CHANGE: join weights to expose n_peptides next to feature
        wq_name <- gsub('^\"|\"$', '', wq)

        out <- dplyr::tbl(con, register_name) |>
          dplyr::left_join(dplyr::tbl(con, wq_name), by = c("rank","feature")) |>
          dplyr::arrange(rank, feature, group_col, group1, group2) |>
          dplyr::collect()

        out <- out |>
          dplyr::select(
            tidyselect::any_of("view"),
            rank, feature, n_peptides,   # <- added here, right after feature
            group_col, group1, group2,
            n1, N1 = n1_tot, prop1, percent1,
            n2, N2 = n2_tot, prop2, percent2,
            tidyselect::any_of(c("ratio","delta_ratio")),
            p_raw,
            p_adj_rank, passed_rank_bh, category_rank_bh,
            p_adj_rank_wbh, passed_rank_wbh, category_rank_wbh
          )

        .ph_log_ok(
          "Done (per-rank FDR)",
          bullets = c(
            paste0("rows: ", nrow(out)),
            paste0("ranks: ", paste(unique(out$rank), collapse = ", ")),
            paste0("POOL: ", paste(paste0(pool_tbl$rank, "=", pool_tbl$POOL), collapse = "; ")),
            paste0("pairs per feature across universes: ", sum_pairs)
          )
        )
        return(out)
      } else {
        .ph_log_ok("Materialization complete (TABLE; p/q per-rank ready)", bullets = register_name)
        return(dplyr::tbl(con, register_name))
      }
    }
  )
}
