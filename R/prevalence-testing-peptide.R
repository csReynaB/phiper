# ==============================================================================
# Global shift in peptide-level prevalence via subject-level permutation
# - NEW: rank_feature_keep — optional named list of values to keep per rank.
# - NEW: weight_mode "n_eff_sqrt" (sqrt of expected positives n1*p1+n2*p2)
# - NEW: prev_strat = "decile" => prevalence-stratified Stouffer (deciles of pooled prev)
# - NEW: stat_mode  = "diff" or "asin" for per-peptide z
# ------------------------------------------------------------------------------

ph_prevalence_shift_permutation <- function(
    x,
    rank_cols,
    group_cols,                     # universes; pairwise contrasts built inside
    exist_col           = "exist",  # 0/1 peptide presence per sample
    interaction         = FALSE,
    combine_cols        = NULL,
    interaction_sep     = "::",
    # ---- permutation/test options ----
    B_permutations      = 2000L,    # number of permutations
    seed                = 1L,       # RNG seed
    smooth_eps_num      = 0.5,      # Laplace smoothing: (x+eps)/(n+2*eps)
    smooth_eps_den_mult = 2.0,
    min_max_prev        = 0.0,      # keep peptides with max(p1,p2) >= this (in [0,1])
    weight_mode         = c("equal","se_invvar","n_eff_sqrt"),  # Stouffer weights
    stat_mode           = c("diff","asin"),                     # per-peptide z
    prev_strat          = c("none","decile"),                   # prevalence stratification
    winsor_z            = 4.0,      # winsorize z at +/- this threshold
    parallel            = NULL,     # future.apply compatible (optional; not used inside)
    collect             = TRUE,     # keep TRUE (materialization not needed here)
    register_name       = NULL,     # reserved (no-op here)
    progress_perm_every = 1L,       # log every k permutations (1 = every perm)
    # ---- keep only these (rank,feature) levels ------------------------------
    rank_feature_keep   = NULL      # named list: rank -> vector of features to keep
) {
  weight_mode <- match.arg(weight_mode)
  stat_mode   <- match.arg(stat_mode)
  prev_strat  <- match.arg(prev_strat)

  .q   <- function(con, nm) as.character(DBI::dbQuoteIdentifier(con, nm))
  .sym <- rlang::sym

  .ph_with_timing(
    headline = "prevalence_shift_permutation (global peptide-level shift)",
    step     = NULL, bullets = NULL,
    expr     = {

      # ---- validate / normalize input ---------------------------------------
      tryCatch({
        chk::chk_character(rank_cols);  chk::chk_true(length(rank_cols) >= 1)
        chk::chk_character(group_cols); chk::chk_true(length(group_cols) >= 1)
        chk::chk_string(exist_col)
        chk::chk_number(B_permutations); chk::chk_true(B_permutations >= 100)
        if (!is.null(rank_feature_keep)) {
          if (!is.list(rank_feature_keep) || is.null(names(rank_feature_keep)))
            .ph_abort("rank_feature_keep must be a *named* list: rank -> values to keep.")
        }
      }, error = function(e) .ph_abort("Invalid arguments", bullets = e$message))

      # long data
      df_long <- try({
        if (inherits(x, "phip_data")) {
          x$data_long |>
            dplyr::select(tidyselect::any_of(c("sample_id","subject_id","peptide_id", exist_col, group_cols)))
        } else {
          chk::chk_data(x)
          need <- c("sample_id","peptide_id", exist_col, group_cols, "subject_id")
          miss <- setdiff(need, colnames(x))
          if (length(miss) && !"subject_id" %in% miss) {
            tibble::as_tibble(x) |>
              dplyr::select(tidyselect::any_of(c("sample_id","peptide_id", exist_col, group_cols, "subject_id")))
          } else if (!length(miss)) {
            tibble::as_tibble(x) |>
              dplyr::select(tidyselect::any_of(need))
          } else {
            .ph_abort("Missing required columns", bullets = paste("-", miss))
          }
        }
      }, silent = TRUE)
      if (inherits(df_long, "try-error")) .ph_abort("Could not prepare input data.")

      df_long <- df_long |>
        dplyr::mutate(!! .sym(exist_col) := dplyr::coalesce(as.integer(!! .sym(exist_col)), 0L))

      # connection (needed only to access peptide_library columns efficiently)
      con <- NULL
      if (inherits(x, "phip_data")) con <- tryCatch(x$meta$con, error = function(...) NULL)
      if (is.null(con) && inherits(df_long, "tbl_sql")) con <- dbplyr::remote_con(df_long)

      # ---- attach ranks from peptide_library as needed (same-src, robust) ---------
      ranks_needing_lib <- setdiff(rank_cols, "peptide_id")
      if (length(ranks_needing_lib)) {
        if (!inherits(x, "phip_data") || is.null(x$peptide_library))
          .ph_abort("Peptide library required for non-peptide ranks.")

        miss_tax <- setdiff(ranks_needing_lib, colnames(x$peptide_library))
        if (length(miss_tax))
          .ph_abort("Requested taxonomy/annotation columns not in peptide_library.",
                    bullets = paste("-", miss_tax))

        # minimal slice of the library
        lib_min <- x$peptide_library |>
          dplyr::select(tidyselect::any_of(c("peptide_id", ranks_needing_lib))) |>
          dplyr::distinct()

        # figure out df_long's connection (must exist for DB-join)
        con_long <- NULL
        if (inherits(x, "phip_data")) con_long <- tryCatch(x$meta$con, error = function(...) NULL)
        if (is.null(con_long) && inherits(df_long, "tbl_sql")) con_long <- dbplyr::remote_con(df_long)
        if (is.null(con_long))
          .ph_abort("No DuckDB connection found for joining the peptide library.")

        # helper: check if two tbl_sql share the same src
        same_src <- function(a, b) {
          inherits(a, "tbl_sql") && inherits(b, "tbl_sql") &&
            identical(dbplyr::remote_con(a), dbplyr::remote_con(b))
        }

        # ensure the library is on the SAME src as df_long
        if (inherits(lib_min, "tbl_sql") && same_src(df_long, lib_min)) {
          lib_src <- lib_min
          drop_stmt <- NULL
        } else {
          # if it's local or on a different src, copy into a unique TEMP table on df_long's con
          ts_suffix <- format(Sys.time(), "%Y%m%d_%H%M%S")
          rnd_suffix <- sprintf("%06d", sample.int(1e6, 1))
          tmp_name <- paste0("ph_tmp_lib_", ts_suffix, "_", rnd_suffix)
          tmp_q    <- as.character(DBI::dbQuoteIdentifier(con_long, tmp_name))

          # defensive drop (just in case)
          try(DBI::dbExecute(con_long, paste0("DROP TABLE IF EXISTS ", tmp_q)), silent = TRUE)

          # write as TEMP (lives only for this connection), overwrite if needed
          if (inherits(lib_min, "tbl_sql")) {
            # collect then write — avoids cross-src auto_copy
            lib_df <- lib_min %>% dplyr::collect()
            DBI::dbWriteTable(con_long, name = tmp_name, value = tibble::as_tibble(lib_df),
                              temporary = TRUE, overwrite = TRUE)
          } else {
            DBI::dbWriteTable(con_long, name = tmp_name, value = tibble::as_tibble(lib_min),
                              temporary = TRUE, overwrite = TRUE)
          }

          # drop on exit
          drop_stmt <- paste0("DROP TABLE IF EXISTS ", tmp_q)
          on.exit(try(DBI::dbExecute(con_long, drop_stmt), silent = TRUE), add = TRUE)

          lib_src <- dplyr::tbl(con_long, tmp_name)
        }

        # (optional) compute df_long to a temp table to stabilize downstream joins
        if (inherits(df_long, "tbl_sql")) df_long <- df_long %>% dplyr::compute()

        # now both are on the same src — safe to join without copy=TRUE
        df_long <- df_long %>% dplyr::left_join(lib_src, by = "peptide_id")
      }
      available_ranks <- intersect(rank_cols, colnames(df_long))
      if (!length(available_ranks)) .ph_abort("None of the requested rank_cols are available.")

      # ---- pivot long on ranks (rank, feature) -------------------------------
      df_ranked_long <- df_long %>%
        tidyr::pivot_longer(
          cols      = tidyselect::all_of(available_ranks),
          names_to  = "rank",
          values_to = "feature"
        ) %>%
        dplyr::mutate(feature = as.character(feature)) %>%
        dplyr::filter(!is.na(feature))

      # ---- optional filter of (rank,feature) (DB-safe: do it after pivot) ----
      if (!is.null(rank_feature_keep)) {
        clauses <- purrr::imap(rank_feature_keep, function(vals, rk) {
          vals_chr <- as.character(vals)
          rlang::expr( (rank == !!rk) & (feature %in% !!vals_chr) )
        })
        filtexpr <- if (length(clauses) == 1L) clauses[[1L]] else purrr::reduce(clauses, function(a,b) rlang::expr( (!!a) | (!!b) ))
        df_ranked_long <- df_ranked_long %>% dplyr::filter(!!filtexpr)
      }

      # ---- build universes (per-column vs interaction) -----------------------
      make_interaction <- function(tbl, c1, c2, sep) {
        comb_name <- paste(c1, c2, sep = " + ")
        tbl |>
          dplyr::filter(!is.na(.data[[c1]]) & !is.na(.data[[c2]])) |>
          dplyr::mutate(
            group_col   = comb_name,
            group_value = paste0(.data[[c1]], sep, .data[[c2]])
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
        .ph_log_info("Grouping universes", bullets = paste0("- ONLY interaction of: ", paste(combine_cols, collapse = " + ")))
      } else if (isTRUE(interaction)) {
        if (length(group_cols) < 2L) .ph_abort("interaction=TRUE needs at least two group_cols.")
        inter_view <- make_interaction(df_ranked_long, group_cols[1], group_cols[2], interaction_sep)
        gs_view <- dplyr::bind_rows(per_column, inter_view)
        .ph_log_info("Grouping universes", bullets = c(
          paste0("- per-column: ", paste(group_cols, collapse = ", ")),
          paste0("- PLUS interaction: ", paste(group_cols[1:2], collapse = " + "))
        ))
      } else {
        gs_view <- per_column
        .ph_log_info("Grouping universes", bullets = paste0("- per-column only: ", paste(group_cols, collapse = ", ")))
      }

      # ---- helper: combined statistic builder --------------------------------
      .combine_T <- function(p1, p2, n1, n2, winsor_z, weight_mode, stat_mode, prev_strat) {
        # per-peptide effect & z
        delta <- p2 - p1
        if (identical(stat_mode, "asin")) {
          z <- (asin(sqrt(p2)) - asin(sqrt(p1))) /
            sqrt(1/(4*pmax(n1,1)) + 1/(4*pmax(n2,1)))
        } else {
          se <- sqrt( (p1*(1-p1)/pmax(n1,1)) + (p2*(1-p2)/pmax(n2,1)) )
          z  <- delta / pmax(se, 1e-12)
        }
        z <- pmax(pmin(z, winsor_z), -winsor_z)

        # weights
        if (identical(weight_mode, "se_invvar") && identical(stat_mode, "asin")) {
          # Asin z already variance-stabilized; still allow mild inv-var using same denominator
          w <- 1 / sqrt(1/(4*pmax(n1,1)) + 1/(4*pmax(n2,1)))
        } else if (identical(weight_mode, "se_invvar")) {
          se <- sqrt( (p1*(1-p1)/pmax(n1,1)) + (p2*(1-p2)/pmax(n2,1)) )
          w  <- 1/pmax(se, 1e-6)
        } else if (identical(weight_mode, "n_eff_sqrt")) {
          w  <- sqrt(pmax(n1*p1 + n2*p2, 1e-12))
        } else {
          w  <- rep(1.0, length(z))
        }

        # combine
        if (identical(prev_strat, "decile")) {
          p_pool <- (n1*p1 + n2*p2) / pmax(n1+n2, 1)
          br <- stats::quantile(p_pool, probs = seq(0,1,length.out = 11), na.rm = TRUE)
          bin <- cut(p_pool, breaks = unique(br), include.lowest = TRUE)
          z_bin <- tapply(seq_along(z), bin, function(ii) {
            if (!length(ii)) return(0)
            wi <- w[ii]; zi <- z[ii]
            sum(wi*zi)/sqrt(sum(wi^2))
          })
          T_obs <- mean(unlist(z_bin), na.rm = TRUE)
        } else {
          T_obs <- sum(w * z) / sqrt(sum(w^2))
        }

        # presentation weights & summaries (independent of prev_strat)
        w_norm <- w / pmax(sum(w), 1e-12)
        list(T_obs = as.numeric(T_obs),
             delta = delta,
             w_norm = w_norm)
      }

      # ---- helper: compute combined Z and permutation p for one contrast -----
      .one_contrast <- function(dat_uv, weight_mode, stat_mode, prev_strat, B, seed,
                                smooth_eps_num, smooth_eps_den_mult,
                                min_max_prev, winsor_z,
                                allow_parallel = FALSE,
                                progress_perm_every = 1L) {

        if (inherits(dat_uv, "tbl_sql") || inherits(dat_uv, "tbl_lazy")) {
          dat_uv <- dat_uv %>% dplyr::collect()
        }

        lv <- sort(unique(dat_uv$group_value))
        if (length(lv) != 2L) return(NULL)
        g1 <- lv[1]; g2 <- lv[2]

        has_subj <- "subject_id" %in% names(dat_uv)
        paired_ok <- FALSE
        if (has_subj) {
          subj_g1 <- dat_uv %>% dplyr::filter(group_value == g1) %>% dplyr::distinct(subject_id) %>% dplyr::pull()
          subj_g2 <- dat_uv %>% dplyr::filter(group_value == g2) %>% dplyr::distinct(subject_id) %>% dplyr::pull()
          paired_ids <- intersect(subj_g1, subj_g2)
          paired_ok <- length(paired_ids) > 0L
        }

        if (paired_ok) {
          dat_a <- dat_uv %>% dplyr::filter(subject_id %in% paired_ids)
          design <- "paired"
          n_subj <- length(paired_ids)
        } else {
          dat_a <- dat_uv
          design <- "unpaired"
          n_subj <- NA_integer_
        }

        prev_tbl <- dat_a %>%
          dplyr::group_by(peptide_id, group_value) %>%
          dplyr::summarise(
            x = sum(exist > 0L),
            n = dplyr::n_distinct(if ("subject_id" %in% names(dat_a)) subject_id else sample_id),
            .groups = "drop"
          ) %>%
          tidyr::pivot_wider(names_from = group_value, values_from = c(x,n), values_fill = 0L)

        if (!all(c(paste0("x_", g1), paste0("x_", g2),
                   paste0("n_", g1), paste0("n_", g2)) %in% names(prev_tbl))) {
          return(NULL)
        }

        eps <- smooth_eps_num
        den_mult <- smooth_eps_den_mult
        p1 <- (prev_tbl[[paste0("x_", g1)]] + eps) / (prev_tbl[[paste0("n_", g1)]] + den_mult*eps)
        p2 <- (prev_tbl[[paste0("x_", g2)]] + eps) / (prev_tbl[[paste0("n_", g2)]] + den_mult*eps)

        keep <- (pmax(p1, p2) >= min_max_prev)
        if (!any(keep)) return(NULL)

        p1 <- p1[keep]; p2 <- p2[keep]
        n1 <- prev_tbl[[paste0("n_", g1)]][keep]; n2 <- prev_tbl[[paste0("n_", g2)]][keep]

        comb <- .combine_T(p1, p2, n1, n2, winsor_z, weight_mode, stat_mode, prev_strat)
        T_obs <- comb$T_obs
        delta <- comb$delta
        w_norm <- comb$w_norm

        mean_delta        <- mean(delta)
        frac_delta_pos    <- mean(delta > 0)
        mean_delta_w      <- sum(w_norm * delta)
        frac_delta_pos_w  <- sum(w_norm * as.numeric(delta > 0))

        set.seed(seed)
        B <- as.integer(B)
        T_perm <- numeric(B)

        if (design == "paired") {
          wide <- dat_a %>%
            dplyr::select(subject_id, peptide_id, group_value, exist) %>%
            tidyr::pivot_wider(
              names_from = c(group_value), values_from = exist, values_fill = 0L,
              id_cols = c(subject_id, peptide_id)
            )
          mat  <- wide %>% dplyr::arrange(subject_id, peptide_id)
          subj <- unique(mat$subject_id)
          Jmap <- split(seq_len(nrow(mat)), mat$peptide_id)

          .T_from_assignment <- function(flip_vec) {
            x1 <- numeric(length(Jmap)); n1p <- numeric(length(Jmap))
            x2 <- numeric(length(Jmap)); n2p <- numeric(length(Jmap))
            names(x1) <- names(Jmap)

            fset <- as.logical(flip_vec); names(fset) <- subj
            for (jj in seq_along(Jmap)) {
              rows <- Jmap[[jj]]
              sids <- mat$subject_id[rows]
              g1v  <- mat[[g1]][rows]; g2v <- mat[[g2]][rows]
              swap <- fset[as.character(sids)]
              a1 <- ifelse(swap, g2v, g1v)
              a2 <- ifelse(swap, g1v, g2v)
              x1[jj]  <- sum(a1); x2[jj] <- sum(a2)
              n1p[jj] <- length(a1); n2p[jj] <- length(a2)
            }
            p1p <- (x1 + eps) / (n1p + den_mult*eps)
            p2p <- (x2 + eps) / (n2p + den_mult*eps)
            keep_p <- (pmax(p1p, p2p) >= min_max_prev)
            if (!any(keep_p)) return(0)
            .combine_T(p1p[keep_p], p2p[keep_p], n1p[keep_p], n2p[keep_p],
                       winsor_z, weight_mode, stat_mode, prev_strat)$T_obs
          }

          for (b in seq_len(B)) {
            flip <- sample(c(FALSE, TRUE), size = length(subj), replace = TRUE)
            T_perm[b] <- .T_from_assignment(flip)

            if (progress_perm_every > 0L && (b %% progress_perm_every) == 0L) {
              .ph_log_info(
                sprintf("[perm] %d/%d (paired; %s vs %s)", b, B, g1, g2),
                bullets = c(
                  paste0("kept peptides: ", length(delta)),
                  paste0("T_obs: ", sprintf("%.3f", T_obs))
                )
              )
            }
          }

        } else {
          subj_g1 <- dat_a %>% dplyr::filter(group_value == g1) %>% dplyr::distinct(sample_id) %>% dplyr::pull()
          subj_g2 <- dat_a %>% dplyr::filter(group_value == g2) %>% dplyr::distinct(sample_id) %>% dplyr::pull()
          n1s <- length(subj_g1); n2s <- length(subj_g2)
          all_ids <- c(subj_g1, subj_g2)

          obs <- dat_a %>%
            dplyr::select(sample_id, peptide_id, group_value, exist)

          .T_from_split <- function(idx1) {
            lab <- rep(NA_character_, length(all_ids)); names(lab) <- all_ids
            lab[idx1] <- g1
            lab[setdiff(seq_along(all_ids), idx1)] <- g2
            lab_df <- tibble::tibble(sample_id = all_ids, group_value_perm = lab)

            tmp <- obs %>%
              dplyr::left_join(lab_df, by = "sample_id") %>%
              dplyr::group_by(peptide_id, group_value_perm) %>%
              dplyr::summarise(x = sum(exist > 0L), n = dplyr::n(), .groups = "drop") %>%
              tidyr::pivot_wider(names_from = group_value_perm, values_from = c(x,n), values_fill = 0L)

            if (!all(c(paste0("x_", g1), paste0("x_", g2),
                       paste0("n_", g1), paste0("n_", g2)) %in% names(tmp))) return(0)

            p1p <- (tmp[[paste0("x_", g1)]] + eps) / (tmp[[paste0("n_", g1)]] + den_mult*eps)
            p2p <- (tmp[[paste0("x_", g2)]] + eps) / (tmp[[paste0("n_", g2)]] + den_mult*eps)
            keep_p <- (pmax(p1p, p2p) >= min_max_prev)
            if (!any(keep_p)) return(0)
            n1p  <- tmp[[paste0("n_", g1)]][keep_p]; n2p <- tmp[[paste0("n_", g2)]][keep_p]
            .combine_T(p1p[keep_p], p2p[keep_p], n1p, n2p,
                       winsor_z, weight_mode, stat_mode, prev_strat)$T_obs
          }

          for (b in seq_len(B)) {
            idx1 <- sample(seq_along(all_ids), size = n1s, replace = FALSE)
            T_perm[b] <- .T_from_split(idx1)

            if (progress_perm_every > 0L && (b %% progress_perm_every) == 0L) {
              .ph_log_info(
                sprintf("[perm] %d/%d (unpaired; %s vs %s)", b, B, g1, g2),
                bullets = c(
                  paste0("kept peptides: ", length(delta)),
                  paste0("T_obs: ", sprintf("%.3f", T_obs))
                )
              )
            }
          }
        }

        p_two <- (1 + sum(abs(T_perm) >= abs(T_obs))) / (1 + B)
        list(
          group1 = g1, group2 = g2,
          design = design,
          n_subjects_paired = if (design=="paired") n_subj else NA_integer_,
          n_peptides_used   = length(delta),
          T_obs = as.numeric(T_obs),
          p_perm = as.numeric(p_two),
          mean_delta = mean_delta,
          frac_delta_pos = frac_delta_pos,
          mean_delta_w = as.numeric(mean_delta_w),
          frac_delta_pos_w = as.numeric(frac_delta_pos_w)
        )
      } # end .one_contrast

      # ---- assemble universes and run contrasts per (rank, feature) ----------
      dat_all <- gs_view %>%
        dplyr::select(sample_id, dplyr::any_of("subject_id"),
                      peptide_id, !! .sym(exist_col), rank, feature, group_col, group_value) %>%
        dplyr::rename(exist = !! .sym(exist_col))

      # keep peptides with at least one event somewhere (speed)
      pep_keep <- dat_all %>%
        dplyr::group_by(peptide_id) %>%
        dplyr::summarise(any_exist = any(exist > 0L), .groups = "drop") %>%
        dplyr::filter(any_exist) %>%
        dplyr::pull(peptide_id)
      dat_all <- dat_all %>% dplyr::filter(peptide_id %in% pep_keep)

      # enumerate strata (materialize in R, then optionally filter by rank_feature_keep)
      strata <- dat_all %>%
        dplyr::distinct(rank, feature, group_col) %>%
        dplyr::collect() %>%
        dplyr::arrange(rank, feature, group_col)

      if (!is.null(rank_feature_keep) && length(rank_feature_keep)) {
        allow <- purrr::imap_dfr(rank_feature_keep, function(val, rk) {
          tibble::tibble(rank = rk, feature = as.character(val))
        })
        strata <- dplyr::semi_join(strata, allow, by = c("rank","feature"))
      }

      if (nrow(strata) == 0L) {
        .ph_abort("No valid strata to test (check group levels / ranks).")
      }

      results <- vector("list", nrow(strata))

      .ph_log_info("Running permutation shift test per (rank, feature, group_col)",
                   bullets = c(paste0("B: ", B_permutations),
                               paste0("weight_mode: ", weight_mode),
                               paste0("stat_mode: ", stat_mode),
                               paste0("prev_strat: ", prev_strat),
                               paste0("winsor_z: ", winsor_z),
                               paste0("min_max_prev: ", min_max_prev)))

      for (i in seq_len(nrow(strata))) {
        st <- strata[i, ]

        dat_uv <- dat_all %>%
          dplyr::filter(
            rank      == st$rank,
            feature   == st$feature,
            group_col == st$group_col
          ) %>%
          dplyr::collect()

        lv <- dat_uv %>%
          dplyr::distinct(group_value) %>%
          dplyr::arrange(group_value) %>%
          dplyr::pull(group_value)

        if (length(lv) < 2L) { results[[i]] <- NULL; next }
        if (length(lv) > 2L) { results[[i]] <- NULL; next }

        out <- .one_contrast(
          dat_uv              = dat_uv,
          weight_mode         = weight_mode,
          stat_mode           = stat_mode,
          prev_strat          = prev_strat,
          B                   = B_permutations,
          seed                = seed + i - 1L,
          smooth_eps_num      = smooth_eps_num,
          smooth_eps_den_mult = smooth_eps_den_mult,
          min_max_prev        = min_max_prev,
          winsor_z            = winsor_z,
          allow_parallel      = isTRUE(parallel),
          progress_perm_every = progress_perm_every
        )

        if (is.null(out)) { results[[i]] <- NULL; next }
        results[[i]] <- dplyr::bind_cols(st, tibble::as_tibble(out))
      }

      res <- results |> purrr::compact() |> dplyr::bind_rows()

      if (!nrow(res)) {
        .ph_abort("No valid contrasts to test (check group levels / ranks).")
      }

      # multiplicity (per rank)
      res <- res %>%
        dplyr::group_by(rank) %>%
        dplyr::mutate(p_adj_rank = p.adjust(p_perm, method = "BH")) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          category_rank_bh = dplyr::case_when(
            p_adj_rank < 0.05 ~ "significant (BH, per rank)",
            p_perm    < 0.05 ~ "nominal only",
            TRUE ~ "not significant"
          )
        )

      .ph_log_ok("Done (permutation shift)",
                 bullets = c(
                   paste0("rows: ", nrow(res)),
                   paste0("ranks: ", paste(unique(res$rank), collapse = ", ")),
                   paste0("B: ", B_permutations)
                 ))

      res %>%
        dplyr::select(
          rank, feature, group_col, group1, group2, design,
          n_subjects_paired, n_peptides_used,
          T_obs, p_perm, p_adj_rank,
          mean_delta, frac_delta_pos,
          mean_delta_w, frac_delta_pos_w,
          category_rank_bh
        ) %>%
        dplyr::arrange(rank, feature, group_col, group1, group2)
    }
  )
}
