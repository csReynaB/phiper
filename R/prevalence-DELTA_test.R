# ==============================================================================
#' Global shift in peptide-level prevalence via subject-level permutation
#'
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom doParallel registerDoParallel
#'
#' @description
#' Tests for a **global** (genome/antigen-wide) shift in peptide-level prevalence
#' between two groups within each `(rank, feature, group_col)` stratum by
#' aggregating per-peptide effects into a single Stouffer-type statistic and
#' assessing significance via **label permutation** (paired or unpaired, chosen
#' automatically from the data).
#'
#' @details
#' **Overview.** For each stratum `(rank, feature, group_col)` that yields a
#' unique ordered pair of groups `(g1, g2)`, the method:
#'
#' 1) **Peptide-level prevalence with smoothing.**
#'    For each peptide, compute smoothed prevalences in `g1` and `g2`:
#'    \deqn{ \hat p_k = \frac{x_k + \varepsilon}{n_k + \lambda \varepsilon}, \; k \in \{g1,g2\} , }
#'    where `x_k` is the number of samples with presence (`exist > 0`), `n_k`
#'    is the number of (distinct) samples (or paired subjects), `\varepsilon`
#'    is `smooth_eps_num`, and `\lambda` is `smooth_eps_den_mult`. Peptides with
#'    `max( \hat p_{g1}, \hat p_{g2} ) < min_max_prev` are discarded.
#'
#' 2) **Per-peptide effect and z-score (`stat_mode`).**
#'    - `"diff"`: effect \eqn{\Delta = \hat p_{g2} - \hat p_{g1}} with
#'      \eqn{ SE = \sqrt{\hat p_{g1}(1-\hat p_{g1})/n_{g1} + \hat p_{g2}(1-\hat p_{g2})/n_{g2}} } and
#'      \eqn{ z = \Delta / SE } (guarded by a tiny floor).
#'    - `"asin"`: variance-stabilized difference of angles:
#'      \eqn{ z = \frac{\arcsin\sqrt{\hat p_{g2}} - \arcsin\sqrt{\hat p_{g1}}}{ \sqrt{ 1/(4n_{g1}) + 1/(4n_{g2}) } } .}
#'    Each \eqn{z} is **winsorized** to \eqn{\pm} `winsor_z`.
#'
#' 3) **Weights (`weight_mode`).**
#'    - `"equal"`: all peptides weight 1.
#'    - `"se_invvar"`: inverse standard error (or the analogous denominator for
#'      `"asin"`).
#'    - `"n_eff_sqrt"`: \eqn{ \sqrt{ n_{g1}\hat p_{g1} + n_{g2}\hat p_{g2} } } (proxy for expected positives).
#'
#' 4) **Combine into one test statistic (Stouffer).**
#'    Let \eqn{w_i} be weights and \eqn{z_i} peptide z-scores.
#'    - If `prev_strat = "none"`:
#'      \deqn{ T_{\mathrm{obs}} = \frac{ \sum_i w_i z_i }{ \sqrt{ \sum_i w_i^2 } } . }
#'    - If `prev_strat = "decile"`: bin peptides by deciles of pooled prevalence
#'      \eqn{ (\;n_{g1}\hat p_{g1} + n_{g2}\hat p_{g2}\;)/(\;n_{g1}+n_{g2}\;) } and combine a
#'      Stouffer \eqn{T_b} within each bin; report the **mean** of bin-level z’s
#'      as \eqn{T_{\mathrm{obs}}} (stabilizes across prevalence regimes).
#'
#' 5) **Permutation scheme and p-value.**
#'    We compute a **two-sided** permutation p-value for the global shift:
#'
#'    - **Paired design:** when both groups have measurements for the same
#'      `subject_id`, we independently flip each subject’s labels (`g1`↔`g2`)
#'      with 50/50 probability, recompute all steps (1–4) to get \eqn{T_b}.
#'    - **Unpaired design:** we shuffle sample labels while **preserving group
#'      sizes** (draw a random split of samples into `n1`/`n2`), recompute (1–4)
#'      to get \eqn{T_b}.
#'
#'    With \eqn{B} permutations, the p-value is
#'    \deqn{ p = \frac{1 + \sum_{b=1}^{B} \mathbf{1}\{ |T_b| \ge |T_{\mathrm{obs}}| \} }{1 + B} , }
#'    which is the standard *add-one* estimate ensuring non-zero p under the
#'    global null.
#'
#' 6) **Multiplicity.** For each `rank`, apply BH to `p_perm` → `p_adj_rank`.
#'
#' **Notes.**
#' - If a stratum doesn’t contain exactly two group levels, it’s skipped.
#' - Progress can be printed every `progress_perm_every` permutations.
#'
#' @param x A `phip_data` (recommended; must carry `peptide_library` when
#'   `rank_cols` include non-`peptide_id`) or a data.frame with the necessary
#'   columns.
#' @param rank_cols Character vector of rank columns (e.g. `"peptide_id"`,
#'   `"species"`, ...). Non-peptide ranks are joined from `x$peptide_library`.
#' @param group_cols Character vector of grouping columns that define universes;
#'   pairwise contrasts are constructed within each universe.
#' @param exist_col Name of the 0/1 presence column. Default `"exist"`.
#' @param interaction Logical; if `TRUE`, also create the interaction of the
#'   first two `group_cols`. Ignored if `combine_cols` is set.
#' @param combine_cols Optional length-2 vector; if provided, only that
#'   interaction is used as `group_col`.
#' @param interaction_sep Separator for interaction values. Default `"::"`.
#' @param B_permutations Number of permutations `B`. Default `2000L`.
#' @param seed RNG seed. Default `1L`.
#' @param smooth_eps_num Laplace numerator epsilon. Default `0.5`.
#' @param smooth_eps_den_mult Multiplier for denominator epsilon. Default `2.0`.
#' @param min_max_prev Keep peptides with `max(p1, p2) >=` this. Default `0.0`.
#' @param weight_mode One of `c("equal","se_invvar","n_eff_sqrt")`.
#' @param stat_mode One of `c("diff","asin")`.
#' @param prev_strat One of `c("none","decile")` for prevalence-stratified
#'   combining.
#' @param winsor_z Winsorization threshold for per-peptide z. Default `4.0`.
#' @param parallel Reserved (not used internally). Default `NULL`.
#' @param collect Keep `TRUE`; materialization is not required here.
#' @param register_name Reserved; no-op here.
#' @param progress_perm_every Print a small log every k permutations. Default `1L`.
#' @param rank_feature_keep Optional **named list** `rank -> values to keep`;
#'   filters `(rank, feature)` after pivot.
#'
#' @return A tibble with one row per tested stratum:
#'   `rank, feature, group_col, group1, group2, design, n_subjects_paired,
#'   n_peptides_used, T_obs, p_perm, p_adj_rank, mean_delta, frac_delta_pos,
#'   mean_delta_w, frac_delta_pos_w, category_rank_bh`.
#'
#' @export
ph_prevalence_shift <- function(
    x, rank_cols, group_cols,
    exist_col            = "exist",
    interaction          = FALSE,
    combine_cols         = NULL,
    interaction_sep      = "::",
    # permutation options
    B_permutations       = 2000L,
    seed                 = 1L,
    smooth_eps_num       = 0.5,
    smooth_eps_den_mult  = 2.0,
    min_max_prev         = 0.0,
    weight_mode          = c("equal","se_invvar","n_eff_sqrt"),
    stat_mode            = c("diff","asin"),
    prev_strat           = c("none","decile"),
    winsor_z             = 4.0,
    # parallel
    parallel             = NULL,   # integer cores; NULL => auto (1 by default)
    rank_feature_keep    = NULL,
    peptide_library      = NULL,
    auto_fetch_library   = FALSE,
    # logging
    log                  = FALSE,
    log_file             = "ph_prevalence_shift.log",
    # STREAM/LOG — new arguments
    stream_path          = NULL,          # path to RDS multi-object stream (serialize() appends)
    return_results       = TRUE,          # read back the stream and return tibble
    append_stream        = FALSE          # if TRUE, don't re-init header in stream
) {
  weight_mode <- match.arg(weight_mode)
  stat_mode   <- match.arg(stat_mode)
  prev_strat  <- match.arg(prev_strat)

  # ---- Basic argument validation ---------------------------------------------
  tryCatch({
    chk::chk_character(rank_cols);  chk::chk_true(length(rank_cols) >= 1)
    chk::chk_character(group_cols); chk::chk_true(length(group_cols) >= 1)
    chk::chk_string(exist_col)
    chk::chk_number(B_permutations); chk::chk_true(B_permutations >= 100)
    if (!is.null(rank_feature_keep)) {
      if (!is.list(rank_feature_keep) || is.null(names(rank_feature_keep))) {
        .ph_abort("rank_feature_keep must be a *named* list: rank -> values to keep.")
      }
    }
    if (!is.null(stream_path)) chk::chk_string(stream_path)
  }, error = function(e) .ph_abort("Invalid arguments", bullets = e$message))

  # ---- Prepare long data with required columns --------------------------------
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
        tibble::as_tibble(x) |> dplyr::select(tidyselect::any_of(need))
      } else {
        .ph_abort("Missing required columns", bullets = paste("-", miss))
      }
    }
  }, silent = TRUE)
  if (inherits(df_long, "try-error")) {
    .ph_abort("Could not prepare input data.")
  }

  # Ensure exist_col is integer 0/1
  df_long <- df_long |>
    dplyr::mutate(!!rlang::sym(exist_col) := dplyr::coalesce(as.integer(!!rlang::sym(exist_col)), 0L))

  # ---- Attach peptide library for non-peptide ranks if needed -----------------
  ranks_needing_lib <- setdiff(rank_cols, "peptide_id")
  if (length(ranks_needing_lib)) {
    if (inherits(x, "phip_data")) {
      if (is.null(x$peptide_library)) {
        .ph_abort("Peptide library required for non-peptide ranks.")
      }
      miss_tax <- setdiff(ranks_needing_lib, colnames(x$peptide_library))
      if (length(miss_tax)) {
        .ph_abort("Requested taxonomy/annotation columns not in peptide_library.",
                  bullets = paste("-", miss_tax))
      }
      lib_src <- x$peptide_library %>%
        dplyr::select(tidyselect::any_of(c("peptide_id", ranks_needing_lib))) %>%
        dplyr::distinct()
      if (inherits(df_long, "tbl_sql")) df_long <- df_long %>% dplyr::collect()
      if (inherits(lib_src, "tbl_sql")) lib_src <- lib_src %>% dplyr::collect()
      df_long <- df_long %>% dplyr::left_join(lib_src, by = "peptide_id")
    } else {
      have_all <- all(ranks_needing_lib %in% colnames(df_long))
      if (!have_all) {
        lib_src <- NULL
        if (!is.null(peptide_library)) {
          miss_tax <- setdiff(ranks_needing_lib, colnames(peptide_library))
          if (length(miss_tax)) {
            .ph_abort("Provided `peptide_library` misses required columns:",
                      bullets = paste("-", miss_tax))
          }
          lib_src <- peptide_library %>%
            dplyr::select(tidyselect::any_of(c("peptide_id", ranks_needing_lib))) %>%
            dplyr::distinct()
        } else if (isTRUE(auto_fetch_library)) {
          if (!exists("get_peptide_meta", where = asNamespace("phiper"), inherits = FALSE)) {
            .ph_abort("auto_fetch_library=TRUE but `get_peptide_meta()` not available in {phiper}.")
          }
          .ph_log_info("Retrieving peptide metadata via get_peptide_meta() (auto_fetch_library=TRUE)")
          lib_fetched <- phiper::get_peptide_meta()
          if (inherits(lib_fetched, "tbl_sql") || inherits(lib_fetched, "tbl_lazy")) {
            lib_fetched <- lib_fetched %>% dplyr::collect()
          }
          miss_tax <- setdiff(ranks_needing_lib, colnames(lib_fetched))
          if (length(miss_tax)) {
            .ph_abort("Fetched peptide_meta lacks required columns:",
                      bullets = paste("-", miss_tax))
          }
          lib_src <- lib_fetched %>%
            dplyr::select(tidyselect::any_of(c("peptide_id", ranks_needing_lib))) %>%
            dplyr::distinct()
        } else {
          .ph_abort(
            "Peptide library required for non-peptide ranks.",
            bullets = c(
              "- Pass `peptide_library` with the needed columns, or",
              "- Set `auto_fetch_library = TRUE`, or",
              "- Provide those rank columns already in `x`."
            )
          )
        }
        if (inherits(df_long, "tbl_sql")) df_long <- df_long %>% dplyr::collect()
        if (inherits(lib_src, "tbl_sql")) lib_src <- lib_src %>% dplyr::collect()
        df_long <- df_long %>% dplyr::left_join(lib_src, by = "peptide_id")
      }
    }
  }
  available_ranks <- intersect(rank_cols, colnames(df_long))
  if (!length(available_ranks)) {
    .ph_abort("None of the requested rank_cols are available.")
  }

  # ---- Pivot into (rank, feature) long format --------------------------------
  df_ranked_long <- df_long %>%
    tidyr::pivot_longer(
      cols      = tidyselect::all_of(available_ranks),
      names_to  = "rank",
      values_to = "feature"
    ) %>%
    dplyr::mutate(feature = as.character(feature)) %>%
    dplyr::filter(!is.na(feature))
  if (!"peptide_id" %in% colnames(df_ranked_long)) {
    df_ranked_long <- df_ranked_long %>%
      dplyr::mutate(peptide_id = dplyr::if_else(
        .data$rank == "peptide_id",
        as.character(.data$feature),
        NA_character_
      ))
  }

  # ---- Optional filtering of specific rank-feature pairs ---------------------
  if (!is.null(rank_feature_keep)) {
    clauses <- purrr::imap(rank_feature_keep, function(vals, rk) {
      vals_chr <- as.character(vals)
      rlang::expr((rank == !!rk) & (feature %in% !!vals_chr))
    })
    filtexpr <- if (length(clauses) == 1L) {
      clauses[[1L]]
    } else {
      purrr::reduce(clauses, function(a, b) rlang::expr((!!a) | (!!b)))
    }
    df_ranked_long <- df_ranked_long %>% dplyr::filter(!!filtexpr)
  }

  # ---- Construct grouping universes ------------------------------------------
  make_interaction <- function(tbl, c1, c2, sep) {
    comb_name <- paste(c1, c2, sep = " + ")
    tbl %>%
      dplyr::filter(!is.na(.data[[c1]]) & !is.na(.data[[c2]])) %>%
      dplyr::mutate(
        group_col   = comb_name,
        group_value = paste0(.data[[c1]], sep, .data[[c2]])
      )
  }
  per_column <- df_ranked_long %>%
    tidyr::pivot_longer(
      cols      = tidyselect::all_of(group_cols),
      names_to  = "group_col",
      values_to = "group_value"
    ) %>%
    dplyr::filter(!is.na(group_value))
  if (!is.null(combine_cols)) {
    inter_view <- make_interaction(df_ranked_long, combine_cols[1], combine_cols[2], interaction_sep)
    if (log) {
      .ph_log_info("Grouping universes", bullets = paste0("- ONLY interaction of: ", paste(combine_cols, collapse = " + ")))
    }
  } else if (isTRUE(interaction)) {
    if (length(group_cols) < 2L) {
      .ph_abort("interaction=TRUE needs at least two group_cols.")
    }
    inter_view <- make_interaction(df_ranked_long, group_cols[1], group_cols[2], interaction_sep)
    if (log) {
      .ph_log_info("Grouping universes", bullets = c(
        paste0("- per-column: ", paste(group_cols, collapse = ", ")),
        paste0("- PLUS interaction: ", paste(group_cols[1:2], collapse = " + "))
      ))
    }
  } else {
    inter_view <- NULL
    if (log) {
      .ph_log_info("Grouping universes", bullets = paste0("- per-column only: ", paste(group_cols, collapse = ", ")))
    }
  }

  # ---- Normalize and combine group-wise views --------------------------------
  col_order <- c("sample_id","subject_id","peptide_id", exist_col,
                 "rank","feature","group_col","group_value")
  select_norm <- function(tbl) {
    tbl2 <- tbl
    if (!"peptide_id" %in% colnames(tbl2)) {
      tbl2 <- tbl2 %>%
        dplyr::mutate(peptide_id = dplyr::if_else(
          .data$rank == "peptide_id",
          as.character(.data$feature), NA_character_
        ))
    }
    cols_present <- intersect(col_order, colnames(tbl2))
    tbl2 %>%
      dplyr::select(tidyselect::all_of(cols_present)) %>%
      dplyr::mutate(
        group_col   = as.character(.data$group_col),
        group_value = as.character(.data$group_value)
      )
  }
  per_norm <- per_column %>% select_norm() %>% dplyr::compute()
  if (!is.null(inter_view)) {
    inter_norm <- inter_view %>% select_norm() %>% dplyr::compute()
    gs_view   <- dplyr::union_all(per_norm, inter_norm)
  } else {
    gs_view   <- per_norm
  }

  # ---- Materialize full long data for analysis -------------------------------
  dat_all <- gs_view %>%
    dplyr::select(sample_id, dplyr::any_of("subject_id"), peptide_id,
                  !!rlang::sym(exist_col), rank, feature, group_col, group_value) %>%
    dplyr::rename(exist = !!rlang::sym(exist_col)) %>%
    dplyr::collect()

  # Drop peptides with no presence in any sample
  pep_keep <- dat_all %>%
    dplyr::group_by(peptide_id) %>%
    dplyr::summarise(any_exist = any(exist > 0L), .groups = "drop") %>%
    dplyr::filter(any_exist) %>%
    dplyr::pull(peptide_id)
  dat_all <- dat_all %>% dplyr::filter(peptide_id %in% pep_keep)

  # Define all contrast strata
  strata <- dat_all %>%
    dplyr::distinct(rank, feature, group_col) %>%
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

  n_contrasts <- nrow(strata)

  # ---- STREAM: init stream header --------------------------------------------
  log_to_file <- isTRUE(log) && (!is.null(log_file) && nzchar(log_file))

  if (!is.null(stream_path) && !isTRUE(append_stream)) {
    header <- list(
      type = "ph_prevalence_shift_stream_header",
      version = 1L,
      created = as.character(Sys.time()),
      args = list(
        weight_mode = weight_mode, stat_mode = stat_mode, prev_strat = prev_strat,
        B_permutations = B_permutations, min_max_prev = min_max_prev, winsor_z = winsor_z
      ),
      n_peptides_total = length(unique(dat_all$peptide_id)),
      n_contrasts = n_contrasts
    )
    .stream_init(stream_path, header)
  }

  # ---- Weighted work total precompute (m_eff_i estimate, no permutations) ----
  compute_m_eff_for_stratum <- function(st_row, dat_all_local) {
    gdf <- dat_all_local[
      dat_all_local$rank == st_row$rank &
        dat_all_local$feature == st_row$feature &
        dat_all_local$group_col == st_row$group_col,
      , drop = FALSE
    ]
    lv <- sort(unique(gdf$group_value)); if (length(lv) != 2L) return(0)
    g1 <- lv[1]; g2 <- lv[2]

    prev_tbl <- gdf %>%
      dplyr::group_by(peptide_id, group_value) %>%
      dplyr::summarise(
        x = sum(exist > 0L),
        n = dplyr::n_distinct(if ("subject_id" %in% names(gdf)) subject_id else sample_id),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(names_from = group_value, values_from = c(x,n), values_fill = 0L)

    need <- c(paste0("x_", g1), paste0("x_", g2), paste0("n_", g1), paste0("n_", g2))
    if (!all(need %in% names(prev_tbl))) return(0)

    eps <- smooth_eps_num; den_mult <- smooth_eps_den_mult
    p1 <- (prev_tbl[[paste0("x_", g1)]] + eps) / (prev_tbl[[paste0("n_", g1)]] + den_mult*eps)
    p2 <- (prev_tbl[[paste0("x_", g2)]] + eps) / (prev_tbl[[paste0("n_", g2)]] + den_mult*eps)
    keep <- (pmax(p1, p2) >= min_max_prev)
    if (!any(keep)) return(0)

    p1 <- p1[keep]; p2 <- p2[keep]
    n1 <- prev_tbl[[paste0("n_", g1)]][keep]; n2 <- prev_tbl[[paste0("n_", g2)]][keep]

    comb <- .combine_T_internal(
      p1, p2, n1, n2, winsor_z = winsor_z,
      weight_mode = weight_mode, stat_mode = stat_mode, prev_strat = prev_strat
    )
    m_eff_i <- 1 / pmax(sum(comb$w_norm^2), 1e-12)
    as.numeric(m_eff_i)
  }

  m_eff_est <- vapply(seq_len(n_contrasts), function(i) compute_m_eff_for_stratum(strata[i,], dat_all), numeric(1))
  weight_i  <- as.numeric(B_permutations) * m_eff_est
  work_total <- sum(weight_i)

  # init progress file & log header (with peptides+contrasts+weighted total)
  progress_path <- .progress_file(log_file, stream_path)
  if (isTRUE(log)) .progress_init(progress_path)

  # ----- FAST availability header line ----------------------------------------
  fast_possible <- identical(prev_strat, "none") &&
    weight_mode %in% c("equal","se_invvar","n_eff_sqrt") &&
    stat_mode %in% c("asin","diff") &&
    exists("perm_bitset_T_parallel", mode = "function")
  hdr_fast <- if (fast_possible) "CPP available (bitset+TBB)" else "CPP not available (fallback R)"

  if (isTRUE(log)) {
    hdr_bullets <- c(
      sprintf("total contrasts: %d", n_contrasts),
      sprintf("B per contrast: %d", B_permutations),
      sprintf("peptides after filtering: %d", length(unique(dat_all$peptide_id))),
      sprintf("weighted work total (sum B*m_eff): %.3f", work_total),
      sprintf("mode: %s", if (!is.null(parallel) && parallel > 1L) "parallel (requested)" else "auto/sequential"),
      sprintf("engine: %s", hdr_fast)  # <-- NEW
    )
    if (log_to_file) .ph_log_info_file(log_file, "Starting permutation computations", bullets = hdr_bullets)
    else             .ph_log_info("Starting permutation computations", bullets = hdr_bullets)
  }

  # Prepare result container only if not streaming
  result_rows <- if (is.null(stream_path)) list() else NULL

  # ---- Determine number of workers -------------------------------------------
  plan_changed <- FALSE; original_plan <- NULL
  n_workers <- 1L
  if (!is.null(parallel)) {
    parallel <- as.integer(parallel)
    if (!is.na(parallel) && parallel > 1L) {
      n_workers <- min(parallel, n_contrasts)
    }
  } else {
    if (requireNamespace("future", quietly = TRUE)) {
      if (future::nbrOfWorkers() > 1L) {
        n_workers <- min(future::nbrOfWorkers(), n_contrasts)
      }
    }
  }
  if (n_workers > 1L) {
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
      .ph_warn("Parallel execution requested but {future}/{future.apply} not available; falling back to sequential.")
      n_workers <- 1L
    }
  }

  # ---- Execute: sequential or parallel ---------------------------------------
  if (n_workers == 1L) {
    ## SEQUENTIAL
    done_w <- 0
    for (i in seq_len(n_contrasts)) {
      st <- strata[i, , drop = FALSE]
      dat_uv <- dat_all[
        dat_all$rank == st$rank &
          dat_all$feature == st$feature &
          dat_all$group_col == st$group_col,
        , drop = FALSE]
      res <- .one_contrast_internal(
        dat_uv              = dat_uv,
        weight_mode         = weight_mode,
        stat_mode           = stat_mode,
        prev_strat          = prev_strat,
        B                   = as.integer(B_permutations),
        seed                = as.integer(seed + i - 1L),
        smooth_eps_num      = smooth_eps_num,
        smooth_eps_den_mult = smooth_eps_den_mult,
        min_max_prev        = min_max_prev,
        winsor_z            = winsor_z,
        progress_perm_every = 0L
      )
      if (!is.null(res)) {
        row <- dplyr::bind_cols(st, tibble::as_tibble(res))

        # per-contrast phiper-style detail (with FAST marker)
        if (isTRUE(log)) {
          eng_mark <- if (!is.null(row$engine) && row$engine == "CPP") "[C]" else "[-]"
          msg <- sprintf(
            "%s %s / %s | n_pep=%d m_eff=%.2f b=%d p=%.6g",
            eng_mark, row$rank, row$feature,
            row$n_peptides_used, row$m_eff, row$b, row$p_perm
          )
          if (log_to_file) .ph_log_info_file(log_file, msg) else .ph_log_info(msg)
        }

        # stream or collect
        if (is.null(stream_path)) {
          result_rows[[length(result_rows) + 1L]] <- row
        } else {
          .stream_write_object(stream_path, row)
        }

        # weighted progress update (with FAST marker)
        if (isTRUE(log)) {
          delta_w <- as.numeric(B_permutations) * as.numeric(row$m_eff)
          done_w  <- .progress_add(progress_path, delta_w)
          prog    <- if (work_total > 0) 100 * (done_w / work_total) else 100
          eng_mark <- if (!is.null(row$engine) && row$engine == "CPP") "[C]" else "[-]"
          if (log_to_file) .ph_log_info_file(log_file, sprintf("progress=%.1f%% %s", prog, eng_mark))
          else             .ph_log_info(sprintf("progress=%.1f%% %s", prog, eng_mark))
        }
      }
      if (isTRUE(log) && !log_to_file) {
        cat(sprintf("[prev_shift] %d/%d\n", i, n_contrasts))
        flush.console()
      }
    }

  } else {
    ## PARALLEL via future.apply
    if (requireNamespace("future", quietly = TRUE)) {
      original_plan <- future::plan()  # capture current plan
      if (future::nbrOfWorkers() < 2L || future::nbrOfWorkers() != n_workers) {
        if (.Platform$OS.type == "windows") {
          future::plan(future::multisession, workers = n_workers)
        } else {
          future::plan(future::multicore, workers = n_workers)
        }
        plan_changed <- TRUE
      }
    }
    # split into batches
    k <- min(n_workers, n_contrasts)
    batches <- split(seq_len(n_contrasts), cut(seq_len(n_contrasts), breaks = k, labels = FALSE))
    batch_list <- lapply(batches, function(idx_vec) {
      list(
        data    = dplyr::semi_join(dat_all, strata[idx_vec, ], by = c("rank","feature","group_col")),
        strata  = tibble::as_tibble(strata[idx_vec, , drop = FALSE]),
        offset  = min(idx_vec) - 1L,
        # pass stream/log/progress to workers
        stream_path    = stream_path,
        log_file       = if (isTRUE(log)) log_file else NULL,
        progress_path  = if (isTRUE(log)) .progress_file(log_file, stream_path) else NULL,
        B_permutations = as.integer(B_permutations)
      )
    })

    result_lists <- future.apply::future_lapply(
      X = batch_list,
      FUN = phiper:::.do_batch,         # updated version supports stream/log/progress + FAST marker in logs
      weight_mode = weight_mode, stat_mode = stat_mode, prev_strat = prev_strat,
      B = as.integer(B_permutations), seed_base = as.integer(seed),
      smooth_eps_num = smooth_eps_num, smooth_eps_den_mult = smooth_eps_den_mult,
      min_max_prev = min_max_prev, winsor_z = winsor_z,
      log = isTRUE(log),
      log_path = if (isTRUE(log)) log_file else NULL,
      n_contrasts = n_contrasts,
      future.seed = TRUE,
      future.packages = c("phiper","dplyr","tibble","tidyr","purrr","magrittr","stats","filelock")
    )
    if (is.null(stream_path)) {
      result_rows <- purrr::flatten(result_lists)
    } else {
      result_rows <- NULL  # already streamed by workers
    }

    # master: log final percentage
    if (isTRUE(log) && work_total > 0) {
      conr <- file(.progress_file(log_file, stream_path), open = "rb")
      done_w <- tryCatch(unserialize(conr), error = function(e) 0); close(conr)
      prog <- 100 * (done_w / work_total)
      if (log_to_file) .ph_log_ok_file(log_file, sprintf("progress=%.1f%% (final)", prog))
      else             .ph_log_ok(sprintf("progress=%.1f%% (final)", prog))
    }
  }

  # ---- Final log --------------------------------------------------------------
  if (isTRUE(log)) {
    headline <- "Done (permutation shift)"
    bullets <- c(
      sprintf("ranks: %s", paste(unique(strata$rank), collapse = ", ")),
      sprintf("B: %d", B_permutations),
      sprintf("workers: %d", ifelse(n_workers > 1L, n_workers, 1L))
    )
    if (log_to_file) .ph_log_ok_file(log_file, headline, bullets = bullets) else .ph_log_ok(headline, bullets = bullets)
  }

  # Restore original future plan (shuts down clusters we started)
  if (exists("original_plan") && isTRUE(plan_changed) && !is.null(original_plan)) {
    future::plan(original_plan)
  }

  # ---- Return results: read stream or bind in-memory --------------------------
  res <- if (is.null(stream_path)) {
    dplyr::bind_rows(result_rows)
  } else if (isTRUE(return_results)) {
    .stream_read_all(stream_path) %>% dplyr::filter(!.data$type %in% "ph_prevalence_shift_stream_header")
  } else {
    return(invisible(NULL))
  }

  # ---- p.adjust + final select/arrange (unchanged) ----------------------------
  res <- res %>%
    dplyr::group_by(rank) %>%
    dplyr::mutate(p_adj_rank = p.adjust(p_perm, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      category_rank_bh = dplyr::case_when(
        p_adj_rank < 0.05 ~ "significant (BH, per rank)",
        p_perm     < 0.05 ~ "nominal only",
        TRUE              ~ "not significant"
      )
    ) %>%
    dplyr::select(
      rank, feature, group_col, group1, group2, design,
      n_subjects_paired, n_peptides_used, m_eff,
      T_obs, p_perm, b,
      p_adj_rank,
      mean_delta, frac_delta_pos, mean_delta_w, frac_delta_pos_w,
      category_rank_bh
    ) %>%
    dplyr::arrange(rank, feature, group_col, group1, group2)

  return(res)
}
