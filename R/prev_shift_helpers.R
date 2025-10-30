#' @keywords internal
.combine_T_internal <- function(p1, p2, n1, n2, winsor_z,
                                weight_mode = c("equal","se_invvar","n_eff_sqrt"),
                                stat_mode   = c("diff","asin"),
                                prev_strat  = c("none","decile")) {
  weight_mode <- match.arg(weight_mode)
  stat_mode   <- match.arg(stat_mode)
  prev_strat  <- match.arg(prev_strat)

  delta <- p2 - p1

  if (identical(stat_mode, "asin")) {
    z <- (asin(sqrt(p2)) - asin(sqrt(p1))) /
      sqrt(1/(4*pmax(n1,1)) + 1/(4*pmax(n2,1)))
  } else {
    se <- sqrt(p1*(1-p1)/pmax(n1,1) + p2*(1-p2)/pmax(n2,1))
    z  <- delta / pmax(se, 1e-12)
  }
  z <- pmax(pmin(z, winsor_z), -winsor_z)

  if (identical(weight_mode, "se_invvar") && identical(stat_mode, "asin")) {
    w <- 1 / sqrt(1/(4*pmax(n1,1)) + 1/(4*pmax(n2,1)))
  } else if (identical(weight_mode, "se_invvar")) {
    se <- sqrt(p1*(1-p1)/pmax(n1,1) + p2*(1-p2)/pmax(n2,1))
    w  <- 1/pmax(se, 1e-6)
  } else if (identical(weight_mode, "n_eff_sqrt")) {
    w  <- sqrt(pmax(n1*p1 + n2*p2, 1e-12))
  } else {
    w  <- rep(1.0, length(z))
  }

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

  # Normalized weights used for m_eff downstream.
  w_norm <- w / pmax(sum(w), 1e-12)
  list(T_obs = as.numeric(T_obs), delta = delta, w_norm = w_norm)
}

#' @keywords internal
.one_contrast_internal <- function(dat_uv, weight_mode, stat_mode, prev_strat, B, seed,
                                   smooth_eps_num, smooth_eps_den_mult,
                                   min_max_prev, winsor_z,
                                   progress_perm_every = 1L) {
  # Ensure we have a local data.frame
  if (inherits(dat_uv, "tbl_sql") || inherits(dat_uv, "tbl_lazy")) {
    dat_uv <- dat_uv %>% dplyr::collect()
  }

  # Two groups required
  lv <- sort(unique(dat_uv$group_value))
  if (length(lv) != 2L) return(NULL)
  g1 <- lv[1]; g2 <- lv[2]

  has_subj  <- "subject_id" %in% names(dat_uv)
  paired_ok <- FALSE
  if (has_subj) {
    subj_g1 <- dat_uv %>% dplyr::filter(group_value == g1) %>% dplyr::distinct(subject_id) %>% dplyr::pull()
    subj_g2 <- dat_uv %>% dplyr::filter(group_value == g2) %>% dplyr::distinct(subject_id) %>% dplyr::pull()
    paired_ok <- length(intersect(subj_g1, subj_g2)) > 0L
  }

  dat_a <- if (paired_ok) dat_uv %>% dplyr::filter(subject_id %in% intersect(subj_g1, subj_g2)) else dat_uv
  design <- if (paired_ok) "paired" else "unpaired"

  # ----------------------------------------------------------------------------
  # FAST ENGINE (Rcpp bitset + TBB) — only when prev_strat == "none".
  # Supports stat_mode: "asin" | "diff"; weight_mode: "equal" | "se_invvar" | "n_eff_sqrt".
  # Falls back to R code for other cases (e.g. prev_strat="decile").
  # Requires: perm_bitset_T_parallel() compiled from src/perm_bitset_parallel.cpp
  # ----------------------------------------------------------------------------
  if (identical(prev_strat, "none") &&
      stat_mode %in% c("asin","diff") &&
      weight_mode %in% c("equal","se_invvar","n_eff_sqrt") &&
      exists("perm_bitset_T_parallel", mode = "function")) {

    if (identical(design, "unpaired")) {
      # subjects per group (use sample_id as subject proxy when no subject_id)
      subj_g1 <- dat_a %>%
        dplyr::filter(group_value == g1) %>%
        dplyr::distinct(sample_id) %>% dplyr::pull()
      subj_g2 <- dat_a %>%
        dplyr::filter(group_value == g2) %>%
        dplyr::distinct(sample_id) %>% dplyr::pull()
      all_ids <- c(subj_g1, subj_g2)
      nA <- length(subj_g1); nB <- length(subj_g2); N <- nA + nB
      if (N == 0L || nA == 0L || nB == 0L) return(NULL)

      id_map <- setNames(seq_len(N), all_ids)

      # Build list of 1-based subject indices with exist==1 for each peptide
      hits_by_peptide <- dat_a %>%
        dplyr::mutate(idx = id_map[.data$sample_id]) %>%
        dplyr::group_by(peptide_id) %>%
        dplyr::summarise(idx_hits = list(as.integer(idx[exist > 0L])), .groups = "drop") %>%
        dplyr::arrange(peptide_id)

      if (nrow(hits_by_peptide) == 0L) return(NULL)

      ker <- perm_bitset_T_parallel(
        paired                   = FALSE,
        # unpaired inputs
        hits_by_peptide_unpaired = hits_by_peptide$idx_hits,
        N_unpaired               = as.integer(N),
        g1_idx_unpaired          = as.integer(seq_len(nA)),  # observed split: first nA
        nA                       = as.integer(nA),
        nB                       = as.integer(nB),
        # paired inputs (unused)
        hits_g1_paired           = list(),
        hits_g2_paired           = list(),
        P                        = 0L,
        # common
        B                        = as.integer(B),
        eps                      = smooth_eps_num,
        den_mult                 = smooth_eps_den_mult,
        min_max_prev             = min_max_prev,
        winsor_z                 = winsor_z,
        stat_mode                = stat_mode,
        weight_mode              = weight_mode,
        seed                     = as.integer(seed)
        # grain -> default 512; expose if needed
      )

      return(list(
        group1 = g1, group2 = g2,
        design = "unpaired",
        n_subjects_paired = NA_integer_,
        n_peptides_used   = as.integer(ker$n_peptides_used),
        m_eff             = as.numeric(ker$m_eff),
        T_obs             = as.numeric(ker$T_obs),
        p_perm            = as.numeric(ker$p_perm),
        b                 = as.integer(ker$b),
        mean_delta        = as.numeric(ker$mean_delta),
        frac_delta_pos    = as.numeric(ker$frac_delta_pos),
        mean_delta_w      = as.numeric(ker$mean_delta_w),
        frac_delta_pos_w  = as.numeric(ker$frac_delta_pos_w),
        engine            = "CPP"
      ))
    }

    if (identical(design, "paired")) {
      # paired subjects: complete intersection & stable order
      subj_g1 <- dat_uv %>% dplyr::filter(group_value == g1) %>% dplyr::distinct(subject_id) %>% dplyr::pull()
      subj_g2 <- dat_uv %>% dplyr::filter(group_value == g2) %>% dplyr::distinct(subject_id) %>% dplyr::pull()
      subj_pair <- sort(intersect(subj_g1, subj_g2))
      P <- length(subj_pair)
      if (P == 0L) return(NULL)

      id_map <- setNames(seq_len(P), subj_pair)

      # Build per-peptide hit lists for both groups in the SAME subject order (1..P)
      hits_g1 <- dat_uv %>%
        dplyr::filter(group_value == g1, subject_id %in% subj_pair) %>%
        dplyr::mutate(idx = id_map[.data$subject_id]) %>%
        dplyr::group_by(peptide_id) %>%
        dplyr::summarise(idx_hits = list(as.integer(idx[exist > 0L])), .groups = "drop") %>%
        dplyr::arrange(peptide_id)

      hits_g2 <- dat_uv %>%
        dplyr::filter(group_value == g2, subject_id %in% subj_pair) %>%
        dplyr::mutate(idx = id_map[.data$subject_id]) %>%
        dplyr::group_by(peptide_id) %>%
        dplyr::summarise(idx_hits = list(as.integer(idx[exist > 0L])), .groups = "drop") %>%
        dplyr::arrange(peptide_id)

      # ensure identical peptide order
      if (!identical(hits_g1$peptide_id, hits_g2$peptide_id)) {
        hits_g2 <- hits_g2 %>% dplyr::semi_join(hits_g1, by = "peptide_id") %>% dplyr::arrange(peptide_id)
        hits_g1 <- hits_g1 %>% dplyr::semi_join(hits_g2, by = "peptide_id") %>% dplyr::arrange(peptide_id)
      }
      if (nrow(hits_g1) == 0L) return(NULL)

      ker <- perm_bitset_T_parallel(
        paired                   = TRUE,
        # unpaired inputs (unused)
        hits_by_peptide_unpaired = list(),
        N_unpaired               = 0L,
        g1_idx_unpaired          = integer(),
        nA                       = 0L, nB = 0L,
        # paired inputs
        hits_g1_paired           = hits_g1$idx_hits,
        hits_g2_paired           = hits_g2$idx_hits,
        P                        = as.integer(P),
        # common
        B                        = as.integer(B),
        eps                      = smooth_eps_num,
        den_mult                 = smooth_eps_den_mult,
        min_max_prev             = min_max_prev,
        winsor_z                 = winsor_z,
        stat_mode                = stat_mode,
        weight_mode              = weight_mode,
        seed                     = as.integer(seed)
      )

      return(list(
        group1 = g1, group2 = g2,
        design = "paired",
        n_subjects_paired = as.integer(P),
        n_peptides_used   = as.integer(ker$n_peptides_used),
        m_eff             = as.numeric(ker$m_eff),
        T_obs             = as.numeric(ker$T_obs),
        p_perm            = as.numeric(ker$p_perm),
        b                 = as.integer(ker$b),
        mean_delta        = as.numeric(ker$mean_delta),
        frac_delta_pos    = as.numeric(ker$frac_delta_pos),
        mean_delta_w      = as.numeric(ker$mean_delta_w),
        frac_delta_pos_w  = as.numeric(ker$frac_delta_pos_w),
        engine            = "CPP"
      ))
    }
  }
  # ----------------------- END FAST ENGINE; FALLBACK TO R ----------------------

  prev_tbl <- dat_a %>%
    dplyr::group_by(peptide_id, group_value) %>%
    dplyr::summarise(
      x = sum(exist > 0L),
      n = dplyr::n_distinct(if ("subject_id" %in% names(dat_a)) subject_id else sample_id),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = group_value, values_from = c(x,n), values_fill = 0L)

  need <- c(paste0("x_", g1), paste0("x_", g2), paste0("n_", g1), paste0("n_", g2))
  if (!all(need %in% names(prev_tbl))) return(NULL)

  eps <- smooth_eps_num; den_mult <- smooth_eps_den_mult
  p1 <- (prev_tbl[[paste0("x_", g1)]] + eps) / (prev_tbl[[paste0("n_", g1)]] + den_mult*eps)
  p2 <- (prev_tbl[[paste0("x_", g2)]] + eps) / (prev_tbl[[paste0("n_", g2)]] + den_mult*eps)

  keep <- (pmax(p1, p2) >= min_max_prev)
  if (!any(keep)) return(NULL)

  p1 <- p1[keep]; p2 <- p2[keep]
  n1 <- prev_tbl[[paste0("n_", g1)]][keep]; n2 <- prev_tbl[[paste0("n_", g2)]][keep]

  comb <- .combine_T_internal(p1, p2, n1, n2, winsor_z, weight_mode, stat_mode, prev_strat)
  T_obs  <- comb$T_obs
  delta  <- comb$delta
  w_norm <- comb$w_norm

  # Moment summaries
  mean_delta        <- mean(delta)
  frac_delta_pos    <- mean(delta > 0)
  mean_delta_w      <- sum(w_norm * delta)
  frac_delta_pos_w  <- sum(w_norm * as.numeric(delta > 0))

  # Effective number of peptides under normalized weights (sum w_i = 1)
  m_eff <- 1 / pmax(sum(w_norm^2), 1e-12)

  set.seed(seed)
  B <- as.integer(B)
  T_perm <- numeric(B)

  if (design == "paired") {
    wide <- dat_a %>%
      dplyr::select(subject_id, peptide_id, group_value, exist) %>%
      tidyr::pivot_wider(names_from = c(group_value), values_from = exist, values_fill = 0L,
                         id_cols = c(subject_id, peptide_id))
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
      .combine_T_internal(p1p[keep_p], p2p[keep_p], n1p[keep_p], n2p[keep_p],
                          winsor_z, weight_mode, stat_mode, prev_strat)$T_obs
    }

    for (b in seq_len(B)) {
      flip <- sample(c(FALSE, TRUE), size = length(subj), replace = TRUE)
      T_perm[b] <- .T_from_assignment(flip)
      if (progress_perm_every > 0L && (b %% progress_perm_every == 0L)) {
        # silent in this internal; outer caller handles logging
        invisible(NULL)
      }
    }
  } else {
    subj_g1 <- dat_a %>% dplyr::filter(group_value == g1) %>% dplyr::distinct(sample_id) %>% dplyr::pull()
    subj_g2 <- dat_a %>% dplyr::filter(group_value == g2) %>% dplyr::distinct(sample_id) %>% dplyr::pull()
    n1s <- length(subj_g1); n2s <- length(subj_g2)
    all_ids <- c(subj_g1, subj_g2)

    obs <- dat_a %>% dplyr::select(sample_id, peptide_id, group_value, exist)

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

      need <- c(paste0("x_", g1), paste0("x_", g2), paste0("n_", g1), paste0("n_", g2))
      if (!all(need %in% names(tmp))) return(0)

      p1p <- (tmp[[paste0("x_", g1)]] + eps) / (tmp[[paste0("n_", g1)]] + den_mult*eps)
      p2p <- (tmp[[paste0("x_", g2)]] + eps) / (tmp[[paste0("n_", g2)]] + den_mult*eps)
      keep_p <- (pmax(p1p, p2p) >= min_max_prev)
      if (!any(keep_p)) return(0)
      n1p  <- tmp[[paste0("n_", g1)]][keep_p]; n2p <- tmp[[paste0("n_", g2)]][keep_p]
      .combine_T_internal(p1p[keep_p], p2p[keep_p], n1p, n2p,
                          winsor_z, weight_mode, stat_mode, prev_strat)$T_obs
    }

    for (b in seq_len(B)) {
      idx1 <- sample(seq_along(all_ids), size = n1s, replace = FALSE)
      T_perm[b] <- .T_from_split(idx1)
      if (progress_perm_every > 0L && (b %% progress_perm_every == 0L)) {
        invisible(NULL)
      }
    }
  }

  # Two-sided permutation p and hits
  b_hits <- sum(abs(T_perm) >= abs(T_obs))
  p_two  <- (1 + b_hits) / (1 + B)

  list(
    group1 = g1, group2 = g2,
    design = design,
    n_subjects_paired = if (design == "paired") length(unique(dat_a$subject_id)) else NA_integer_,
    n_peptides_used   = if (exists("delta")) length(delta) else NA_integer_,
    m_eff             = if (exists("m_eff")) as.numeric(m_eff) else NA_real_,
    T_obs             = as.numeric(T_obs),
    p_perm            = as.numeric(p_two),
    b                 = as.integer(b_hits),
    mean_delta        = if (exists("mean_delta")) mean_delta else NA_real_,
    frac_delta_pos    = if (exists("frac_delta_pos")) frac_delta_pos else NA_real_,
    mean_delta_w      = if (exists("mean_delta_w")) as.numeric(mean_delta_w) else NA_real_,
    frac_delta_pos_w  = if (exists("frac_delta_pos_w")) as.numeric(frac_delta_pos_w) else NA_real_,
    engine            = "R"
  )
}


#' @keywords internal
.do_batch <- function(batch,
                      weight_mode, stat_mode, prev_strat,
                      B, seed_base, smooth_eps_num, smooth_eps_den_mult,
                      min_max_prev, winsor_z,
                      log = FALSE, log_path = NULL, n_contrasts = NA_integer_) {
  # ensure magrittr on workers (prints suppressed)
  suppressPackageStartupMessages(requireNamespace("magrittr"))

  # Unpack batch payload
  dat_b         <- batch$data
  strata_b      <- batch$strata
  offset        <- batch$offset %||% 0L
  stream_path   <- batch$stream_path %||% NULL
  log_file      <- batch$log_file %||% NULL
  progress_path <- batch$progress_path %||% NULL
  B_here        <- batch$B_permutations %||% B

  # Fallback for historical arg name: if log_path (fn arg) is given, prefer that
  file_path <- log_file %||% log_path

  n_strata <- nrow(strata_b)
  out <- if (is.null(stream_path)) vector("list", n_strata) else NULL
  counter <- offset

  for (ii in seq_len(n_strata)) {
    st <- strata_b[ii, , drop = FALSE]
    dat_uv <- dat_b[
      dat_b$rank      == st$rank &
        dat_b$feature   == st$feature &
        dat_b$group_col == st$group_col, , drop = FALSE
    ]

    # --- compute one contrast (use B_here; seed offset by global index) -------
    res <- phiper:::.one_contrast_internal(
      dat_uv,
      weight_mode = weight_mode,
      stat_mode   = stat_mode,
      prev_strat  = prev_strat,
      B           = as.integer(B_here),
      seed        = as.integer(seed_base + offset + ii - 1L),
      smooth_eps_num      = smooth_eps_num,
      smooth_eps_den_mult = smooth_eps_den_mult,
      min_max_prev        = min_max_prev,
      winsor_z            = winsor_z,
      progress_perm_every = 0L
    )

    if (!is.null(res)) {
      row <- dplyr::bind_cols(st, tibble::as_tibble(res))

      # ----- weighted progress increment (whole contrast finished) ------------
      if (isTRUE(log) && !is.null(progress_path)) {
        delta_w <- as.numeric(B_here) * as.numeric(row$m_eff)
        done_w  <- .progress_add(progress_path, delta_w)  # cumulative B*m_eff
        # Log absolute cumulative work; percentage is computed on master
        if (!is.null(file_path)) {
          .ph_log_info_file(file_path,
                            sprintf("progress_weight_done=%.3f (B*m_eff cum)", done_w))
        }
      }

      # ----- per-contrast detail line (with FAST marker) ----------------------
      if (isTRUE(log) && !is.null(file_path)) {
        eng_mark <- if (!is.null(row$engine) && identical(row$engine, "CPP")) "[C]" else "[-]"
        msg <- sprintf(
          "%s %s / %s | n_pep=%s m_eff=%.2f b=%s p=%.6g",
          eng_mark,
          row$rank %||% NA_character_,
          row$feature %||% NA_character_,
          row$n_peptides_used %||% NA_integer_,
          row$m_eff %||% NA_real_,
          row$b %||% NA_integer_,
          row$p_perm %||% NA_real_
        )
        .ph_log_info_file(file_path, msg)
      }

      # ----- stream or collect ------------------------------------------------
      if (is.null(stream_path)) {
        out[[ii]] <- row
      } else {
        .stream_write_object(stream_path, row)
      }
    }

    # legacy simple counter line (kept for compatibility)
    if (isTRUE(log) && !is.null(file_path)) {
      counter <- counter + 1L
      cat(sprintf("[prev_shift] %d/%d\n", counter, n_contrasts),
          file = file_path, append = TRUE)
    }
  }

  if (is.null(stream_path)) purrr::compact(out) else NULL
}


# ====== STREAM HELPERS (RDS multi-object stream + locking) ====================

.stream_lock_path <- function(path) paste0(path, ".lock")

.stream_write_object <- function(path, obj) {
  # Use a file lock to make appends from multiple workers safe.
  if (!requireNamespace("filelock", quietly = TRUE)) {
    .ph_abort("`filelock` package required for streaming. Install it or disable streaming.")
  }
  lock_path <- .stream_lock_path(path)
  lk <- filelock::lock(lock_path, timeout = 60000) # up to 60s wait
  on.exit(try(filelock::unlock(lk), silent = TRUE), add = TRUE)

  con <- file(path, open = "ab") # append-binary
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  serialize(obj, con, xdr = FALSE)
  invisible(TRUE)
}

# Read back the whole stream into one tibble (used if return_results=TRUE).
.stream_read_all <- function(path) {
  if (!file.exists(path)) return(tibble::tibble())
  con <- file(path, open = "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  out <- list()
  repeat {
    obj <- try(unserialize(con), silent = TRUE)
    if (inherits(obj, "try-error")) break
    out[[length(out) + 1L]] <- obj
  }
  if (!length(out)) return(tibble::tibble())
  dplyr::bind_rows(out)
}

# Initialize stream with a small header (for sanity/debug).
.stream_init <- function(path, header) {
  if (file.exists(path)) {
    # Overwrite if exists to avoid mixing runs by accident.
    unlink(path)
  }
  .stream_write_object(path, header)
}

# ====== LOG-TO-FILE HELPERS (phiper-style, same layout) =======================

.ph_log_to_file <- function(path, level = "INFO", headline,
                            step = NULL, bullets = NULL) {
  lines <- .ph_compose_lines(level, headline, step, bullets)
  cat(paste0(lines, collapse = "\n"), "\n", sep = "", file = path, append = TRUE)
  invisible(lines)
}

.ph_log_info_file <- function(path, headline, step = NULL, bullets = NULL)
  .ph_log_to_file(path, "INFO", headline, step, bullets)

.ph_log_ok_file <- function(path, headline, step = NULL, bullets = NULL)
  .ph_log_to_file(path, "OK", headline, step, bullets)

.ph_log_warn_file <- function(path, headline, step = NULL, bullets = NULL)
  .ph_log_to_file(path, "WARN", headline, step, bullets)

# ===== Weighted progress (shared across workers) ==============================

.progress_file <- function(log_file = NULL, stream_path = NULL) {
  base <- if (!is.null(log_file) && nzchar(log_file)) log_file else stream_path %||% "ph_prev_shift"
  paste0(base, ".progress.bin")
}
.progress_lock_path <- function(p) paste0(p, ".lock")

.progress_init <- function(path) {
  if (!requireNamespace("filelock", quietly = TRUE)) {
    .ph_abort("`filelock` package required for progress logging. Install it or disable logging.")
  }
  if (file.exists(path)) unlink(path)
  con <- file(path, open = "wb"); on.exit(try(close(con), silent = TRUE), add = TRUE)
  serialize(0, con, xdr = FALSE)  # start at 0
  invisible(TRUE)
}

.progress_add <- function(path, delta) {
  lk <- filelock::lock(.progress_lock_path(path), timeout = 60000)
  on.exit(try(filelock::unlock(lk), silent = TRUE), add = TRUE)

  cur <- 0
  if (file.exists(path)) {
    conr <- file(path, open = "rb"); on.exit(try(close(conr), silent = TRUE), add = TRUE)
    cur <- tryCatch(unserialize(conr), error = function(e) 0)
  }
  newv <- cur + as.numeric(delta)
  conw <- file(path, open = "wb"); on.exit(try(close(conw), silent = TRUE), add = TRUE)
  serialize(newv, conw, xdr = FALSE)
  newv
}
