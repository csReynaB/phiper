# ==============================================================================
# Pairwise repertoire similarity across categorical time (subjects × subjects)
# ==============================================================================

#' @title Pairwise repertoire similarity across categorical time
#' @description
#' Compute subject×subject similarities between repertoires at pairs of
#' categorical timepoints. Supports Jaccard, Sorensen, Kulczynski (binary),
#' Bray–Curtis (abundance), and Pearson (zero-fill union).
#'
#' Modes with group_col:
#' - "within_group": all timepoint pairs inside each group
#' - "within_time" : at each time, all group pairs present at that time
#' - "mixed"       : label = paste(group,time); all label pairs
#'
#' @param x            <phip_data>; must be longitudinal.
#' @param time_col     Categorical time column name in `x$data_long`.
#' @param group_col    Optional column in `x$data_long` defining groups (or NULL).
#' @param value_col    One of "exist","relative","fold_change".
#' @param methods      Any of: "jaccard","sorensen","kulczynski","braycurtis","pearson".
#' @param compare_mode When `group_col` is set: "within_group","within_time","mixed".
#' @param timepairs    Optional tibble/data.frame with cols `time1`,`time2` (character).
#' @param min_features Integer ≥0; require ≥ this many nonzero features per (subject,time).
#' @param zero_fill_for_pearson Logical; zero-fill union for Pearson (default TRUE).
#' @param write_duckdb Logical; if TRUE, write result to DuckDB via `x$meta$con`.
#' @param register_name Table name to write when `write_duckdb=TRUE`.
#' @param collect      Logical; collect local tibble at the end (default TRUE).
#' @return Long tibble: method,time1,time2,group1,group2,subject1,subject2,similarity
#' @export
ph_compute_pairwise_similarity_cattime <- function(
    x,
    time_col,
    group_col    = NULL,
    value_col    = c("exist","relative","fold_change"),
    methods      = c("jaccard"),
    compare_mode = c("within_group","within_time","mixed"),
    timepairs    = NULL,
    min_features = 1L,
    zero_fill_for_pearson = TRUE,
    write_duckdb = FALSE,
    register_name = NULL,
    collect = TRUE
) {
  stopifnot(inherits(x, "phip_data"))
  if (!isTRUE(x$meta$longitudinal)) {
    .ph_abort("Longitudinal data required.", step = "meta$longitudinal check",
              bullets = "x$meta$longitudinal is FALSE")
  }

  value_col    <- match.arg(value_col)
  compare_mode <- match.arg(compare_mode)
  methods      <- unique(tolower(methods))

  valid_methods <- c("jaccard","sorensen","kulczynski","braycurtis","pearson")
  bad <- setdiff(methods, valid_methods)
  .chk_cond(length(bad) > 0, sprintf("Unknown methods: %s", paste(bad, collapse = ", ")))

  need <- c("subject_id","sample_id","peptide_id", time_col, value_col)
  if (!is.null(group_col)) need <- c(need, group_col)
  miss <- setdiff(need, colnames(x$data_long))
  .chk_cond(length(miss) > 0,
            sprintf("Missing required columns in data_long: %s", paste(miss, collapse = ", ")))

  .ph_with_timing(
    headline = "Pairwise repertoire similarity (categorical time)",
    step     = sprintf("value=%s; methods=%s", value_col, paste(methods, collapse=",")),
    expr = {
      # syms for tidy-eval
      tcol <- rlang::sym(time_col)
      vcol <- rlang::sym(value_col)
      gsym <- if (is.null(group_col)) rlang::sym(".phip__group") else rlang::sym(group_col)

      sid  <- rlang::sym("subject_id")
      samp <- rlang::sym("sample_id")
      pid  <- rlang::sym("peptide_id")

      grp  <- rlang::sym("group")
      tim  <- rlang::sym("time_cat")   # avoid 'time'
      val  <- rlang::sym("value")
      pres <- rlang::sym("present")
      nft  <- rlang::sym("nfeat")

      # --- Standardize to: subject_id | sample_id | peptide_id | group | time_cat | value
      base_tbl <- if (is.null(group_col)) {
        x$data_long |>
          dplyr::select(!!sid, !!samp, !!pid, !!tcol, !!vcol) |>
          dplyr::mutate(!!grp := "all")
      } else {
        x$data_long |>
          dplyr::select(!!sid, !!samp, !!pid, !!gsym, !!tcol, !!vcol) |>
          dplyr::rename(!!grp := !!gsym)
      } |>
        dplyr::transmute(
          !!sid, !!samp, !!pid,
          !!grp,
          !!tim := as.character(!!tcol),
          !!val := !!vcol
        )

      # presence & sparse-profile filter
      base_tbl <- base_tbl |>
        dplyr::mutate(!!pres := dplyr::coalesce(!!val, 0) > 0)

      # robust nfeat without distinct()/count() chain
      prof_sizes <- base_tbl |>
        dplyr::filter(!!pres) |>
        dplyr::group_by(!!sid, !!grp, !!tim) |>
        dplyr::summarise(!!nft := dplyr::n_distinct(!!pid), .groups = "drop")

      base_tbl <- base_tbl |>
        dplyr::inner_join(prof_sizes, by = c("subject_id","group","time_cat")) |>
        dplyr::filter(!!nft >= !!min_features)

      # time levels & pairs
      time_levels <- base_tbl |>
        dplyr::distinct(!!tim) |>
        dplyr::arrange(!!tim) |>
        dplyr::collect() |>
        dplyr::pull(!!tim)

      if (is.null(timepairs)) {
        .chk_cond(length(time_levels) < 2L, "Need at least 2 distinct time levels to form pairs.")
        tp <- utils::combn(time_levels, 2, simplify = FALSE)
        timepairs <- tibble::tibble(
          time1 = vapply(tp, `[[`, character(1), 1),
          time2 = vapply(tp, `[[`, character(1), 2)
        )
      } else {
        timepairs <- tibble::as_tibble(timepairs) |>
          dplyr::transmute(time1 = as.character(time1),
                           time2 = as.character(time2))
      }

      # build the pairing map
      make_pair_map <- function(mode) {
        if (is.null(group_col)) {
          timepairs |> dplyr::mutate(group1 = "all", group2 = "all")

        } else if (mode == "within_group") {
          glevels <- base_tbl |> dplyr::distinct(!!grp) |> dplyr::collect() |> dplyr::pull(!!grp)
          tidyr::crossing(timepairs, group1 = glevels) |>
            dplyr::mutate(group2 = group1)

        } else if (mode == "within_time") {
          groups <- base_tbl |> dplyr::distinct(!!grp) |> dplyr::collect() |> dplyr::pull(!!grp)
          gpairs <- utils::combn(groups, 2, simplify = FALSE)
          gpairs <- tibble::tibble(group1 = vapply(gpairs, `[[`, character(1), 1),
                                   group2 = vapply(gpairs, `[[`, character(1), 2))
          tidyr::crossing(time_cat = time_levels, gpairs) |>
            dplyr::transmute(time1 = time_cat, time2 = time_cat, group1, group2)

        } else { # mixed  ---- build labels in R (after collect)
          seg <- base_tbl |>
            dplyr::distinct(!!grp, !!tim) |>
            dplyr::collect() |>
            dplyr::mutate(label = paste(!!grp, !!tim, sep = "__"))

          labels <- seg$label

          lpairs <- utils::combn(labels, 2, simplify = FALSE)
          lpairs <- tibble::tibble(
            lab1 = vapply(lpairs, `[[`, character(1), 1),
            lab2 = vapply(lpairs, `[[`, character(1), 2)
          )

          split_lab <- function(z) {
            m <- strsplit(z, "__", fixed = TRUE)
            tibble::tibble(
              group    = vapply(m, `[`, character(1), 1),
              time_cat = vapply(m, `[`, character(1), 2)
            )
          }

          dplyr::bind_cols(
            split_lab(lpairs$lab1) |> dplyr::rename(group1 = group, time1 = time_cat),
            split_lab(lpairs$lab2) |> dplyr::rename(group2 = group, time2 = time_cat)
          )
        }
      }
      pair_map <- make_pair_map(if (is.null(group_col)) "nogroup" else compare_mode)

      # helpers to key to requested sets (string 'by' is fine here)
      key_left  <- function(L, pm)  dplyr::inner_join(L, dplyr::distinct(pm, group1, time1),
                                                      by = c("group" = "group1", "time_cat" = "time1"))
      key_right <- function(R, pm)  dplyr::inner_join(R, dplyr::distinct(pm, group2, time2),
                                                      by = c("group" = "group2", "time_cat" = "time2"))

      restrict_to_pairs <- function(df, pm) {
        dplyr::inner_join(df, dplyr::distinct(pm, group1, time1, group2, time2),
                          by = c("group1","time1","group2","time2"))
      }

      # ----- Binary methods ---------------------------------------------------
      compute_binary <- function(seg, pm) {
        L <- seg |>
          dplyr::filter(!!pres) |>
          dplyr::select(!!grp, !!tim, subject1 = !!sid, !!pid)
        R <- seg |>
          dplyr::filter(!!pres) |>
          dplyr::select(!!grp, !!tim, subject2 = !!sid, !!pid)

        Lsz <- L |> dplyr::group_by(!!grp, !!tim, subject1) |>
          dplyr::summarise(n1 = dplyr::n(), .groups = "drop")
        Rsz <- R |> dplyr::group_by(!!grp, !!tim, subject2) |>
          dplyr::summarise(n2 = dplyr::n(), .groups = "drop")

        Lkey <- key_left(L, pm)  |> dplyr::select(group1 = !!grp, time1 = !!tim, subject1, !!pid)
        Rkey <- key_right(R, pm) |> dplyr::select(group2 = !!grp, time2 = !!tim, subject2, !!pid)

        I <- dplyr::inner_join(Lkey, Rkey, by = rlang::set_names("peptide_id")) |>
          dplyr::group_by(group1, time1, subject1, group2, time2, subject2) |>
          dplyr::summarise(n_int = dplyr::n(), .groups = "drop") |>
          restrict_to_pairs(pm)

        out <- I |>
          dplyr::left_join(dplyr::rename(Lsz, group1 = !!grp, time1 = !!tim), by = c("group1","time1","subject1")) |>
          dplyr::left_join(dplyr::rename(Rsz, group2 = !!grp, time2 = !!tim), by = c("group2","time2","subject2")) |>
          dplyr::mutate(
            jaccard    = n_int / pmax(n1 + n2 - n_int, 1L),
            sorensen   = (2 * n_int) / pmax(n1 + n2, 1L),
            kulczynski = 0.5 * (n_int / pmax(n1, 1L) + n_int / pmax(n2, 1L))
          ) |>
          dplyr::select(group1, time1, subject1, group2, time2, subject2,
                        jaccard, sorensen, kulczynski)

        out
      }

      # ----- Abundance methods (Bray–Curtis & Pearson) -----------------------
      compute_abund_corr <- function(seg, pm) {
        L <- seg |> dplyr::select(!!grp, !!tim, subject1 = !!sid, !!pid, value1 = !!val)
        R <- seg |> dplyr::select(!!grp, !!tim, subject2 = !!sid, !!pid, value2 = !!val)

        Lkey <- key_left(L, pm)  |> dplyr::select(group1 = !!grp, time1 = !!tim, subject1, !!pid, value1)
        Rkey <- key_right(R, pm) |> dplyr::select(group2 = !!grp, time2 = !!tim, subject2, !!pid, value2)

        sumL <- Lkey |> dplyr::group_by(group1, time1, subject1) |>
          dplyr::summarise(s1 = sum(dplyr::coalesce(value1, 0)), .groups = "drop")
        sumR <- Rkey |> dplyr::group_by(group2, time2, subject2) |>
          dplyr::summarise(s2 = sum(dplyr::coalesce(value2, 0)), .groups = "drop")

        LR <- dplyr::inner_join(Lkey, Rkey, by = rlang::set_names("peptide_id")) |>
          dplyr::group_by(group1, time1, subject1, group2, time2, subject2) |>
          dplyr::summarise(
            sum_min = sum(pmin(dplyr::coalesce(value1, 0), dplyr::coalesce(value2, 0))),
            sum_xy  = sum(dplyr::coalesce(value1, 0) * dplyr::coalesce(value2, 0)),
            sum_x2  = sum((dplyr::coalesce(value1, 0))^2),
            sum_y2  = sum((dplyr::coalesce(value2, 0))^2),
            .groups = "drop"
          ) |>
          restrict_to_pairs(pm)

        out <- LR |>
          dplyr::left_join(sumL, by = c("group1","time1","subject1")) |>
          dplyr::left_join(sumR, by = c("group2","time2","subject2")) |>
          dplyr::mutate(
            braycurtis = 1 - (2 * sum_min) / pmax(s1 + s2, 1e-12)
          )

        if ("pearson" %in% methods) {
          if (isTRUE(zero_fill_for_pearson)) {
            nzL <- Lkey |> dplyr::filter(dplyr::coalesce(value1, 0) != 0) |>
              dplyr::group_by(group1, time1, subject1) |>
              dplyr::summarise(n_a = dplyr::n(), .groups = "drop")
            nzR <- Rkey |> dplyr::filter(dplyr::coalesce(value2, 0) != 0) |>
              dplyr::group_by(group2, time2, subject2) |>
              dplyr::summarise(n_b = dplyr::n(), .groups = "drop")
            n_int <- Lkey |> dplyr::filter(dplyr::coalesce(value1, 0) != 0) |>
              dplyr::select(group1, time1, subject1, !!pid) |>
              dplyr::inner_join(
                Rkey |> dplyr::filter(dplyr::coalesce(value2, 0) != 0) |>
                  dplyr::select(group2, time2, subject2, !!pid),
                by = rlang::set_names("peptide_id")
              ) |>
              dplyr::group_by(group1, time1, subject1, group2, time2, subject2) |>
              dplyr::summarise(n_int = dplyr::n(), .groups = "drop") |>
              restrict_to_pairs(pm)

            out <- out |>
              dplyr::left_join(nzL, by = c("group1","time1","subject1")) |>
              dplyr::left_join(nzR, by = c("group2","time2","subject2")) |>
              dplyr::left_join(n_int, by = c("group1","time1","subject1","group2","time2","subject2")) |>
              dplyr::mutate(
                n_a   = dplyr::coalesce(n_a, 0L),
                n_b   = dplyr::coalesce(n_b, 0L),
                n_int = dplyr::coalesce(n_int, 0L),
                N     = pmax(n_a + n_b - n_int, 1L),
                mx    = s1 / N,
                my    = s2 / N,
                cov_xy = (sum_xy - N * mx * my),
                var_x  = (sum_x2 - N * (mx^2)),
                var_y  = (sum_y2 - N * (my^2)),
                pearson = dplyr::if_else(var_x <= 0 | var_y <= 0, NA_real_, cov_xy / sqrt(var_x * var_y))
              )
          } else {
            .ph_warn("pearson with zero_fill_for_pearson=FALSE not fully supported; set TRUE for stable results.")
          }
        }

        out |>
          dplyr::select(group1, time1, subject1, group2, time2, subject2,
                        tidyselect::any_of(c("braycurtis","pearson")))
      }

      out_list <- list()
      if (length(intersect(methods, c("jaccard","sorensen","kulczynski"))) > 0) {
        out_list <- c(out_list, list(
          compute_binary(base_tbl, pair_map) |>
            tidyr::pivot_longer(
              cols = tidyselect::any_of(c("jaccard","sorensen","kulczynski")),
              names_to = "method", values_to = "similarity"
            ) |>
            dplyr::filter(.data$method %in% methods)
        ))
      }
      if (length(intersect(methods, c("braycurtis","pearson"))) > 0) {
        out_list <- c(out_list, list(
          compute_abund_corr(base_tbl, pair_map) |>
            tidyr::pivot_longer(
              cols = tidyselect::any_of(c("braycurtis","pearson")),
              names_to = "method", values_to = "similarity"
            ) |>
            dplyr::filter(.data$method %in% methods)
        ))
      }

      res <- if (length(out_list)) dplyr::bind_rows(out_list) else {
        .ph_abort("No methods selected yielded results.")
      }

      res <- res |>
        dplyr::mutate(
          time1 = as.character(.data$time1),
          time2 = as.character(.data$time2),
          group1 = ifelse(is.null(group_col), NA_character_, as.character(.data$group1)),
          group2 = ifelse(is.null(group_col), NA_character_, as.character(.data$group2))
        ) |>
        dplyr::select(.data$method, .data$time1, .data$time2,
                      .data$group1, .data$group2,
                      .data$subject1, .data$subject2, .data$similarity)

      # write via x$meta$con
      if (isTRUE(write_duckdb)) {
        con <- tryCatch(x$meta$con, error = function(e) NULL)
        if (is.null(con)) {
          .ph_abort("write_duckdb=TRUE but x$meta$con is NULL or unavailable.",
                    step = "write_duckdb guard")
        }
        if (is.null(register_name) || !length(register_name) ||
            is.na(register_name)   || !nzchar(register_name)) {
          .ph_abort("write_duckdb=TRUE requires a non-empty `register_name`.",
                    step = "write_duckdb guard")
        }
        res_local <- if (inherits(res, "tbl_lazy")) dplyr::collect(res) else res
        DBI::dbWriteTable(con, register_name, res_local, overwrite = TRUE)
        .ph_log_ok("Wrote pairwise similarity table to DuckDB",
                   bullets = c(paste("table:", register_name),
                               paste("rows:", format(nrow(res_local), big.mark=",")),
                               paste("methods:", paste(methods, collapse=", "))))
      }

      res <- if (isTRUE(collect) && inherits(res, "tbl_lazy")) dplyr::collect(res) else res

      attr(res, "meta") <- list(
        time_levels  = time_levels,
        group_levels = if (!is.null(group_col)) base_tbl |> dplyr::distinct(!!grp) |> dplyr::collect() |> dplyr::pull(!!grp) else "all",
        value_col    = value_col,
        methods      = methods,
        compare_mode = if (is.null(group_col)) "nogroup" else compare_mode
      )
      class(res) <- c("phip_pairwise_similarity", class(res))
      res
    },
    verbose = .ph_opt("verbose", TRUE)
  )
}
