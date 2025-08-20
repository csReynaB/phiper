# ==============================================================================
# phiper logging utilities (ASCII only; based on the chk and cli packages)
# ==============================================================================
# ---- user-tweakable globals (set via options() in .onLoad or zzz.R) ----------
# options(
#   phiper.log.verbose   = TRUE,
#   phiper.log.time_fmt  = "%Y-%m-%d %H:%M:%S",
#   phiper.log.width     = getOption("width", 80)
# )

.ph_opt <- function(key,
                    default = NULL) {
  getOption(paste0("phiper.log.", key), default)
}

.ph_now <- function() {
  format(Sys.time(), .ph_opt("time_fmt", "%H:%M:%S"))
}

.ph_base_prefix <- function(level = "INFO") {
  sprintf("[%s] %-5s ", .ph_now(), toupper(level)[1])
}

# wraps the text nicely, regardless of the console width
.ph_wrap <- function(text,
                     prefix) {
  w <- .ph_opt("width", getOption("width", 80))
  # strwrap: 'initial' for first line, 'prefix' for continuations
  strwrap(text, width = w, initial = prefix, prefix = strrep(" ",
                                                             nchar(prefix)))
}

# Compose multi-depth message lines
# currently the maximal supported log depth is 3:
# headline (required), step (optional), bullets (optional chr vec)
.ph_compose_lines <- function(level,
                              headline,
                              step = NULL,
                              bullets = NULL) {
  base  <- .ph_base_prefix(level)
  stepP <- paste0(strrep(" ", nchar(base)), "-> ")
  bullP <- paste0(strrep(" ", nchar(base)), "  - ")

  out <- character(0)
  if (!is.null(headline) && nzchar(headline)) {
    out <- c(out, .ph_wrap(headline, base))
  }
  if (!is.null(step) && nzchar(step)) {
    out <- c(out, .ph_wrap(step, stepP))
  }
  if (length(bullets)) {
    for (b in bullets) {
      if (isTRUE(is.na(b)) || !nzchar(b)) next
      out <- c(out, .ph_wrap(b, bullP))
    }
  }
  out
}

# ---- Public logging helpers --------------------------------------------------
## monitor task progress
.ph_log_info <- function(headline,
                         step = NULL,
                         bullets = NULL,
                         verbose = .ph_opt("verbose", TRUE)) {
  if (!isTRUE(verbose)) return(invisible(character()))
  lines <- .ph_compose_lines("INFO", headline, step, bullets)
  cat(paste0(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}

## monitor task progression
.ph_log_ok <- function(headline,
                       step = NULL,
                       bullets = NULL,
                       verbose = .ph_opt("verbose", TRUE)) {
  if (!isTRUE(verbose)) return(invisible(character()))
  lines <- .ph_compose_lines("OK", headline, step, bullets)
  cat(paste0(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}

# Warnings/errors via chk, but formatted to match the style of the logger
.ph_warn <- function(headline, step = NULL, bullets = NULL, ...) {
  lines <- .ph_compose_lines("WARN", headline, step, bullets)
  msg   <- paste(lines, collapse = "\n")
  if (requireNamespace("chk", quietly = TRUE)) {
    chk::wrn(msg, ...)                 # single string -> respects \n
  } else {
    warning(msg, call. = FALSE, ...)   # fallback if chk not installed
  }
  invisible(lines)
}

.ph_abort <- function(headline, step = NULL, bullets = NULL, ...) {
  lines <- .ph_compose_lines("ERROR", headline, step, bullets)
  msg   <- paste(lines, collapse = "\n")
  if (requireNamespace("chk", quietly = TRUE)) {
    chk::abort_chk(msg, ...)           # single string -> respects \n
  } else {
    stop(msg, call. = FALSE, ...)      # fallback if chk not installed
  }
}

# ---- Back-compaibility (vecmatch) --------------------------------------------
# drop-in replacement for older .msg() funct
.msg <- function(verbose, ...) {
  if (isTRUE(verbose)) .ph_log_info(headline = paste(...))
  invisible(NULL)
}

# original conditional helper to not break down older code
# upgraded to the unified phiper style
.chk_cond <- function(condition,
                      error_message,
                      error   = TRUE,
                      step    = NULL,
                      bullets = NULL,
                      ...) {
  # log nopthing
  if (!isTRUE(condition)) return(invisible(FALSE))

  # print error and abort exec
  if (isTRUE(error)) {
    .ph_abort(headline = error_message, step = step, bullets = bullets, ...)
  } else {
  # print warning and go on
    .ph_warn(headline = error_message, step = step, bullets = bullets, ...)
  }
  invisible(TRUE)
}

# ---- timing helper for sections ----------------------------------------------
## many tasks in phiper can be long/take a while; it was important to have the
## infos on timing - this func wraps a task to get a start/end pair in the same
## style
.ph_with_timing <- function(headline,
                            step = NULL,
                            bullets = NULL,
                            expr,
                            verbose = .ph_opt("verbose", TRUE)) {
  t0 <- proc.time()[["elapsed"]] # register start

  # log infos with the logger
  .ph_log_info(headline = headline,
               step = step,
               bullets = bullets,
               verbose = verbose)

  # return elapsed time
  on.exit({
    dt <- round(proc.time()[["elapsed"]] - t0, 3)
    .ph_log_ok(headline = paste0(headline, " - done"),
               step = sprintf("elapsed: %ss", dt),
               verbose = verbose)
  }, add = TRUE)

  force(expr)
}

# ==============================================================================
# phiper checks + additional helpers (ASCII-only, unified with phiper logger)
# it depends on: .ph_abort(), .ph_warn(), .chk_cond(), word_list(), add_quotes()
# ==============================================================================

# -- check if filename has given extension ------------------------------------
# comes in handy when loading .csv or .parquet; provide filename and vector of
# extensions to check (eg c(".csv", ".parquet"))
.chk_extension <- function(name,
                           x_name,
                           ext_vec) {
  if (is.null(ext_vec) || !length(ext_vec)) return(invisible(TRUE))

  base  <- basename(name %||% "")  # extracting filename from paths
  parts <- strsplit(base, "\\.", fixed = FALSE)[[1]] # names + complex ext

  # taking last ext after . (eg .tar.gz --> .gz)
  ext   <- if (length(parts) > 1L){
    tolower(paste(parts[-1L], collapse = "."))
  } else {
    ""
  }

  norm <- function(x) sub("^\\.+", "", tolower(x)) # normalize ext
  got  <- if (nzchar(ext)) norm(ext) else "<none>"
  ok   <- nzchar(ext) && got %in% norm(ext_vec) # final ext check

  if (!ok) {
    .ph_abort(
      headline = sprintf("Invalid file extension for `%s`.", x_name),
      step     = sprintf("validating path: %s", name),
      bullets  = c(
        sprintf("got: %s", add_quotes(got, 2L)),
        sprintf("allowed: %s",
                word_list(add_quotes(norm(ext_vec), 2L), and_or = "or"))
      )
    )
  }
  invisible(TRUE)
}

# -- check if NULL and replace with default when TRUE (warn in unified style) --
.chk_null_default <- function(x,
                              x_name,
                              method,
                              default) {
  if (is.null(x)) {
    # format the default for print
    fmt <- function(v) {
      if (is.character(v) && length(v) == 1L) return(add_quotes(v, 2L))
      if (is.atomic(v)   && length(v) == 1L) return(as.character(v))
      sprintf("<%s>", paste(class(v), collapse = "/"))
    }

    # generate warning and the replace
    .ph_warn(
      headline = sprintf("Argument `%s` missing; using default.", x_name),
      step     = sprintf("method: %s", add_quotes(method, 2L)),
      bullets  = sprintf("default: %s", fmt(default))
    )
    x <- default
  }
  x
}

# -- validate path to file -----------------------------------------------------
.chk_path <- function(path,
                      arg_name,
                      extension) {

  ## error when path not a string
  .chk_cond(
    !chk::vld_string(path),
    sprintf("`%s` must be a character scalar.", arg_name),
    step    = "path validation",
    bullets = sprintf("got class: %s", paste(class(path), collapse = "/"))
  )

  ## error when path does not exist
  .chk_cond(
    !chk::vld_file(path),
    sprintf("File for `%s` does not exist.", arg_name),
    step    = "path validation",
    bullets = sprintf("path: %s", path)
  )

  # optionally extension check if provided
  if (!missing(extension) && length(extension)) {
    .chk_extension(path,
                   arg_name,
                   extension)
  }

  invisible(TRUE)
}

# -- clean wordlists for message generation ------------------------------------
# for multiple arguments/values
word_list <- function(word_list = NULL,
                      and_or = "and",
                      is_are = FALSE,
                      quotes = FALSE) {

  # Make "a and b" / "a, b, and c"; optionally append "is/are".
  word_list <- setdiff(word_list, c(NA_character_, ""))

  if (is.null(word_list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word_list <- add_quotes(word_list, quotes)

  len_wl <- length(word_list)

  if (len_wl == 1L) {
    out <- word_list
    if (is_are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is.null(and_or) || isFALSE(and_or)) {
    out <- paste(word_list, collapse = ", ")
  } else {
    and_or <- match.arg(and_or, c("and", "or"))
    if (len_wl == 2L) {
      out <- sprintf("%s %s %s", word_list[1L], and_or, word_list[2L])
    } else {
      out <- sprintf("%s, %s %s",
                     paste(word_list[-len_wl], collapse = ", "),
                     and_or, word_list[len_wl])
    }
  }

  if (is_are) out <- sprintf("%s are", out)
  attr(out, "plural") <- TRUE
  out
}

# -- quoting helper (unified error style) --------------------------------------
# define number of quotes you want --> for printing logs/messages/warnings
# or define the quotes itself as a string
add_quotes <- function(x,
                       quotes = 2L) {
  if (isFALSE(quotes)) return(x)
  if (isTRUE(quotes))  quotes <- '"'

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, quotes))
  }

  if (!chk::vld_count(quotes) || quotes > 2) {
    .ph_abort(
      headline = "Invalid `quotes` argument.",
      step     = "formatting add_quotes()",
      bullets  = c(
        "allowed: FALSE, TRUE, 0, 1, 2, or a single-character string",
        sprintf("got class: %s", paste(class(quotes), collapse = "/"))
      )
    )
  }

  if (quotes == 0L) return(x)
  if (quotes == 1L) return(sprintf("'%s'", x))
  sprintf('"%s"', x)
}

# -- not-in operator -----------------------------------------------------------
`%nin%` <- function(x, inx) {
  !(x %in% inx)
}

# -- NULL-coalescing helper ----------------------------------------------------
`%||%` <- function(x, y) if (!is.null(x)) x else y
