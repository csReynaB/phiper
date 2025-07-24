.rename_to_standard <- function(tbl, colname_map) {
  # colname_map is a named list:
  #   std_name = c("old1", "old2", ...), ...
  # We want: old -> std
  old_to_new <- stats::setNames(names(colname_map), unlist(colname_map))

  dplyr::rename_with(
    tbl,
    function(nm) {
      # vectorised: change only those present in the mapping
      replace <- old_to_new[nm] # returns NA where not mapped
      ifelse(is.na(replace), nm, replace) # keep originals where no match
    }
  )
}

# Rename columns in place - handles both tables and views
.rename_to_standard_inplace <- function(tbl, con, colname_map) {
  ## --- checks -------------------------------------------------------------
  stopifnot(is.character(tbl) && length(tbl) == 1L)
  stopifnot(DBI::dbIsValid(con))
  stopifnot(is.list(colname_map) && length(colname_map) > 0)

  q <- function(x) DBI::dbQuoteString(con, x) # 'string'
  qi <- function(x) DBI::dbQuoteIdentifier(con, x) # "identifier"

  ## --- map old --> new ------------------------------------------------------
  old_to_new <- stats::setNames(names(colname_map), unlist(colname_map))

  ## --- object type --------------------------------------------------------
  meta <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT table_type
         FROM information_schema.tables
        WHERE table_name = %s
          AND table_schema = current_schema",
      q(tbl)
    )
  )
  if (nrow(meta) == 0) {
    stop(sprintf("Object `%s` does not exist in current schema.", tbl),
      call. = FALSE
    )
  }

  is_view <- identical(toupper(meta$table_type[1]), "VIEW")

  ## --- existing columns ---------------------------------------------------
  existing_cols <- DBI::dbGetQuery(
    con,
    sprintf("PRAGMA table_info(%s)", qi(tbl))
  )$name

  old_to_new <- old_to_new[names(old_to_new) %in% existing_cols]
  if (length(old_to_new) == 0L) {
    message("No matching columns found - nothing to rename.")
    return(invisible(tbl))
  }

  ## --- TABLE: simple ALTER TABLE -----------------------------------------
  if (!is_view) {
    for (old in names(old_to_new)) {
      new <- unname(old_to_new[[old]])
      DBI::dbExecute(
        con,
        sprintf(
          "ALTER TABLE %s RENAME COLUMN %s TO %s",
          qi(tbl), qi(old), qi(new)
        )
      )
    }
    return(invisible(tbl))
  }

  ## --- VIEW: recreate with aliases ---------------------------------------
  ## 1. try information_schema.views
  defn <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT view_definition
         FROM information_schema.views
        WHERE table_name = %s
          AND table_schema = current_schema",
      q(tbl)
    )
  )$view_definition

  ## 2. fallback to SHOW CREATE VIEW
  if (is.na(defn) || trimws(defn) == "") {
    defn <- DBI::dbGetQuery(
      con, sprintf("SHOW CREATE VIEW %s", qi(tbl))
    )$sql[1]
  }

  if (is.na(defn) || trimws(defn) == "") {
    stop("Cannot fetch view definition.", call. = FALSE)
  }

  ## 3. strip leading 'CREATE ... AS' and trailing ';'
  defn <- sub(
    pattern = "(?is)^\\s*create\\s+(or\\s+replace\\s+)?view\\s+[^ ]+\\s+as\\s+",
    replacement = "",
    x = defn,
    perl = TRUE
  )
  defn <- sub(";\\s*$", "", defn)

  ## 4. build SELECT list with aliases
  sel_list <- vapply(
    existing_cols,
    function(col) {
      if (col %in% names(old_to_new)) {
        sprintf("%s AS %s", qi(col), qi(old_to_new[[col]]))
      } else {
        qi(col)
      }
    },
    FUN.VALUE = character(1)
  )

  create_sql <- sprintf(
    "CREATE OR REPLACE VIEW %s AS
     SELECT %s
       FROM (%s) AS src;",
    qi(tbl),
    paste(sel_list, collapse = ", "),
    defn
  )
  DBI::dbExecute(con, create_sql)

  invisible(tbl)
}
