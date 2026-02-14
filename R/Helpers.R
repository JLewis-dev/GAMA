# HELPERS =====================================================================

.GAMA_VERSION <- '0.1.2'

# NCBI best practice:

# - set options(ENTREZ_EMAIL=...) for responsible use
# - set rentrez::set_entrez_key() to increase rate limits
# - throttle requests, especially in loops and batch operations
# Configure NCBI access parameters (email and API key)

# NCBI etiquette

.gama_ncbi_config <- function(email = NULL, api_key = NULL) {
  if (!is.null(email) && nzchar(email)) {
    options(ENTREZ_EMAIL = email)
  }
  if (!is.null(api_key) && nzchar(api_key)) {
    rentrez::set_entrez_key(api_key)
  }
  invisible(TRUE)
}

.ncbi_has_key <- function() {
  key <- Sys.getenv('ENTREZ_KEY', unset = NA_character_)
  !is.na(key) && nzchar(key)
}

.ncbi_throttle <- function(keyed = .ncbi_has_key()) {
  Sys.sleep(if (isTRUE(keyed)) 0.12 else 0.35)
}

# Null-coalescing operator

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Safe rentrez wrappers

.safe_search <- purrr::safely(rentrez::entrez_search)

.ncbi_search <- function(db, sp) {
  .ncbi_throttle()
  .safe_search(
    db          = db,
        term        = paste0(sp, '[Organism]'),
    retmax      = 999999,
    use_history = TRUE
  )$result
}

.safe_entrez_summary <- function(db, id, retries = 10, wait = 0.5) {
  for (i in seq_len(retries)) {
    .ncbi_throttle()
    out <- try(rentrez::entrez_summary(db = db, id = id), silent = TRUE)
    if (!inherits(out, 'try-error') && !is.null(out)) return(out)
    Sys.sleep(wait * i)
  }
  stop('NCBI esummary repeatedly failed for: ', paste(id, collapse = ', '))
}

.fetch_esummary_batched <- function(db, ids, batch_size = 100) {
  ids <- ids %||% character()
  if (!length(ids)) return(stats::setNames(vector('list', 0), character()))
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  out <- vector('list', length(ids))
  names(out) <- ids
  for (b in batches) {
    res <- .safe_entrez_summary(db, b)
    res <- .normalise_esummary_list(res, b)
    out[b] <- res[b]
  }
  out
}

.normalise_esummary_list <- function(SUMS, IDS) {
  if (is.null(SUMS)) {
    out <- vector('list', length(IDS))
    names(out) <- IDS
    return(out)
  }
  if (inherits(SUMS, 'esummary')) {
    out <- vector('list', length(IDS))
    names(out) <- IDS
    out[[1]] <- SUMS
    return(out)
  }
  if (!all(names(SUMS) %in% IDS) || length(SUMS) != length(IDS)) {
    tmp <- vector('list', length(IDS))
    names(tmp) <- IDS
    for (id in IDS) tmp[[id]] <- SUMS[[id]] %||% NULL
    return(tmp)
  }
  SUMS
}

# Nomenclature

.shorten_species <- function(x) {
  vapply(
    strsplit(x, ' '),
    function(parts) {
      if (length(parts) < 2) return(x)
      paste0(substr(parts[1], 1, 1), '. ', parts[2])
    },
    character(1)
  )
}

# Globals

utils::globalVariables(c(
'.data',
'A', 'B', 'S', 'SRA',
'chromatin', 'count', 'database', 'epigenomic', 'genomic',
'label', 'other', 'prop', 'score', 'score_component', 'segment',
'species_label', 'status', 'subclass', 'transcriptomic', 'unknown'
))

# Provenance

.flatten_to_char <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.atomic(x)) {
    if (length(x) == 0) return(NA_character_)
    return(as.character(x)[1])
  }
  if (is.list(x)) {
    flat <- unlist(x, recursive = TRUE, use.names = FALSE)
    if (length(flat) == 0) return(NA_character_)
    return(as.character(flat)[1])
  }
  NA_character_
}

.make_query_info <- function(species, dbs, retmax) {
  list(
    tool_version   = .GAMA_VERSION,
    query_time_utc = format(as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'),
    databases      = dbs,
    terms          = stats::setNames(as.list(paste0(species, '[Organism]')), species)
  )
}

.as_gdt_table <- function(tbl, results) {
  qi <- attr(results, 'query_info')
  if (is.null(qi)) warning('No query_info found on results object.')
  attr(tbl, 'query_info') <- qi
  class(tbl) <- c('gdt_tbl', class(tbl))
  tbl
}

#' Print GAMA summary tibble
#'
#' Prints a `gdt_tbl` object with attached query provenance. When present, the
#' `query_info` attribute is displayed above the table, showing the tool
#' version, query timestamp (UTC), and queried databases.
#'
#' @details
#' This method extends the default tibble printing behaviour by prepending a
#' single-line provenance header. If no provenance is attached, a warning
#' message is printed instead.
#'
#' @param x An object of class `gdt_tbl`.
#' @param ... Further arguments passed to the next print method.
#'
#' @return The input object, invisibly.
#'
#' @seealso [summarise_availability()], [summarise_sra_availability()]
#'
#' @examples
#' \dontrun{
#' RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
#' SUMMARY <- summarise_availability(RESULTS)
#' print(SUMMARY)
#' }
#' @export
print.gdt_tbl <- function(x, ...) {
  qi <- attr(x, 'query_info')
  if (!is.null(qi)) {
    cat(
      paste0(
        'GAMA v', qi$tool_version,
        ' | Query time (UTC): ', qi$query_time_utc,
        ' | Databases: ', paste(qi$databases, collapse = ', '),
        '\n'
      )
    )
  } else {
    cat('GAMA | WARNING: no provenance attached\n')
  }
  NextMethod()
}
