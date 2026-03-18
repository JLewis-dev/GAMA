# HELPERS =====================================================================

.GAMA_VERSION <- '0.2.4'

# NCBI configuration

.gama_ncbi_config <- function(api_key = NULL) {
  if (!is.null(api_key) && nzchar(api_key)) {
    rentrez::set_entrez_key(api_key)
  }
  invisible(TRUE)
}

# NCBI defaults

.NCBI_RETMAX <- 999999L

# Null-coalescing operator

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Progress bars

.pb_init <- function(n) {
  utils::txtProgressBar(min = 0, max = n, style = 3)
}

.pb_tick <- function(pb, i) {
  utils::setTxtProgressBar(pb, i)
}

.pb_close <- function(pb) {
  close(pb)
}

# Safe rentrez wrappers

.safe_entrez_search <- function(db, term, retmax = .NCBI_RETMAX, use_history = TRUE, retries = 10, wait = 0.5) {
  for (i in seq_len(retries)) {
    out <- try(
    rentrez::entrez_search(
      db          = db,
      term        = term,
      retmax      = retmax,
      use_history = use_history
    ),
    silent = TRUE
    )
    if (!inherits(out, 'try-error') && !is.null(out)) return(out)
    Sys.sleep(wait * i)
  }
  .gama_stop('NCBI esearch repeatedly failed for ', db, ': ', term)
}

.ncbi_search <- function(db, sp) {
  .safe_entrez_search(
  db          = db,
  term        = paste0(sp, '[Organism]'),
  retmax      = .NCBI_RETMAX,
  use_history = TRUE
  )
}

.safe_entrez_summary <- function(db, ..., retries = 10, wait = 0.5) {
  for (i in seq_len(retries)) {
    out <- try(rentrez::entrez_summary(db = db, ...), silent = TRUE)
    if (!inherits(out, 'try-error') && !is.null(out)) return(out)
    Sys.sleep(wait * i)
  }
  .gama_stop('NCBI esummary repeatedly failed for db = ', db, '.')
}

.search_count <- function(x) {
  as.integer(x$count %||% 0L)
}

.search_ids <- function(x) {
  x$ids %||% character()
}

.search_has_history <- function(x) {
  !is.null(x$web_history)
}

.search_is_truncated <- function(x) {
  .search_count(x) > length(.search_ids(x))
}

.fetch_esummary_batched <- function(db, ids, batch_size = 100) {
  ids <- ids %||% character()
  if (!length(ids)) return(stats::setNames(vector('list', 0), character()))
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  out <- vector('list', length(ids))
  names(out) <- ids
  for (b in batches) {
    res <- .safe_entrez_summary(db, id = b)
    res <- .normalise_esummary_list(res, b)
    out[b] <- res[b]
  }
  out
}

.esummary_uid <- function(x, fallback = NA_character_) {
  uid <- .flatten_to_char(x$uid %||% x$id %||% x$Id)
  if (is.na(uid) || !nzchar(uid)) fallback else uid
}

.normalise_esummary_history_list <- function(SUMS) {
  if (is.null(SUMS)) return(stats::setNames(vector('list', 0), character()))
  out <- if (inherits(SUMS, 'esummary')) list(SUMS) else as.list(SUMS)
  if (!length(out)) return(stats::setNames(vector('list', 0), character()))
  nms <- names(out) %||% rep('', length(out))
  nms <- vapply(seq_along(out), function(i) {
    if (!is.na(nms[[i]]) && nzchar(nms[[i]])) return(nms[[i]])
    .esummary_uid(out[[i]], fallback = as.character(i))
  }, character(1))
  names(out) <- nms
  out
}

.fetch_esummary_history_batched <- function(db, web_history, count, batch_size = 100) {
  count <- as.integer(count %||% 0L)
  if (!count || is.null(web_history)) return(stats::setNames(vector('list', 0), character()))
  starts <- seq.int(0L, count - 1L, by = batch_size)
  out <- list()
  for (start in starts) {
    size <- min(batch_size, count - start)
    res <- .safe_entrez_summary(
    db          = db,
    web_history = web_history,
    retstart    = start,
    retmax      = size
    )
    out <- c(out, .normalise_esummary_history_list(res))
  }
  out
}

.fetch_search_summaries <- function(db, search, batch_size = 100) {
  count <- .search_count(search)
  ids   <- .search_ids(search)
  if (count == 0L) return(stats::setNames(vector('list', 0), character()))
  if (length(ids) && count <= length(ids)) return(.fetch_esummary_batched(db, ids, batch_size = batch_size))
  if (.search_has_history(search)) return(.fetch_esummary_history_batched(db, search$web_history, count, batch_size = batch_size))
  if (length(ids)) return(.fetch_esummary_batched(db, ids, batch_size = batch_size))
  stats::setNames(vector('list', 0), character())
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
  '.env',
  'A',
  'B',
  'BioProject',
  'bioproject',
  'BioSample',
  'biosample',
  'chromatin',
  'class',
  'count',
  'database',
  'denom_total',
  'eff',
  'eff_units',
  'entrez_uid',
  'epigenomic',
  'genomic',
  'geo_linked',
  'geo_linked_denom',
  'geo_prop',
  'gse_ids',
  'gsm_ids',
  'lab',
  'label',
  'label_y',
  'max',
  'med',
  'min',
  'other',
  'prop',
  'q25',
  'q75',
  'S',
  'SRA',
  'score',
  'score_component',
  'segment',
  'species',
  'species_label',
  'status',
  'subclass',
  'transcriptomic',
  'uids_per_unit_max',
  'uids_per_unit_median',
  'uids_per_unit_min',
  'uids_per_unit_q25',
  'uids_per_unit_q75',
  'unknown',
  'x',
  'xmax',
  'xmin',
  'ymax',
  'ymin'
))

# Messaging

.gama_prefix <- function() {
  paste0('GAMA v', .GAMA_VERSION, ' | ')
}

.gama_msg <- function(..., verbose = TRUE) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  message(.gama_prefix(), paste0(..., collapse = ''))
  invisible(NULL)
}

.gama_warn <- function(..., call. = FALSE) {
  warning(.gama_prefix(), paste0(..., collapse = ''), call. = call.)
}

.gama_stop <- function(..., call. = FALSE) {
  stop(.gama_prefix(), paste0(..., collapse = ''), call. = call.)
}

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

.make_query_info <- function(species, dbs) {
  list(
    tool_version   = .GAMA_VERSION,
    query_time_utc = format(as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'),
    databases      = dbs,
    terms          = stats::setNames(as.list(paste0(species, '[Organism]')), species)
  )
}

.as_gdt_table <- function(tbl, results) {
  qi <- attr(results, 'query_info')
  if (is.null(qi)) .gama_warn('No query_info found on results object.')
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
#' Extends the default tibble printing behaviour by prepending a single-line
#' provenance header. If no provenance is attached, a warning message is
#' printed instead.
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
    cat(.gama_prefix(), 'WARNING: no provenance attached.\n', sep = '')
  }
  NextMethod()
}
