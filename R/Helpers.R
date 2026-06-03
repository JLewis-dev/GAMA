# HELPERS =====================================================================

.GAMA_VERSION <- '0.3.2'

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

.pb_state <- function(n, enabled = interactive()) {
  n <- as.integer(n %||% 0L)
  if (!isTRUE(enabled) || n <= 0L) return(NULL)
  state <- new.env(parent = emptyenv())
  state$total <- n
  state$current <- 0L
  state$pb <- .pb_init(n)
  state
}

.pb_advance <- function(state, n = 1L) {
  if (is.null(state) || is.null(state$pb)) return(invisible(NULL))
  n <- as.integer(n %||% 1L)
  n <- max(0L, n)
  state$current <- min(state$total, state$current + n)
  .pb_tick(state$pb, state$current)
  invisible(NULL)
}

.pb_close_state <- function(state, complete = TRUE) {
  if (is.null(state) || is.null(state$pb)) return(invisible(NULL))
  if (isTRUE(complete)) .pb_tick(state$pb, state$total)
  .pb_close(state$pb)
  state$pb <- NULL
  invisible(NULL)
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

.search_batch_count <- function(search, batch_size = 100) {
  count <- .search_count(search)
  ids <- .search_ids(search)
  if (count == 0L) return(0L)
  if (length(ids) && count <= length(ids)) return(as.integer(ceiling(length(ids) / batch_size)))
  if (.search_has_history(search)) return(as.integer(ceiling(count / batch_size)))
  if (length(ids)) return(as.integer(ceiling(length(ids) / batch_size)))
  0L
}

.fetch_esummary_batched <- function(db, ids, batch_size = 100, progress_state = NULL, progress_by = 'batch') {
  progress_by <- match.arg(progress_by, c('batch', 'record'))
  ids <- ids %||% character()
  if (!length(ids)) return(stats::setNames(vector('list', 0), character()))
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  out <- vector('list', length(ids))
  names(out) <- ids
  for (b in batches) {
    res <- .safe_entrez_summary(db, id = b)
    res <- .normalise_esummary_list(res, b)
    out[b] <- res[b]
    .pb_advance(progress_state, if (progress_by == 'record') length(b) else 1L)
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

.as_clean_species <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[!is.na(x) & nzchar(x)]
}

.as_synonym_list <- function(synonyms) {
  if (is.null(synonyms)) return(list())
  if (is.list(synonyms)) {
    nms <- names(synonyms) %||% character()
    if (!length(synonyms) || length(nms) != length(synonyms) || any(is.na(nms) | !nzchar(trimws(nms)))) {
      .gama_stop('`synonyms` must be a named list or named character vector.')
    }
    out <- lapply(synonyms, .as_clean_species)
    names(out) <- trimws(nms)
    return(out)
  }
  vals <- as.character(synonyms)
  nms <- names(synonyms) %||% character()
  keep <- !is.na(vals) & nzchar(trimws(vals))
  vals <- trimws(vals[keep])
  nms <- trimws(nms[keep])
  if (!length(vals) || length(nms) != length(vals) || any(is.na(nms) | !nzchar(nms))) {
    .gama_stop('`synonyms` must be a named list or named character vector.')
  }
  canon <- unique(nms)
  out <- stats::setNames(vector('list', length(canon)), canon)
  for (i in seq_along(canon)) out[[i]] <- vals[nms == canon[[i]]]
  out
}

.prepare_synonym_map <- function(species, synonyms = NULL) {
  species <- unique(.as_clean_species(species))
  if (!length(species)) .gama_stop('`species` must contain at least one valid species name.')
  syn <- .as_synonym_list(synonyms)
  if (!length(syn)) {
    query_terms <- stats::setNames(as.list(species), species)
    synonym_map <- stats::setNames(vector('list', length(species)), species)
    return(list(
      species     = species,
      query_terms = query_terms,
      synonyms    = synonym_map
    ))
  }
  syn <- lapply(seq_along(syn), function(i) {
    canonical <- names(syn)[i]
    aliases <- unique(setdiff(.as_clean_species(syn[[i]]), canonical))
    aliases
  })
  names(syn) <- names(.as_synonym_list(synonyms))
  members <- unlist(Map(function(canonical, aliases) c(canonical, aliases), names(syn), syn), use.names = FALSE)
  dup <- unique(members[duplicated(members)])
  if (length(dup) > 0L) {
    .gama_stop('Species names cannot appear in more than one synonym group: ', paste(dup, collapse = ', '), '.')
  }
  active <- vapply(seq_along(syn), function(i) {
    members_i <- c(names(syn)[i], syn[[i]])
    any(members_i %in% species)
  }, logical(1))
  inactive <- names(syn)[!active]
  if (length(inactive) > 0L) {
    .gama_warn('Ignoring synonym groups with no matching entry in `species`: ', paste(inactive, collapse = ', '), '.')
  }
  syn <- syn[active]
  mapped <- species
  for (i in seq_along(syn)) {
    canonical <- names(syn)[i]
    members_i <- c(canonical, syn[[i]])
    mapped[mapped %in% members_i] <- canonical
  }
  species_out <- unique(mapped)
  query_terms <- stats::setNames(as.list(species_out), species_out)
  synonym_map <- stats::setNames(vector('list', length(species_out)), species_out)
  syn_lookup <- syn
  for (sp in species_out) {
    if (sp %in% names(syn_lookup)) {
      query_terms[[sp]] <- unique(c(sp, syn_lookup[[sp]]))
      synonym_map[[sp]] <- syn_lookup[[sp]]
    }
  }
  list(
    species     = species_out,
    query_terms = query_terms,
    synonyms    = synonym_map
  )
}

.normalise_search_result <- function(x) {
  list(
    count       = as.integer(.search_count(x)),
    ids         = .search_ids(x),
    web_history = x$web_history %||% NULL
  )
}

.collapse_searches <- function(searches) {
  if (!length(searches)) {
    return(list(
      count       = 0L,
      ids         = character(),
      web_history = NULL
    ))
  }
  searches <- lapply(searches, .normalise_search_result)
  if (length(searches) == 1L) return(searches[[1]])
  ids <- unique(unlist(lapply(searches, .search_ids), use.names = FALSE))
  list(
    count       = as.integer(length(ids)),
    ids         = ids,
    web_history = NULL
  )
}

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
  'Assembly',
  'B',
  'BioProject',
  'bioproject',
  'best_n50',
  'BioSample',
  'biosample',
  'chromatin',
  'chromosome',
  'class',
  'complete',
  'count',
  'contig',
  'database',
  'denom_total',
  'eff',
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
  'level',
  'level_class',
  'max',
  'med',
  'min',
  'n50',
  'other',
  'prop',
  'q25',
  'q75',
  'S',
  'SRA',
  'score',
  'score_component',
  'scaffold',
  'segment',
  'species',
  'species_label',
  'status',
  'subclass',
  'transcriptomic',
  'unknown',
  'x',
  'xmax',
  'xmin',
  'ymax',
  'ymin',
  '.interaction_term',
  'aerial',
  'anatomy_class',
  'anatomy_subclass',
  'anatomy_term',
  'col_total',
  'end',
  'expected',
  'ground',
  'in_vitro',
  'linked',
  'mixed',
  'modality_class',
  'n',
  'operable',
  'prop_linked',
  'prop_not',
  'reproductive',
  'residual',
  'row_total',
  'term_order',
  'total_n',
  'whole',
  'x_term'
))

# Messaging

.gama_prefix <- function() {
  paste0('GAMA ', .GAMA_VERSION, ' | ')
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

.gama_format_count <- function(x) {
  format(as.numeric(x), big.mark = ',', scientific = FALSE, trim = TRUE)
}

.skew_id_recovery_table <- function(species, unit, class, records, included) {
  species <- as.character(species)
  records <- as.integer(records)
  included <- as.integer(included)
  excluded <- records - included
  id_recovery_prop <- rep(NA_real_, length(records))
  has_records <- records > 0L
  id_recovery_prop[has_records] <- included[has_records] / records[has_records]
  tibble::tibble(
    species = species,
    unit = rep(as.character(unit), length.out = length(species)),
    class = rep(as.character(class), length.out = length(species)),
    records = records,
    included = included,
    excluded = excluded,
    id_recovery_prop = id_recovery_prop
  )
}

# Object validation

.set_gama_object <- function(x, object_name) {
  attr(x, 'gama_object') <- object_name
  x
}

.detect_gama_object <- function(x) {
  tag <- attr(x, 'gama_object', exact = TRUE)
  if (is.character(tag) && length(tag) == 1L && nzchar(tag)) return(tag)
  if (is.list(x) && !inherits(x, c('data.frame', 'tbl_df', 'tbl'))) {
    ok <- length(x) > 0L && all(vapply(x, function(el) {
      is.list(el) && all(c('assembly', 'sra', 'biosample') %in% names(el))
    }, logical(1)))
    if (ok) return('query_species')
  }
  if (inherits(x, 'data.frame')) {
    nms <- names(x)
    if (all(c('species', 'Assembly', 'SRA', 'BioSample', 'A', 'S', 'B', 'score') %in% nms)) {
      return('summarise_availability')
    }
    if (all(c('species', 'Assembly', 'complete', 'chromosome', 'scaffold', 'contig', 'best_n50') %in% nms)) {
      return('summarise_assembly_availability')
    }
    if (all(c('species', 'SRA', 'genomic', 'transcriptomic', 'epigenomic', 'chromatin', 'other', 'unknown') %in% nms)) {
      return('summarise_sra_availability')
    }
    if (all(c('species', 'BioSample', 'operable', 'aerial', 'ground', 'reproductive', 'whole', 'in_vitro', 'other', 'mixed', 'unknown') %in% nms)) {
      return('summarise_biosample_availability')
    }
    if (all(c('species', 'BioSample', 'expected', 'residual') %in% nms) && any(c('class', 'modality_class') %in% nms) && any(c('anatomy_class', 'anatomy_subclass') %in% nms)) {
      return('summarise_interaction')
    }
    if (all(c('species', 'entrez_uid', 'level', 'n50', 'coverage', 'biosample', 'bioproject', 'submitter', 'release_date', 'ftp_path') %in% nms)) {
      return('extract_assembly_metadata')
    }
    if (all(c('species', 'entrez_uid', 'biosample', 'bioproject', 'strategy_raw', 'strategy_norm', 'class', 'subclass', 'geo_linked', 'gse_ids', 'gsm_ids') %in% nms)) {
      return('extract_sra_metadata')
    }
    if (all(c('species', 'entrez_uid', 'biosample', 'bioproject', 'tissue_raw', 'tissue_norm', 'anatomy_class', 'anatomy_subclass', 'anatomy_term') %in% nms)) {
      return('extract_biosample_metadata')
    }
    if (all(c('species', 'BioProject', 'anatomy_class', 'min', 'q25', 'med', 'q75', 'max', 'eff') %in% nms)) {
      return('summarise_biosample_skew')
    }
    if (all(c('species', 'class', 'min', 'q25', 'med', 'q75', 'max', 'eff') %in% nms) && any(c('BioProject', 'BioSample') %in% nms)) {
      return('summarise_sra_skew')
    }
  }
  'incompatible'
}

.gama_input_error <- function(expected, detected = NULL, detail = NULL) {
  detected <- detected %||% 'incompatible'
  if (identical(detected, 'incompatible')) {
    suffix <- 'detected incompatible object.'
  } else if (!is.null(detail) && nzchar(detail)) {
    suffix <- paste0('detected ', detected, ' object ', detail)
  } else {
    suffix <- paste0('detected ', detected, ' object.')
  }
  .gama_stop('Input must be the output of ', expected, '(); ', suffix)
}

.gama_require_output <- function(x, expected, required_cols = NULL) {
  detected <- .detect_gama_object(x)
  if (!identical(detected, expected)) {
    .gama_input_error(expected, detected = detected)
  }
  if (!is.null(required_cols)) {
    missing <- setdiff(required_cols, names(x))
    if (length(missing) > 0L) {
      .gama_input_error(
        expected,
        detected = detected,
        detail = paste0('missing required columns: ', paste(missing, collapse = ', '), '.')
      )
    }
  }
  x
}

.gama_require_cache <- function(x, attr_name, required_cols = NULL, source) {
  detected <- .detect_gama_object(x)
  if (!identical(detected, source)) {
    .gama_input_error(source, detected = detected)
  }
  cache <- attr(x, attr_name, exact = TRUE)
  if (is.null(cache)) {
    .gama_input_error(
      source,
      detected = detected,
      detail = paste0("missing required cache '", attr_name, "'.")
    )
  }
  if (!is.null(required_cols)) {
    missing <- setdiff(required_cols, names(cache))
    if (length(missing) > 0L) {
      .gama_input_error(
        source,
        detected = detected,
        detail = paste0(
          "cache '", attr_name, "' missing required columns: ",
          paste(missing, collapse = ', '),
          '.'
        )
      )
    }
  }
  cache
}

# Parameter validation

.gama_parameter_key <- function(x) {
  x <- trimws(as.character(x))
  x <- tolower(x)
  gsub('[^a-z0-9]+', '', x)
}

.gama_format_values <- function(x) {
  paste0("'", x, "'", collapse = ', ')
}

.gama_rank_parameters <- c('highest', 'lowest', 'A-Z', 'Z-A', 'input')

.gama_sra_skew_unit_parameters <- c('bioproject', 'biosample')

.gama_validate_logical_parameter <- function(x, arg) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) return(x)
  .gama_stop('`', arg, '` must be TRUE or FALSE.')
}

.gama_suggest_variant_parameter <- function(x, variants) {
  if (is.null(variants) || !length(variants)) return(NA_character_)
  key <- .gama_parameter_key(x)
  hit <- unique(unname(variants[names(variants) == key]))
  if (length(hit) == 1L) hit else NA_character_
}

.gama_allow_fuzzy_parameter <- function(x) {
  key <- .gama_parameter_key(x)
  nchar(key) > 4L
}

.gama_suggest_fuzzy_parameter <- function(x, parameters, max_distance = 2L) {
  if (!.gama_allow_fuzzy_parameter(x)) return(NA_character_)
  parameters <- unique(as.character(parameters))
  if (!length(parameters)) return(NA_character_)
  input_key <- .gama_parameter_key(x)
  parameter_keys <- .gama_parameter_key(parameters)
  d <- as.integer(utils::adist(input_key, parameter_keys))
  best <- min(d)
  if (best > max_distance) return(NA_character_)
  hit <- unique(parameters[d == best])
  if (length(hit) == 1L) hit else NA_character_
}

.gama_suggest_parameter <- function(x, parameters, variants = NULL) {
  parameters <- unique(as.character(parameters))
  input_key <- .gama_parameter_key(x)
  parameter_keys <- .gama_parameter_key(parameters)
  normalised_hit <- unique(parameters[parameter_keys == input_key])
  if (length(normalised_hit) == 1L) return(normalised_hit)
  variant_hit <- .gama_suggest_variant_parameter(x, variants)
  if (!is.na(variant_hit) && variant_hit %in% parameters) return(variant_hit)
  .gama_suggest_fuzzy_parameter(x, parameters)
}

.gama_format_suggestions <- function(invalid, suggestions, pairs = length(invalid) > 1L) {
  suggestions <- unname(suggestions)
  if (!isTRUE(pairs)) return(.gama_format_values(suggestions))
  paste0("'", invalid, "' -> '", suggestions, "'", collapse = '; ')
}

.gama_has_normalised_parameter_match <- function(x, parameters) {
  input_key <- .gama_parameter_key(x)
  parameter_keys <- .gama_parameter_key(parameters)
  input_key %in% parameter_keys
}

.gama_validate_single_parameter <- function(x, arg, parameters) {
  if (length(x) == 1L) return(invisible(NULL))
  has_match <- vapply(
    x,
    .gama_has_normalised_parameter_match,
    parameters = parameters,
    FUN.VALUE = logical(1)
  )
  msg <- paste0('`', arg, '` must be a single value.')
  if (any(!has_match)) {
    msg <- paste0(msg, ' Accepted values are: ', .gama_format_values(parameters), '.')
  }
  .gama_stop(msg)
}

.gama_character_parameter_type_error <- function(arg, multiple = TRUE, allow_null = TRUE) {
  msg <- if (isTRUE(multiple)) {
    paste0('`', arg, '` must be a character vector')
  } else {
    paste0('`', arg, '` must be a single character value')
  }
  if (isTRUE(allow_null)) msg <- paste0(msg, ' or NULL')
  .gama_stop(msg, '.')
}

.gama_validate_parameters <- function(x, arg, parameters, variants = NULL, multiple = TRUE, allow_null = TRUE, show_values = TRUE, no_suggestion_hint = NULL) {
  parameters <- unique(as.character(parameters))
  if (is.null(x)) {
    if (isTRUE(allow_null)) return(invisible(NULL))
    msg <- paste0('Invalid `', arg, '` parameter: NULL.')
    if (isTRUE(show_values)) {
      msg <- paste0(msg, ' Accepted values are: ', .gama_format_values(parameters), '.')
    } else if (!is.null(no_suggestion_hint) && nzchar(no_suggestion_hint)) {
      msg <- paste0(msg, ' ', no_suggestion_hint)
    }
    .gama_stop(msg)
  }
  if (!is.character(x)) .gama_character_parameter_type_error(arg, multiple, allow_null)
  if (!multiple) .gama_validate_single_parameter(x, arg, parameters)
  if (any(is.na(x) | !nzchar(trimws(x)))) .gama_stop('`', arg, '` cannot contain missing or empty values.')
  invalid <- setdiff(unique(x), parameters)
  if (!length(invalid)) return(x)
  suggestions <- vapply(
    invalid,
    .gama_suggest_parameter,
    parameters = parameters,
    variants = variants,
    FUN.VALUE = character(1)
  )
  has_suggestion <- !is.na(suggestions) & nzchar(suggestions)
  if (length(invalid) == 1L && isTRUE(has_suggestion)) {
    msg <- paste0(
      'Invalid `',
      arg,
      '` parameter: ',
      .gama_format_values(invalid),
      '. Did you mean ',
      .gama_format_suggestions(invalid, suggestions, pairs = FALSE),
      '?'
    )
    .gama_stop(msg)
  }
  msg <- paste0(
    'Invalid `',
    arg,
    '` parameter',
    if (length(invalid) > 1L) 's' else '',
    ': ',
    .gama_format_values(invalid),
    '.'
  )
  if (isTRUE(show_values)) {
    msg <- paste0(msg, ' Accepted values are: ', .gama_format_values(parameters), '.')
  } else if (!is.null(no_suggestion_hint) && nzchar(no_suggestion_hint)) {
    msg <- paste0(msg, ' ', no_suggestion_hint)
  }
  .gama_stop(msg)
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

.make_query_info <- function(species, dbs, query_terms = NULL, synonyms = NULL) {
  if (is.null(query_terms)) query_terms <- stats::setNames(as.list(species), species)
  if (is.null(synonyms)) synonyms <- stats::setNames(vector('list', length(species)), species)
  query_terms <- query_terms[species]
  synonyms <- synonyms[species]
  list(
    tool_version   = .GAMA_VERSION,
    query_time_utc = format(as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'),
    databases      = dbs,
    terms          = lapply(query_terms, function(x) paste0(x, '[Organism]')),
    synonyms       = synonyms
  )
}

.as_gdt_table <- function(tbl, results, object_name = NULL) {
  qi <- attr(results, 'query_info', exact = TRUE)
  if (is.null(qi)) .gama_warn('No query_info found on results object.')
  attr(tbl, 'query_info') <- qi
  if (!is.null(object_name)) tbl <- .set_gama_object(tbl, object_name)
  class(tbl) <- c('gdt_tbl', class(tbl))
  tbl
}

#' Print GAMA summary tibble
#'
#' Prints a `gdt_tbl` object with attached query provenance. When present, the
#' `query_info` attribute is displayed above the table, showing the tool
#' version, timestamp (UTC), and databases queried.
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
#' @seealso [query_species()]
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
        'GAMA ', qi$tool_version,
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
