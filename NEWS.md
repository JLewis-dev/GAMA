# GAMA 0.3.2

## Plotting
- Updated Assembly and SRA plot labels to use the word 'records' rather than 'accessions' or 'experiments'

## Documentation
- Updated the roxygen, README, and GAMA user guide terminology to align with the revised wording and converted the user guide from PDF to Markdown

## Testing
- Added a shared test helper for retrieving the active GAMA version in synthetic test fixtures
- Confirmed 501 passing tests and `devtools::check()` with 0 errors, 0 warnings, and 0 notes

---

# GAMA 0.3.1

## Features
- Added BioSample replication-skew analysis via `summarise_biosample_skew()` and `plot_biosample_skew()`
- BioSample skew supports metadata interpretation by quantifying record distribution across BioProjects, optionally for a single `anatomy_class`

## Reliability
- Added BioProject provenance to BioSample profile caches, enabling BioSample skew to run from cached summary outputs without re-querying NCBI
- Added `skew_id_recovery` diagnostics to SRA and BioSample skew outputs, making record inclusion, exclusion, and ID recovery explicit
- Refined skew calculations to exclude records with missing skew-unit IDs after filtering
- Patched provenance retention in `extract_assembly_metadata()` when filtering by `species`

## Refactoring
- Standardised SRA and BioSample skew object naming across fixtures, examples, and tests
- Updated fixture-generation and example-generation workflows for new `BIO_SKEW` objects

## Plotting
- Updated the `plot_sra_skew()` axis label to reference SRA records

## Documentation
- Updated roxygen documentation and the GAMA user guide accordingly

## Testing
- Expanded test coverage for BioSample skew, skew diagnostics, validation, and cache provenance
- Confirmed `devtools::check()` with 0 errors, 0 warnings, and 0 notes

---

# GAMA 0.3.0

## Features
- Added a BioSample tissue workflow via `summarise_biosample_availability()`, `plot_biosample_availability()`, and `extract_biosample_metadata()`, completing the core Assembly/SRA/BioSample function family
- BioSample records can now be classified and filtered by anatomy class, anatomy subclass, and anatomy term
- Added an SRA-BioSample interaction workflow via `summarise_interaction()` and `plot_interaction()`
- Interaction summaries link SRA sequencing modality with BioSample-derived anatomy profiles

## Reliability
- Added cached BioSample anatomy profiles to support downstream extraction and SRA-BioSample interaction analysis
- Strengthened cache validation across SRA, BioSample, and interaction workflows
- Improved object-aware validation across the expanded GAMA API
- Extended parameter validation and suggestion logic for modality, anatomy, plotting, and interaction arguments

## Refactoring
- Added internal BioSample record-level parsing, normalisation, ontology-driven classification, and profile-collapsing logic
- Extended cached profile infrastructure for SRA-BioSample interaction analysis
- Refined shared validation helpers for incompatible objects, missing caches, invalid parameters, and close-match suggestions

## Documentation
- Updated roxygen documentation and the GAMA user guide for the new BioSample and interaction workflows

## Testing
- Added a CRAN-level `testthat` suite covering the full GAMA workflow
- Validated end-to-end object fixtures, core formulae, cache integrity, and parameter/error-message logic

---

# GAMA 0.2.9

## Features
- Added `summarise_assembly_availability()` and `plot_assembly_availability()` to complete the NCBI Assembly workflow, bringing it into line with the SRA workflow

## Reliability
- Expanded object-aware error messaging across summary, plotting, and metadata workflows
- Improved parameter validation and suggestion logic for invalid user inputs

## Refactoring
- Added shared helpers for recognised assembly levels, scoring, and `best_n50` selection
- Extended internal GAMA object tagging to support the new Assembly summary workflow

## Documentation
- Updated roxygen documentation and the GAMA user guide, including the addition of a new API map

## Testing
- Tested the new Assembly summary and plotting workflow using real `query_species()` outputs
- Tested compatible and incompatible input paths to confirm the expanded error messaging

---

# GAMA 0.2.8

## Reliability
- Patched SRA modality subclass capitalisation so `other` and `unknown` are treated consistently as lower-case labels
- Added clearer error messaging for invalid `class`, `subclass`, and `unit` argument parameters, including suggestion logic

## API changes
- Removed the `classes` argument from `plot_sra_geo()`
- `plot_sra_geo()` now uses a fixed GEO-oriented modality display set

## Documentation
- Updated roxygen documentation and the GAMA user guide accordingly

## Testing
- Tested canonical, normalised, fuzzy, and unmatched parameter inputs using real `query_species()` workflows

---

# GAMA 0.2.7

## Reliability
- Standardised input validation and error messaging across user-facing summary and plotting functions
- Added object-aware errors for incompatible GAMA and non-GAMA objects, missing columns, and missing cached profiles

## Refactoring
- Added internal object-type tags to support downstream validation
- Centralised validation and cache-check logic in helpers

## Documentation
- Revised roxygen descriptions and updated the user guide accordingly

## Testing
- Ran checks across compatible and incompatible input paths to confirm the new messaging logic and core functionality

---

# GAMA 0.2.6

## Features
- Added synonym-aware querying to `query_species()` via a new `synonyms` argument
- Query results can now be collapsed under canonical species names across Assembly, SRA, and BioSample using unique database record identifiers
- Synonym-collapsed results are returned as one bucket per canonical species without double counting repeated record IDs

## Refactoring
- Added internal helpers for synonym parsing, validation, canonical mapping, and search-result collapse
- Standardised collapsed search outputs so synonym-merged results retain a consistent internal structure for downstream workflows

## Provenance
- Expanded `query_info` to record queried terms and synonym groups in addition to tool version, query timestamp, and database names

## Documentation
- Expanded roxygen `@seealso` cross-references across user-facing functions for clearer upstream/downstream navigation
- Updated `query_species()` roxygen to describe and demonstrate the new `synonyms` argument
- Updated the GAMA user guide to include the new argument

## Testing
- Confirmed functionality with end-to-end workflow tests

---

# GAMA 0.2.5

## API changes
- Renamed `extract_assembly_metadata()` output column `accession` to `entrez_uid` for semantic correctness

## Documentation
- Updated the user guide accordingly

---

# GAMA 0.2.4

## Reliability
- Replaced safe-search path with retrying `entrez_search()` wrappers
- Added `web_history`-aware summary retrieval for large Entrez result sets

## Refactoring
- Switched metadata retrieval to species-local batching
- Simplified NCBI configuration to API key support only
- Removed GAMA request throttling and related legacy code
- Eliminated redundant Assembly summary refetching during metadata extraction
- Removed redundant legacy SRA batch-size code

## Progress reporting
- Updated `query_species()` to tick once per completed Assembly, SRA, and BioSample search

## Documentation
- Updated the GAMA user guide to reflect revised NCBI configuration and history-aware retrieval behaviour

## Testing
- Confirmed functionality with end-to-end workflow tests

---

# GAMA 0.2.3

## Reliability
- Fixed filtering in `extract_sra_metadata()`, `plot_sra_availability()`, and `plot_sra_skew()`
- This bug occurred because `class`, `subclass`, and `species` were used both as function arguments and as metadata column names
- Updated affected filters to use explicit `.env$...` references

## Testing
- Re-tested `extract_sra_metadata()` and confirmed that `class =` and `subclass =` filters now work correctly
- Re-ran SRA availability, GEO overlay, and skew plotting workflows to confirm expected behaviour

---

# GAMA 0.2.2

## Refactoring
- Renamed `plot_sra_geo_availability()` to `plot_sra_geo()`
- GEO linkage fields are now always cached in `summarise_sra_availability()` output via the attached `sra_profile`, regardless of `include_geo`
- `include_geo = TRUE` now acts only as an output visibility option, appending species-level GEO summary columns without changing what is cached

## Documentation
- Updated roxygen documentation and examples to reflect the revised GEO caching and `plot_sra_geo()` workflow
- Added `docs/GAMA_user_guide.pdf` as a comprehensive reference for GAMA functions and methods

## Testing
- Tested end-to-end SRA availability and GEO overlay workflow

---

# GAMA 0.2.1

## Refactoring
- Centralised user-facing messaging via new helpers (`.gama_msg()`, `.gama_warn()`, `.gama_stop()`) for consistent and informative console outputs
- Standardised info/warning/error messaging style (and user-friendly `call. = FALSE` printing where appropriate) to reduce legacy drift
- Updated provenance and print pathways to use the unified messaging style

## Reliability
- Reduced silent drops in availability workflows: summaries retain requested species and emit explicit messages when species have no data

## Plotting
- `plot_sra_geo_availability()` now matches `plot_sra_availability()` styling more closely:
  - prevents GEO-linked fraction label clipping (margin/placement adjustments)
  - enforces clean 0b
