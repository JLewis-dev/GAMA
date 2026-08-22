# GAMA 0.4.1

## Reliability
- Retained unlinked operable BioSample records in `interaction_profile` for diagnostics

## Refactoring
- Reorganised internal functions and tests to align with the general package workflow

## Documentation
- Updated roxygen documentation and the GAMA user guide accordingly

## Testing
- Updated existing tests and regenerated examples and fixtures
- Confirmed the test suite passes with the updated fixtures

---

# GAMA 0.4.0

## Features
- Added `report_diagnostics()` with `recoverability` and `records` views for SRA and BioSample availability, SRA and BioSample skew, and interaction summaries

## Reliability
- Standardised BioSample record identity on `input_id` across anatomy profiles, metadata extraction, and interaction analysis, retaining BioSample accessions separately and counting operable records per input record
- Revised SRA and BioSample skew caches to retain active filtered records with missing identifiers while excluding those records from skew calculations
- Added `interaction_profile` and revised `interaction_info$match_report` to report `BioSample`, `operable`, `linked`, and `unknown` counts, with joint unknown modality and anatomy cases counted once

## Plotting
- Excluded missing active identifiers from SRA skew point overlays and removed redundant `(log10)` text from SRA and BioSample skew y-axis labels

## Documentation
- Updated roxygen documentation and the GAMA user guide accordingly

## Testing
- Updated testthat coverage for the new diagnostics and revised cached profile structures
- Regenerated examples and fixtures for the revised output structure
- Confirmed 741 passing tests and `devtools::check()` with 0 errors, 0 warnings, and 0 notes

---

# GAMA 0.3.6

## Reliability
- Expanded the BioSample anatomy ontology using vocabulary identified in legume data and aligned canonical terms against the Plant Ontology<sup>TM</sup>
- Expanded BioSample value normalisation to strip `WT`, `wild type`, and `wildtype` descriptors before anatomy classification
- Added numeric affix stripping to recover anatomy terms from encoded values such as `root1`

## Testing
- Regenerated examples and fixtures for the revised behaviour
- Confirmed the test suite passes with the updated fixtures

---

# GAMA 0.3.5

## Reliability
- Added Entrez UIDs to the cached profiles produced by `summarise_biosample_availability()` to further support diagnostics

## Documentation
- Updated roxygen documentation and the GAMA user guide accordingly

## Testing
- Updated testthat coverage for the revised BioSample cached profile structure
- Regenerated examples and fixtures for the revised output structure
- Confirmed the test suite passes with the updated fixtures

---

# GAMA 0.3.4

## Reliability
- Added raw and normalised SRA strategy values to the cached profile produced by `summarise_sra_availability()` to support diagnostics
- `extract_assembly_metadata()` now retains all records tied at the highest assembly level and N50 when `best = TRUE`
- Expanded the BioSample anatomy ontology using terms from parasitic plant data

## Documentation
- Updated roxygen documentation and the GAMA user guide accordingly

## Testing
- Updated testthat coverage for the revised cached profile and tied assembly handling
- Regenerated examples and fixtures for the revised behaviour
- Confirmed the test suite passes with the updated fixtures

---

# GAMA 0.3.3

## API changes
- Renamed `extract_biosample_metadata()` output columns `tissue_raw` and `tissue_norm` to `value_raw` and `value_norm`

## Reliability
- Added raw and normalised BioSample source values to the cached anatomy profile produced by `summarise_biosample_availability()` to support diagnostics
- Expanded the BioSample anatomy ontology using gymnosperm, pteridophyte, and bryophyte data

## Documentation
- Updated roxygen documentation and the GAMA user guide accordingly

## Testing
- Updated testthat coverage for the revised column names and cached anatomy profile
- Regenerated examples and fixtures for the revised output structure
- Confirmed the test suite passes with the updated fixtures

---

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
  - enforces clean 0–1 axis breaks and two-decimal tick labels

## Testing
- Ran end-to-end availability and plotting workflows to confirm changes do not alter core functionality

---

# GAMA 0.2.0

## Features
- Added SRA skew workflow (`summarise_sra_skew()` and `plot_sra_skew()`) to support BioProject/BioSample-level record aggregation
- Diversity summaries include Inverse Simpson index (`eff`, effective number):
  - low values indicate evidence dominated by a small number of projects/samples
  - high values indicate broader, more balanced support

## Refactoring
- `summarise_sra_availability()` now caches parsed SRA profiles on outputs for reuse by downstream summaries/plots
- Centralised progress-bar handling for more consistent reporting
- Added new Imports dependency `{rlang}`

## API changes
- `extract_sra_metadata()` now returns `entrez_uid` (replacing `sra_id`) and adds `biosample` and `bioproject` columns

## Documentation
- Added/expanded roxygen documentation and examples for the SRA skew workflow
- Updated examples to reflect diversity-aware feasibility assessment workflows

## Testing
- Conducted extended end-to-end tests to validate the new workflow

---

# GAMA 0.1.2

## Refactoring
- Removed redundant code
- Improved internal structure and plotting logic

## Documentation
- Revised and expanded roxygen documentation
- Improved help page formatting and cross-references

## Testing
- Conducted end-to-end workflow tests
- Validated plotting and metadata extraction functions

---

# GAMA 0.1.1

## Refactoring
- Simplified internal scripts and removed redundant code
- Improved internal consistency and readability

## Testing
- Validated behaviour against previous version using smoke tests
- Confirmed compatibility with existing saved result objects
