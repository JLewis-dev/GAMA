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
  - enforces clean 0–1 axis breaks and **two-decimal** tick labels

## Testing
- Ran end-to-end availability and plotting workflows to confirm changes do not alter core functionality

---

# GAMA 0.2.0

## Features
- Added SRA skew workflow (`summarise_sra_skew()` and `plot_sra_skew()`) to support BioProject/BioSample-level record aggregation
- Diversity summaries include **Inverse Simpson index** (`eff`, effective number):
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

---

# GAMA 0.1.0

- Initial public release
