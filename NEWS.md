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
