# GAMA
**Genomic Availability & Metadata Analysis Tool**

<!-- badges: start -->
[![R-CMD-check](https://github.com/JLewis-dev/GAMA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JLewis-dev/GAMA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> &#9888;&#65038; **Development Status: Active**<br>
> GAMA is currently in early development. Interfaces, methods, and outputs may change.<br>
> Users are encouraged to validate results independently and report issues.

Public sequencing archives contain an enormous amount of biological information, but their value can only be realised when the data are findable, accessible, interpretable, and reusable. Without tools for rapidly organising and filtering accession metadata, these records risk becoming a form of digital waste: technically available, but difficult to utilise effectively. GAMA addresses this issue by unifying Assembly, SRA, and BioSample species searches into a single R workflow to produce data availability summaries and ontology-based breakdowns of sequencing modality and sample-source anatomy.

---

## Overview

GAMA:

- Unifies queries across NCBI Assembly, SRA, and BioSample
- Computes a data richness score
- Classifies records by:<br>
  assembly level<br>
  sequencing modality<br>
  sample-source anatomy
- Generates publication-ready visualisations
- Enables targeted extraction of accession metadata

<p align='center'>
  <img src='man/figures/API_map.png' alt='GAMA API map' width='800'>
</p>

See the **[GAMA user guide](docs/GAMA_user_guide.md)** for a comprehensive overview of functions and methods.

---

## Installation

Install the development version from GitHub using `pak`:

``` r
install.packages('pak')
pak::pak('JLewis-dev/GAMA')
```

---

## Quick-start

### 1. Load package

``` r
library(GAMA)
```

### 2. Configure NCBI access to improve rate limits (optional)

``` r
#rentrez::set_entrez_key('YOUR_API_KEY')
```

Uncomment and add your API key if you have one.

### 3. Query NCBI databases using a list of species

``` r
RESULTS <- query_species(c('Vigna angularis', 'Vigna vexillata'))
```

### 4. Summarise data richness

``` r
SUMMARY <- summarise_availability(RESULTS)
print(SUMMARY)
```

### 5. Visualise data richness

``` r
plot_availability(SUMMARY)
```

### 6. Summarise SRA modality composition

``` r
SRA_SUMMARY <- summarise_sra_availability(RESULTS)
print(SRA_SUMMARY)
```

### 7. Visualise SRA modality composition

``` r
plot_sra_availability(SRA_SUMMARY)
```

### 8. Extract filtered Assembly accession metadata

``` r
ASM <- extract_assembly_metadata(RESULTS, species = 'Vigna angularis', best = TRUE)
print(ASM)
```

### 9. Extract filtered SRA accession metadata

``` r
SRA <- extract_sra_metadata(RESULTS, species = 'Vigna vexillata', class = 'genomic')
print(SRA)
```

### 10. Cite

``` r
citation('GAMA')
```

---

## Example outputs

### Data richness

<p align='center'>
  <img src='man/figures/Data_richness_plot.png' alt='GAMA data richness plot' width='800'>
</p>

Data richness provides a weighted overview of sequence availability across species, combining genome assembly quality with transformed SRA and BioSample record counts.

<details>
<summary>Data richness table</summary>

| species | Assembly | SRA | BioSample | A | S | B | score |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis thaliana | 381 | 252433 | 252482 | 17.82 | 24.8778 | 12.4391 | 55.1369 |
| Glycine max | 51 | 47278 | 59194 | 15.8777 | 21.5276 | 10.9886 | 48.394 |
| Phaseolus vulgaris | 17 | 9880 | 10171 | 14.5951 | 18.3967 | 9.22739 | 42.2193 |
| Vigna radiata | 9 | 5629 | 5701 | 13.9703 | 17.2717 | 8.64857 | 39.8906 |

</details>

### Sequencing modality

<p align='center'>
  <img src='man/figures/Modality_plot.png' alt='GAMA modality plot' width='800'>
</p>

Sequencing modality summarises the proportional composition of methodological classes. Similar outputs are available for Assembly level and BioSample anatomy composition.

<details>
<summary>Sequencing modality table</summary>

| species | SRA | genomic | transcriptomic | epigenomic | chromatin | other | unknown |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis thaliana | 252433 | 89924 | 133369 | 27979 | 717 | 87 | 357 |
| Glycine max | 47278 | 25814 | 20395 | 1009 | 24 | 0 | 36 |
| Phaseolus vulgaris | 9880 | 6839 | 2838 | 181 | 4 | 0 | 18 |
| Vigna radiata | 5629 | 4761 | 787 | 54 | 27 | 0 | 0 |

</details>

### Replication skew

<p align='center'>
  <img src='man/figures/Skew_plot.png' alt='GAMA skew plot' width='800'>
</p>

Replication skew reveals how broadly records are distributed across BioProjects or BioSamples.

<details>
<summary>Replication skew table</summary>

| species | BioProject | class | min | q25 | med | q75 | max | eff |
| :--- | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis thaliana | 5915 | transcriptomic | 1 | 6 | 12 | 18.5 | 5184 | 235.73 |
| Glycine max | 1320 | transcriptomic | 1 | 1 | 6 | 16 | 801 | 159.992 |
| Phaseolus vulgaris | 150 | transcriptomic | 1 | 5 | 12 | 21 | 200 | 56.0753 |
| Vigna radiata | 66 | transcriptomic | 1 | 4 | 8 | 14.25 | 65 | 30.3152 |

</details>

### Cross-database interaction

<p align='center'>
  <img src='man/figures/Interaction_plot.png' alt='GAMA interaction plot' width='800'>
</p>

Cross-database interaction links sequencing modality with sample-source anatomy to show which experimental approaches have been applied to different biological materials.

<details>
<summary>Cross-database interaction table</summary>

| species | class | anatomy_subclass | BioSample | expected | residual |
| :--- | :--- | :--- | ---: | ---: | ---: |
| Arabidopsis thaliana | genomic | leaf | 39492 | 24002.3 | 99.9805 |
| Arabidopsis thaliana | genomic | shoot_meristem | 14 | 343.156 | -17.7687 |
| Arabidopsis thaliana | genomic | stem | 265 | 2853.44 | -48.4567 |
| Arabidopsis thaliana | genomic | root | 453 | 3509.44 | -51.5938 |
| Arabidopsis thaliana | genomic | root_meristem | 0 | 192.923 | -13.8897 |
| Arabidopsis thaliana | genomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | genomic | flower | 3313 | 3306.09 | 0.120125 |
| Arabidopsis thaliana | genomic | fruit | 15 | 115.363 | -9.34416 |
| Arabidopsis thaliana | genomic | seed | 155 | 1253.35 | -31.0245 |
| Arabidopsis thaliana | genomic | whole | 1930 | 10825.9 | -85.4982 |
| Arabidopsis thaliana | genomic | in_vitro | 525 | 387.476 | 6.98644 |
| Arabidopsis thaliana | genomic | other | 1 | 168.808 | -12.9156 |
| Arabidopsis thaliana | genomic | mixed | 4736 | 3301.2 | 24.972 |
| Arabidopsis thaliana | genomic | unknown | 897 | 1536.54 | -16.3154 |
| Arabidopsis thaliana | transcriptomic | leaf | 27735 | 40565.2 | -63.7026 |
| Arabidopsis thaliana | transcriptomic | shoot_meristem | 968 | 579.952 | 16.1135 |
| Arabidopsis thaliana | transcriptomic | stem | 7579 | 4822.47 | 39.6944 |
| Arabidopsis thaliana | transcriptomic | root | 9829 | 5931.15 | 50.6123 |
| Arabidopsis thaliana | transcriptomic | root_meristem | 564 | 326.051 | 13.1778 |
| Arabidopsis thaliana | transcriptomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | transcriptomic | flower | 5112 | 5587.47 | -6.3609 |
| Arabidopsis thaliana | transcriptomic | fruit | 313 | 194.969 | 8.453 |
| Arabidopsis thaliana | transcriptomic | seed | 3134 | 2118.23 | 22.0704 |
| Arabidopsis thaliana | transcriptomic | whole | 25200 | 18296.3 | 51.0389 |
| Arabidopsis thaliana | transcriptomic | in_vitro | 440 | 654.855 | -8.39601 |
| Arabidopsis thaliana | transcriptomic | other | 474 | 285.294 | 11.1722 |
| Arabidopsis thaliana | transcriptomic | mixed | 3268 | 5579.21 | -30.9424 |
| Arabidopsis thaliana | transcriptomic | unknown | 2922 | 2596.84 | 6.3808 |
| Arabidopsis thaliana | epigenomic | leaf | 6167 | 8751.8 | -27.6298 |
| Arabidopsis thaliana | epigenomic | shoot_meristem | 71 | 125.122 | -4.83849 |
| Arabidopsis thaliana | epigenomic | stem | 830 | 1040.43 | -6.52379 |
| Arabidopsis thaliana | epigenomic | root | 445 | 1279.62 | -23.3319 |
| Arabidopsis thaliana | epigenomic | root_meristem | 28 | 70.3442 | -5.04871 |
| Arabidopsis thaliana | epigenomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | epigenomic | flower | 1650 | 1205.48 | 12.8031 |
| Arabidopsis thaliana | epigenomic | fruit | 26 | 42.0639 | -2.47684 |
| Arabidopsis thaliana | epigenomic | seed | 532 | 457 | 3.50836 |
| Arabidopsis thaliana | epigenomic | whole | 5926 | 3947.36 | 31.493 |
| Arabidopsis thaliana | epigenomic | in_vitro | 217 | 141.283 | 6.37018 |
| Arabidopsis thaliana | epigenomic | other | 43 | 61.5512 | -2.36458 |
| Arabidopsis thaliana | epigenomic | mixed | 2090 | 1203.69 | 25.5461 |
| Arabidopsis thaliana | epigenomic | unknown | 861 | 560.259 | 12.7057 |
| Arabidopsis thaliana | chromatin | leaf | 18 | 53.2911 | -4.83435 |
| Arabidopsis thaliana | chromatin | shoot_meristem | 0 | 0.761891 | -0.872864 |
| Arabidopsis thaliana | chromatin | stem | 29 | 6.33535 | 9.00459 |
| Arabidopsis thaliana | chromatin | root | 4 | 7.79184 | -1.35841 |
| Arabidopsis thaliana | chromatin | root_meristem | 0 | 0.428338 | -0.654475 |
| Arabidopsis thaliana | chromatin | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | chromatin | flower | 3 | 7.34035 | -1.60201 |
| Arabidopsis thaliana | chromatin | fruit | 0 | 0.256134 | -0.506097 |
| Arabidopsis thaliana | chromatin | seed | 6 | 2.78275 | 1.92863 |
| Arabidopsis thaliana | chromatin | whole | 46 | 24.0361 | 4.47999 |
| Arabidopsis thaliana | chromatin | in_vitro | 0 | 0.860293 | -0.92752 |
| Arabidopsis thaliana | chromatin | other | 0 | 0.374796 | -0.612205 |
| Arabidopsis thaliana | chromatin | mixed | 4 | 7.3295 | -1.22982 |
| Arabidopsis thaliana | chromatin | unknown | 5 | 3.41151 | 0.860027 |
| Arabidopsis thaliana | other | leaf | 87 | 40.3159 | 7.35243 |
| Arabidopsis thaliana | other | shoot_meristem | 0 | 0.576387 | -0.759202 |
| Arabidopsis thaliana | other | stem | 0 | 4.79283 | -2.18925 |
| Arabidopsis thaliana | other | root | 0 | 5.8947 | -2.4279 |
| Arabidopsis thaliana | other | root_meristem | 0 | 0.324047 | -0.569251 |
| Arabidopsis thaliana | other | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | other | flower | 0 | 5.55313 | -2.35651 |
| Arabidopsis thaliana | other | fruit | 0 | 0.193771 | -0.440195 |
| Arabidopsis thaliana | other | seed | 0 | 2.10521 | -1.45093 |
| Arabidopsis thaliana | other | whole | 0 | 18.1838 | -4.26425 |
| Arabidopsis thaliana | other | in_vitro | 0 | 0.650831 | -0.806741 |
| Arabidopsis thaliana | other | other | 0 | 0.283541 | -0.532486 |
| Arabidopsis thaliana | other | mixed | 0 | 5.54492 | -2.35477 |
| Arabidopsis thaliana | other | unknown | 0 | 2.58088 | -1.60651 |
| Arabidopsis thaliana | mixed | leaf | 77 | 130.216 | -4.66346 |
| Arabidopsis thaliana | mixed | shoot_meristem | 0 | 1.86166 | -1.36443 |
| Arabidopsis thaliana | mixed | stem | 35 | 15.4803 | 4.96117 |
| Arabidopsis thaliana | mixed | root | 20 | 19.0392 | 0.220198 |
| Arabidopsis thaliana | mixed | root_meristem | 0 | 1.04663 | -1.02305 |
| Arabidopsis thaliana | mixed | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | mixed | flower | 29 | 17.936 | 2.61246 |
| Arabidopsis thaliana | mixed | fruit | 0 | 0.625859 | -0.791112 |
| Arabidopsis thaliana | mixed | seed | 7 | 6.79958 | 0.0768581 |
| Arabidopsis thaliana | mixed | whole | 63 | 58.7317 | 0.55695 |
| Arabidopsis thaliana | mixed | in_vitro | 7 | 2.10211 | 3.37817 |
| Arabidopsis thaliana | mixed | other | 0 | 0.915805 | -0.956977 |
| Arabidopsis thaliana | mixed | mixed | 24 | 17.9095 | 1.43918 |
| Arabidopsis thaliana | mixed | unknown | 19 | 8.33594 | 3.69356 |
| Arabidopsis thaliana | unknown | leaf | 77 | 109.826 | -3.13232 |
| Arabidopsis thaliana | unknown | shoot_meristem | 0 | 1.57016 | -1.25306 |
| Arabidopsis thaliana | unknown | stem | 18 | 13.0563 | 1.36817 |
| Arabidopsis thaliana | unknown | root | 18 | 16.058 | 0.484632 |
| Arabidopsis thaliana | unknown | root_meristem | 0 | 0.882748 | -0.939547 |
| Arabidopsis thaliana | unknown | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | unknown | flower | 38 | 15.1275 | 5.88071 |
| Arabidopsis thaliana | unknown | fruit | 0 | 0.52786 | -0.726539 |
| Arabidopsis thaliana | unknown | seed | 12 | 5.73488 | 2.61618 |
| Arabidopsis thaliana | unknown | whole | 55 | 49.5353 | 0.776442 |
| Arabidopsis thaliana | unknown | in_vitro | 0 | 1.77295 | -1.33152 |
| Arabidopsis thaliana | unknown | other | 0 | 0.772405 | -0.878866 |
| Arabidopsis thaliana | unknown | mixed | 8 | 15.1051 | -1.82814 |
| Arabidopsis thaliana | unknown | unknown | 11 | 7.03067 | 1.49699 |

</details>

---

## Transparency

GAMA outputs retain query provenance, including tool version, timestamp, database sources, search terms, and synonym handling. Full methodological details can be found in the **[GAMA user guide](docs/GAMA_user_guide.md)**.

### Data richness

The data richness score is defined as:

Score = A + S + B

Where A, S, and B are the transformed contributions of Assembly, SRA, and BioSample record counts.

A = best + ln(1 + total − best), with assemblies weighted as:

- Complete = 10
- Chromosome = 8
- Scaffold = 5
- Contig = 2

Here, best is the maximum-weighted assembly, with ties broken by highest N50, and total is the sum of all assembly weights.

S = 2·ln(1 + SRA)<br>
B = ln(1 + BioSample)

This formulation prioritises high-quality assemblies while incorporating diminishing returns for extensively sampled taxa.

### Ontology-driven classification

SRA records are classified into sequencing modality classes and subclasses using a curated ontology.

<p align='center'>
  <img src='man/figures/SRA_classification.png' alt='GAMA SRA modality classification' width='800'>
</p>

A similar ontology-driven approach is used for classifying BioSample records, converting heterogeneous sample-source metadata into standardised classes, subclasses, and terms.

<p align='center'>
  <img src='man/figures/BioSample_classification.png' alt='GAMA BioSample classification' width='800'>
</p>

---

## Intended use

GAMA is designed for:

- Grant and project scoping
- Identification of under-studied taxa
- Strategic prioritisation of existing datasets

It is particularly suited to investigations of underutilised and non-model plant species.

---

## Limitations

- Dependent on NCBI metadata quality
- Runtime increases with species list size and record count
- Novel protocols may not be captured by the modality ontology
- Some taxon-specific terms may not be captured by the anatomy ontology
- Results should be interpreted cautiously during early development

---

## Licence

See the `LICENSE` files for details.

---
