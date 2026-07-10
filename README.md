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
| Arabidopsis thaliana | 378 | 248679 | 249034 | 17.8104 | 24.8478 | 12.4253 | 55.0835 |
| Glycine max | 51 | 46588 | 58528 | 15.8777 | 21.4982 | 10.9773 | 48.3533 |
| Phaseolus vulgaris | 16 | 9633 | 9948 | 14.5109 | 18.3461 | 9.20523 | 42.0622 |
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
| Arabidopsis thaliana | 248679 | 88842 | 131453 | 27255 | 717 | 87 | 325 |
| Glycine max | 46588 | 25509 | 20017 | 1002 | 24 | 0 | 36 |
| Phaseolus vulgaris | 9633 | 6665 | 2765 | 181 | 4 | 0 | 18 |
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
| Arabidopsis thaliana | 5822 | transcriptomic | 1 | 6 | 11 | 18 | 5184 | 229.187 |
| Glycine max | 1306 | transcriptomic | 1 | 1 | 6 | 16 | 801 | 155.337 |
| Phaseolus vulgaris | 146 | transcriptomic | 1 | 5 | 12 | 21 | 200 | 53.6285 |
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
| Arabidopsis thaliana | genomic | leaf | 39095 | 23819.5 | 98.9754 |
| Arabidopsis thaliana | genomic | shoot_meristem | 14 | 346.13 | -17.8521 |
| Arabidopsis thaliana | genomic | stem | 265 | 2815.72 | -48.0693 |
| Arabidopsis thaliana | genomic | root | 452 | 3481.35 | -51.3424 |
| Arabidopsis thaliana | genomic | root_meristem | 0 | 190.651 | -13.8076 |
| Arabidopsis thaliana | genomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | genomic | flower | 3313 | 3328.83 | -0.274417 |
| Arabidopsis thaliana | genomic | fruit | 14 | 116.034 | -9.47224 |
| Arabidopsis thaliana | genomic | seed | 154 | 1144.23 | -29.2739 |
| Arabidopsis thaliana | genomic | whole | 1927 | 10758.3 | -85.1438 |
| Arabidopsis thaliana | genomic | in_vitro | 525 | 390.835 | 6.78647 |
| Arabidopsis thaliana | genomic | other | 1 | 170.271 | -12.9722 |
| Arabidopsis thaliana | genomic | mixed | 4736 | 3312.4 | 24.7353 |
| Arabidopsis thaliana | genomic | unknown | 893 | 1514.69 | -15.9739 |
| Arabidopsis thaliana | transcriptomic | leaf | 27440 | 39872 | -62.2596 |
| Arabidopsis thaliana | transcriptomic | shoot_meristem | 968 | 579.394 | 16.1444 |
| Arabidopsis thaliana | transcriptomic | stem | 7455 | 4713.28 | 39.9357 |
| Arabidopsis thaliana | transcriptomic | root | 9653 | 5827.5 | 50.1126 |
| Arabidopsis thaliana | transcriptomic | root_meristem | 552 | 319.134 | 13.0352 |
| Arabidopsis thaliana | transcriptomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | transcriptomic | flower | 5094 | 5572.19 | -6.40607 |
| Arabidopsis thaliana | transcriptomic | fruit | 313 | 194.232 | 8.52198 |
| Arabidopsis thaliana | transcriptomic | seed | 2782 | 1915.36 | 19.8023 |
| Arabidopsis thaliana | transcriptomic | whole | 24790 | 18008.5 | 50.5341 |
| Arabidopsis thaliana | transcriptomic | in_vitro | 440 | 654.225 | -8.37543 |
| Arabidopsis thaliana | transcriptomic | other | 474 | 285.02 | 11.1938 |
| Arabidopsis thaliana | transcriptomic | mixed | 3221 | 5544.68 | -31.206 |
| Arabidopsis thaliana | transcriptomic | unknown | 2839 | 2535.47 | 6.02805 |
| Arabidopsis thaliana | epigenomic | leaf | 5701 | 8454.5 | -29.9462 |
| Arabidopsis thaliana | epigenomic | shoot_meristem | 71 | 122.855 | -4.67839 |
| Arabidopsis thaliana | epigenomic | stem | 764 | 999.411 | -7.44653 |
| Arabidopsis thaliana | epigenomic | root | 445 | 1235.67 | -22.4928 |
| Arabidopsis thaliana | epigenomic | root_meristem | 28 | 67.6696 | -4.82238 |
| Arabidopsis thaliana | epigenomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | epigenomic | flower | 1650 | 1181.54 | 13.6287 |
| Arabidopsis thaliana | epigenomic | fruit | 26 | 41.1851 | -2.36618 |
| Arabidopsis thaliana | epigenomic | seed | 520 | 406.134 | 5.65012 |
| Arabidopsis thaliana | epigenomic | whole | 5850 | 3818.55 | 32.8743 |
| Arabidopsis thaliana | epigenomic | in_vitro | 217 | 138.723 | 6.64602 |
| Arabidopsis thaliana | epigenomic | other | 43 | 60.436 | -2.24284 |
| Arabidopsis thaliana | epigenomic | mixed | 2084 | 1175.7 | 26.4899 |
| Arabidopsis thaliana | epigenomic | unknown | 841 | 537.624 | 13.0841 |
| Arabidopsis thaliana | chromatin | leaf | 18 | 53.3042 | -4.83554 |
| Arabidopsis thaliana | chromatin | shoot_meristem | 0 | 0.774582 | -0.880103 |
| Arabidopsis thaliana | chromatin | stem | 29 | 6.30111 | 9.04266 |
| Arabidopsis thaliana | chromatin | root | 4 | 7.79069 | -1.35809 |
| Arabidopsis thaliana | chromatin | root_meristem | 0 | 0.426645 | -0.653181 |
| Arabidopsis thaliana | chromatin | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | chromatin | flower | 3 | 7.44937 | -1.63019 |
| Arabidopsis thaliana | chromatin | fruit | 0 | 0.259665 | -0.509573 |
| Arabidopsis thaliana | chromatin | seed | 6 | 2.56061 | 2.14937 |
| Arabidopsis thaliana | chromatin | whole | 46 | 24.0753 | 4.46836 |
| Arabidopsis thaliana | chromatin | in_vitro | 0 | 0.874623 | -0.935213 |
| Arabidopsis thaliana | chromatin | other | 0 | 0.381038 | -0.617283 |
| Arabidopsis thaliana | chromatin | mixed | 4 | 7.41259 | -1.25343 |
| Arabidopsis thaliana | chromatin | unknown | 5 | 3.38962 | 0.874686 |
| Arabidopsis thaliana | other | leaf | 87 | 40.3258 | 7.34998 |
| Arabidopsis thaliana | other | shoot_meristem | 0 | 0.585988 | -0.765498 |
| Arabidopsis thaliana | other | stem | 0 | 4.76693 | -2.18333 |
| Arabidopsis thaliana | other | root | 0 | 5.89382 | -2.42772 |
| Arabidopsis thaliana | other | root_meristem | 0 | 0.322766 | -0.568125 |
| Arabidopsis thaliana | other | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | other | flower | 0 | 5.63561 | -2.37394 |
| Arabidopsis thaliana | other | fruit | 0 | 0.196442 | -0.443218 |
| Arabidopsis thaliana | other | seed | 0 | 1.93715 | -1.39182 |
| Arabidopsis thaliana | other | whole | 0 | 18.2135 | -4.26773 |
| Arabidopsis thaliana | other | in_vitro | 0 | 0.661671 | -0.813432 |
| Arabidopsis thaliana | other | other | 0 | 0.288264 | -0.536902 |
| Arabidopsis thaliana | other | mixed | 0 | 5.60779 | -2.36808 |
| Arabidopsis thaliana | other | unknown | 0 | 2.56432 | -1.60135 |
| Arabidopsis thaliana | mixed | leaf | 77 | 129.321 | -4.60085 |
| Arabidopsis thaliana | mixed | shoot_meristem | 0 | 1.8792 | -1.37084 |
| Arabidopsis thaliana | mixed | stem | 35 | 15.287 | 5.04185 |
| Arabidopsis thaliana | mixed | root | 20 | 18.9009 | 0.252814 |
| Arabidopsis thaliana | mixed | root_meristem | 0 | 1.03508 | -1.01739 |
| Arabidopsis thaliana | mixed | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | mixed | flower | 29 | 18.0728 | 2.57037 |
| Arabidopsis thaliana | mixed | fruit | 0 | 0.62997 | -0.793707 |
| Arabidopsis thaliana | mixed | seed | 7 | 6.21225 | 0.316054 |
| Arabidopsis thaliana | mixed | whole | 61 | 58.4088 | 0.339055 |
| Arabidopsis thaliana | mixed | in_vitro | 7 | 2.12191 | 3.34878 |
| Arabidopsis thaliana | mixed | other | 0 | 0.924432 | -0.961474 |
| Arabidopsis thaliana | mixed | mixed | 24 | 17.9836 | 1.41873 |
| Arabidopsis thaliana | mixed | unknown | 19 | 8.22352 | 3.75793 |
| Arabidopsis thaliana | unknown | leaf | 46 | 95.0205 | -5.02885 |
| Arabidopsis thaliana | unknown | shoot_meristem | 0 | 1.38078 | -1.17506 |
| Arabidopsis thaliana | unknown | stem | 18 | 11.2324 | 2.01928 |
| Arabidopsis thaliana | unknown | root | 17 | 13.8877 | 0.835139 |
| Arabidopsis thaliana | unknown | root_meristem | 0 | 0.760541 | -0.87209 |
| Arabidopsis thaliana | unknown | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | unknown | flower | 38 | 13.2793 | 6.78379 |
| Arabidopsis thaliana | unknown | fruit | 0 | 0.462881 | -0.680354 |
| Arabidopsis thaliana | unknown | seed | 12 | 4.56456 | 3.48022 |
| Arabidopsis thaliana | unknown | whole | 55 | 42.9168 | 1.84445 |
| Arabidopsis thaliana | unknown | in_vitro | 0 | 1.55911 | -1.24864 |
| Arabidopsis thaliana | unknown | other | 0 | 0.679242 | -0.824161 |
| Arabidopsis thaliana | unknown | mixed | 8 | 13.2138 | -1.43429 |
| Arabidopsis thaliana | unknown | unknown | 11 | 6.04237 | 2.01684 |

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
