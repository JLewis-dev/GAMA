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
| Arabidopsis thaliana | 378 | 243573 | 243920 | 17.8104 | 24.8064 | 12.4046 | 55.0213 |
| Glycine max | 51 | 46588 | 58528 | 15.8777 | 21.4982 | 10.9773 | 48.3533 |
| Phaseolus vulgaris | 16 | 9629 | 9948 | 14.5109 | 18.3453 | 9.20523 | 42.0614 |
| Vigna radiata | 9 | 5617 | 5701 | 13.9703 | 17.2675 | 8.64857 | 39.8863 |

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
| Arabidopsis thaliana | 243573 | 88839 | 126366 | 27239 | 717 | 87 | 325 |
| Glycine max | 46588 | 25509 | 20017 | 1002 | 24 | 0 | 36 |
| Phaseolus vulgaris | 9629 | 6661 | 2765 | 181 | 4 | 0 | 18 |
| Vigna radiata | 5617 | 4761 | 775 | 54 | 27 | 0 | 0 |

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
| Arabidopsis thaliana | 5805 | transcriptomic | 1 | 6 | 11 | 18 | 5184 | 220.237 |
| Glycine max | 1306 | transcriptomic | 1 | 1 | 6 | 16 | 801 | 155.337 |
| Phaseolus vulgaris | 146 | transcriptomic | 1 | 5 | 12 | 21 | 200 | 53.6285 |
| Vigna radiata | 65 | transcriptomic | 1 | 4 | 8 | 15 | 65 | 29.6064 |

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
| Arabidopsis thaliana | genomic | leaf | 39092 | 22982.3 | 106.265 |
| Arabidopsis thaliana | genomic | shoot_meristem | 14 | 347.625 | -17.8938 |
| Arabidopsis thaliana | genomic | stem | 265 | 2874.79 | -48.6746 |
| Arabidopsis thaliana | genomic | root | 452 | 3582.61 | -52.3033 |
| Arabidopsis thaliana | genomic | root_meristem | 0 | 197.09 | -14.0389 |
| Arabidopsis thaliana | genomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | genomic | flower | 3313 | 3461.98 | -2.53201 |
| Arabidopsis thaliana | genomic | fruit | 14 | 119.953 | -9.67402 |
| Arabidopsis thaliana | genomic | seed | 154 | 1182.88 | -29.9153 |
| Arabidopsis thaliana | genomic | whole | 1927 | 10975.2 | -86.3684 |
| Arabidopsis thaliana | genomic | in_vitro | 525 | 404.034 | 6.01806 |
| Arabidopsis thaliana | genomic | other | 1 | 170.924 | -12.9973 |
| Arabidopsis thaliana | genomic | mixed | 4736 | 3429.36 | 22.3126 |
| Arabidopsis thaliana | genomic | unknown | 893 | 1657.25 | -18.7734 |
| Arabidopsis thaliana | transcriptomic | leaf | 22618 | 36193.2 | -71.3563 |
| Arabidopsis thaliana | transcriptomic | shoot_meristem | 938 | 547.449 | 16.6919 |
| Arabidopsis thaliana | transcriptomic | stem | 7369 | 4527.29 | 42.2338 |
| Arabidopsis thaliana | transcriptomic | root | 9605 | 5641.99 | 52.7605 |
| Arabidopsis thaliana | transcriptomic | root_meristem | 552 | 310.382 | 13.7146 |
| Arabidopsis thaliana | transcriptomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | transcriptomic | flower | 5135 | 5452.02 | -4.2934 |
| Arabidopsis thaliana | transcriptomic | fruit | 313 | 188.905 | 9.02888 |
| Arabidopsis thaliana | transcriptomic | seed | 2782 | 1862.83 | 21.2967 |
| Arabidopsis thaliana | transcriptomic | whole | 24413 | 17284 | 54.2261 |
| Arabidopsis thaliana | transcriptomic | in_vitro | 440 | 636.282 | -7.78138 |
| Arabidopsis thaliana | transcriptomic | other | 459 | 269.176 | 11.57 |
| Arabidopsis thaliana | transcriptomic | mixed | 3236 | 5400.64 | -29.4553 |
| Arabidopsis thaliana | transcriptomic | unknown | 3064 | 2609.88 | 8.88911 |
| Arabidopsis thaliana | epigenomic | leaf | 5701 | 8150.67 | -27.1338 |
| Arabidopsis thaliana | epigenomic | shoot_meristem | 71 | 123.285 | -4.70892 |
| Arabidopsis thaliana | epigenomic | stem | 744 | 1019.54 | -8.62947 |
| Arabidopsis thaliana | epigenomic | root | 445 | 1270.57 | -23.1609 |
| Arabidopsis thaliana | epigenomic | root_meristem | 28 | 69.8976 | -5.01139 |
| Arabidopsis thaliana | epigenomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | epigenomic | flower | 1670 | 1227.79 | 12.6203 |
| Arabidopsis thaliana | epigenomic | fruit | 26 | 42.5411 | -2.53607 |
| Arabidopsis thaliana | epigenomic | seed | 520 | 419.506 | 4.90647 |
| Arabidopsis thaliana | epigenomic | whole | 5812 | 3892.33 | 30.7695 |
| Arabidopsis thaliana | epigenomic | in_vitro | 217 | 143.29 | 6.15768 |
| Arabidopsis thaliana | epigenomic | other | 43 | 60.6181 | -2.26286 |
| Arabidopsis thaliana | epigenomic | mixed | 2084 | 1216.22 | 24.8831 |
| Arabidopsis thaliana | epigenomic | unknown | 863 | 587.743 | 11.3539 |
| Arabidopsis thaliana | chromatin | leaf | 18 | 51.4336 | -4.66187 |
| Arabidopsis thaliana | chromatin | shoot_meristem | 0 | 0.777972 | -0.882027 |
| Arabidopsis thaliana | chromatin | stem | 29 | 6.43367 | 8.89675 |
| Arabidopsis thaliana | chromatin | root | 4 | 8.01776 | -1.41892 |
| Arabidopsis thaliana | chromatin | root_meristem | 0 | 0.441079 | -0.664138 |
| Arabidopsis thaliana | chromatin | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | chromatin | flower | 3 | 7.74778 | -1.7057 |
| Arabidopsis thaliana | chromatin | fruit | 0 | 0.26845 | -0.518122 |
| Arabidopsis thaliana | chromatin | seed | 6 | 2.64724 | 2.06066 |
| Arabidopsis thaliana | chromatin | whole | 30 | 24.562 | 1.09725 |
| Arabidopsis thaliana | chromatin | in_vitro | 0 | 0.904212 | -0.950901 |
| Arabidopsis thaliana | chromatin | other | 0 | 0.382522 | -0.618484 |
| Arabidopsis thaliana | chromatin | mixed | 4 | 7.67478 | -1.32647 |
| Arabidopsis thaliana | chromatin | unknown | 21 | 3.70887 | 8.97848 |
| Arabidopsis thaliana | other | leaf | 87 | 38.9107 | 7.70929 |
| Arabidopsis thaliana | other | shoot_meristem | 0 | 0.588553 | -0.767172 |
| Arabidopsis thaliana | other | stem | 0 | 4.86721 | -2.20618 |
| Arabidopsis thaliana | other | root | 0 | 6.06561 | -2.46285 |
| Arabidopsis thaliana | other | root_meristem | 0 | 0.333686 | -0.577656 |
| Arabidopsis thaliana | other | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | other | flower | 0 | 5.86137 | -2.42103 |
| Arabidopsis thaliana | other | fruit | 0 | 0.203088 | -0.450653 |
| Arabidopsis thaliana | other | seed | 0 | 2.00269 | -1.41516 |
| Arabidopsis thaliana | other | whole | 0 | 18.5817 | -4.31065 |
| Arabidopsis thaliana | other | in_vitro | 0 | 0.684056 | -0.827077 |
| Arabidopsis thaliana | other | other | 0 | 0.289386 | -0.537946 |
| Arabidopsis thaliana | other | mixed | 0 | 5.80614 | -2.40959 |
| Arabidopsis thaliana | other | unknown | 0 | 2.80584 | -1.67506 |
| Arabidopsis thaliana | mixed | leaf | 77 | 124.782 | -4.27752 |
| Arabidopsis thaliana | mixed | shoot_meristem | 0 | 1.88743 | -1.37384 |
| Arabidopsis thaliana | mixed | stem | 35 | 15.6086 | 4.90824 |
| Arabidopsis thaliana | mixed | root | 20 | 19.4518 | 0.124303 |
| Arabidopsis thaliana | mixed | root_meristem | 0 | 1.0701 | -1.03445 |
| Arabidopsis thaliana | mixed | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | mixed | flower | 29 | 18.7968 | 2.35339 |
| Arabidopsis thaliana | mixed | fruit | 0 | 0.651283 | -0.807021 |
| Arabidopsis thaliana | mixed | seed | 7 | 6.42242 | 0.227908 |
| Arabidopsis thaliana | mixed | whole | 61 | 59.5896 | 0.182705 |
| Arabidopsis thaliana | mixed | in_vitro | 7 | 2.1937 | 3.24506 |
| Arabidopsis thaliana | mixed | other | 0 | 0.928032 | -0.963344 |
| Arabidopsis thaliana | mixed | mixed | 24 | 18.6197 | 1.24687 |
| Arabidopsis thaliana | mixed | unknown | 19 | 8.99804 | 3.33435 |
| Arabidopsis thaliana | unknown | leaf | 40 | 91.6861 | -5.39786 |
| Arabidopsis thaliana | unknown | shoot_meristem | 0 | 1.38682 | -1.17763 |
| Arabidopsis thaliana | unknown | stem | 18 | 11.4687 | 1.92859 |
| Arabidopsis thaliana | unknown | root | 17 | 14.2925 | 0.716161 |
| Arabidopsis thaliana | unknown | root_meristem | 0 | 0.786272 | -0.88672 |
| Arabidopsis thaliana | unknown | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | unknown | flower | 38 | 13.8113 | 6.50873 |
| Arabidopsis thaliana | unknown | fruit | 0 | 0.478541 | -0.691767 |
| Arabidopsis thaliana | unknown | seed | 12 | 4.71899 | 3.35172 |
| Arabidopsis thaliana | unknown | whole | 55 | 43.7845 | 1.69496 |
| Arabidopsis thaliana | unknown | in_vitro | 0 | 1.61186 | -1.26959 |
| Arabidopsis thaliana | unknown | other | 0 | 0.681887 | -0.825765 |
| Arabidopsis thaliana | unknown | mixed | 8 | 13.6811 | -1.53594 |
| Arabidopsis thaliana | unknown | unknown | 17 | 6.61146 | 4.04023 |

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
