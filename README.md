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
| Arabidopsis thaliana | 381 | 252139 | 252453 | 17.82 | 24.8755 | 12.439 | 55.1345 |
| Glycine max | 51 | 47134 | 59194 | 15.8777 | 21.5215 | 10.9886 | 48.3879 |
| Phaseolus vulgaris | 17 | 9856 | 10171 | 14.5951 | 18.3919 | 9.22739 | 42.2144 |
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
| Arabidopsis thaliana | 252139 | 89924 | 133136 | 27918 | 717 | 87 | 357 |
| Glycine max | 47134 | 25814 | 20251 | 1009 | 24 | 0 | 36 |
| Phaseolus vulgaris | 9856 | 6839 | 2814 | 181 | 4 | 0 | 18 |
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
| Arabidopsis thaliana | 5907 | transcriptomic | 1 | 6 | 11 | 18 | 5184 | 234.95 |
| Glycine max | 1316 | transcriptomic | 1 | 1 | 6 | 16 | 801 | 158.41 |
| Phaseolus vulgaris | 149 | transcriptomic | 1 | 5 | 12 | 21 | 200 | 55.3167 |
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
| Arabidopsis thaliana | genomic | leaf | 39492 | 24012.4 | 99.8945 |
| Arabidopsis thaliana | genomic | shoot_meristem | 14 | 343.547 | -17.7797 |
| Arabidopsis thaliana | genomic | stem | 265 | 2854.74 | -48.4699 |
| Arabidopsis thaliana | genomic | root | 453 | 3505.62 | -51.5573 |
| Arabidopsis thaliana | genomic | root_meristem | 0 | 193.143 | -13.8976 |
| Arabidopsis thaliana | genomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | genomic | flower | 3313 | 3307.9 | 0.0885916 |
| Arabidopsis thaliana | genomic | fruit | 15 | 115.494 | -9.35107 |
| Arabidopsis thaliana | genomic | seed | 155 | 1254.78 | -31.0472 |
| Arabidopsis thaliana | genomic | whole | 1930 | 10816 | -85.4424 |
| Arabidopsis thaliana | genomic | in_vitro | 525 | 387.918 | 6.96003 |
| Arabidopsis thaliana | genomic | other | 1 | 169 | -12.9231 |
| Arabidopsis thaliana | genomic | mixed | 4736 | 3304.97 | 24.8923 |
| Arabidopsis thaliana | genomic | unknown | 897 | 1530.46 | -16.1924 |
| Arabidopsis thaliana | transcriptomic | leaf | 27687 | 40524.8 | -63.7719 |
| Arabidopsis thaliana | transcriptomic | shoot_meristem | 968 | 579.79 | 16.1224 |
| Arabidopsis thaliana | transcriptomic | stem | 7573 | 4817.82 | 39.694 |
| Arabidopsis thaliana | transcriptomic | root | 9805 | 5916.28 | 50.557 |
| Arabidopsis thaliana | transcriptomic | root_meristem | 564 | 325.96 | 13.1846 |
| Arabidopsis thaliana | transcriptomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | transcriptomic | flower | 5106 | 5582.62 | -6.37896 |
| Arabidopsis thaliana | transcriptomic | fruit | 313 | 194.915 | 8.45806 |
| Arabidopsis thaliana | transcriptomic | seed | 3134 | 2117.64 | 22.0863 |
| Arabidopsis thaliana | transcriptomic | whole | 25172 | 18253.8 | 51.2058 |
| Arabidopsis thaliana | transcriptomic | in_vitro | 440 | 654.673 | -8.39007 |
| Arabidopsis thaliana | transcriptomic | other | 474 | 285.215 | 11.1784 |
| Arabidopsis thaliana | transcriptomic | mixed | 3268 | 5577.66 | -30.9259 |
| Arabidopsis thaliana | transcriptomic | unknown | 2910 | 2582.9 | 6.4361 |
| Arabidopsis thaliana | epigenomic | leaf | 6162 | 8729.04 | -27.4758 |
| Arabidopsis thaliana | epigenomic | shoot_meristem | 71 | 124.887 | -4.82198 |
| Arabidopsis thaliana | epigenomic | stem | 830 | 1037.76 | -6.44932 |
| Arabidopsis thaliana | epigenomic | root | 445 | 1274.37 | -23.2327 |
| Arabidopsis thaliana | epigenomic | root_meristem | 28 | 70.2119 | -5.03767 |
| Arabidopsis thaliana | epigenomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | epigenomic | flower | 1650 | 1202.5 | 12.9049 |
| Arabidopsis thaliana | epigenomic | fruit | 26 | 41.9848 | -2.46696 |
| Arabidopsis thaliana | epigenomic | seed | 532 | 456.14 | 3.55192 |
| Arabidopsis thaliana | epigenomic | whole | 5886 | 3931.87 | 31.1641 |
| Arabidopsis thaliana | epigenomic | in_vitro | 217 | 141.017 | 6.39856 |
| Arabidopsis thaliana | epigenomic | other | 43 | 61.4354 | -2.35203 |
| Arabidopsis thaliana | epigenomic | mixed | 2090 | 1201.43 | 25.6356 |
| Arabidopsis thaliana | epigenomic | unknown | 849 | 556.358 | 12.4068 |
| Arabidopsis thaliana | chromatin | leaf | 18 | 53.3135 | -4.8364 |
| Arabidopsis thaliana | chromatin | shoot_meristem | 0 | 0.76276 | -0.873361 |
| Arabidopsis thaliana | chromatin | stem | 29 | 6.33822 | 9.0014 |
| Arabidopsis thaliana | chromatin | root | 4 | 7.78334 | -1.3561 |
| Arabidopsis thaliana | chromatin | root_meristem | 0 | 0.428826 | -0.654848 |
| Arabidopsis thaliana | chromatin | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | chromatin | flower | 3 | 7.34437 | -1.60306 |
| Arabidopsis thaliana | chromatin | fruit | 0 | 0.256426 | -0.506386 |
| Arabidopsis thaliana | chromatin | seed | 6 | 2.78592 | 1.92563 |
| Arabidopsis thaliana | chromatin | whole | 46 | 24.0143 | 4.48649 |
| Arabidopsis thaliana | chromatin | in_vitro | 0 | 0.861274 | -0.928048 |
| Arabidopsis thaliana | chromatin | other | 0 | 0.375223 | -0.612554 |
| Arabidopsis thaliana | chromatin | mixed | 4 | 7.33785 | -1.2322 |
| Arabidopsis thaliana | chromatin | unknown | 5 | 3.39801 | 0.869054 |
| Arabidopsis thaliana | other | leaf | 87 | 40.3328 | 7.34822 |
| Arabidopsis thaliana | other | shoot_meristem | 0 | 0.577044 | -0.759634 |
| Arabidopsis thaliana | other | stem | 0 | 4.795 | -2.18975 |
| Arabidopsis thaliana | other | root | 0 | 5.88826 | -2.42657 |
| Arabidopsis thaliana | other | root_meristem | 0 | 0.324416 | -0.569576 |
| Arabidopsis thaliana | other | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | other | flower | 0 | 5.55618 | -2.35715 |
| Arabidopsis thaliana | other | fruit | 0 | 0.193992 | -0.440445 |
| Arabidopsis thaliana | other | seed | 0 | 2.10761 | -1.45176 |
| Arabidopsis thaliana | other | whole | 0 | 18.1673 | -4.26231 |
| Arabidopsis thaliana | other | in_vitro | 0 | 0.651573 | -0.8072 |
| Arabidopsis thaliana | other | other | 0 | 0.283864 | -0.532789 |
| Arabidopsis thaliana | other | mixed | 0 | 5.55124 | -2.35611 |
| Arabidopsis thaliana | other | unknown | 0 | 2.57067 | -1.60333 |
| Arabidopsis thaliana | mixed | leaf | 77 | 130.27 | -4.66727 |
| Arabidopsis thaliana | mixed | shoot_meristem | 0 | 1.86379 | -1.36521 |
| Arabidopsis thaliana | mixed | stem | 35 | 15.4873 | 4.95826 |
| Arabidopsis thaliana | mixed | root | 20 | 19.0184 | 0.225081 |
| Arabidopsis thaliana | mixed | root_meristem | 0 | 1.04783 | -1.02363 |
| Arabidopsis thaliana | mixed | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | mixed | flower | 29 | 17.9458 | 2.60943 |
| Arabidopsis thaliana | mixed | fruit | 0 | 0.626572 | -0.791563 |
| Arabidopsis thaliana | mixed | seed | 7 | 6.80734 | 0.0738431 |
| Arabidopsis thaliana | mixed | whole | 63 | 58.6783 | 0.564174 |
| Arabidopsis thaliana | mixed | in_vitro | 7 | 2.1045 | 3.3746 |
| Arabidopsis thaliana | mixed | other | 0 | 0.916849 | -0.957522 |
| Arabidopsis thaliana | mixed | mixed | 24 | 17.9299 | 1.43354 |
| Arabidopsis thaliana | mixed | unknown | 19 | 8.30297 | 3.71233 |
| Arabidopsis thaliana | unknown | leaf | 77 | 109.872 | -3.13606 |
| Arabidopsis thaliana | unknown | shoot_meristem | 0 | 1.57195 | -1.25377 |
| Arabidopsis thaliana | unknown | stem | 18 | 13.0623 | 1.36622 |
| Arabidopsis thaliana | unknown | root | 18 | 16.0404 | 0.489271 |
| Arabidopsis thaliana | unknown | root_meristem | 0 | 0.883755 | -0.940082 |
| Arabidopsis thaliana | unknown | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | unknown | flower | 38 | 15.1358 | 5.87697 |
| Arabidopsis thaliana | unknown | fruit | 0 | 0.528461 | -0.726953 |
| Arabidopsis thaliana | unknown | seed | 12 | 5.74142 | 2.61196 |
| Arabidopsis thaliana | unknown | whole | 55 | 49.4903 | 0.783198 |
| Arabidopsis thaliana | unknown | in_vitro | 0 | 1.77497 | -1.33228 |
| Arabidopsis thaliana | unknown | other | 0 | 0.773285 | -0.879366 |
| Arabidopsis thaliana | unknown | mixed | 8 | 15.1224 | -1.83153 |
| Arabidopsis thaliana | unknown | unknown | 11 | 7.00286 | 1.51047 |

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
