# GAMA  
**Genomic Availability & Metadata Analysis Tool**

<!-- badges: start -->
[![R-CMD-check](https://github.com/JLewis-dev/GAMA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JLewis-dev/GAMA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> &#9888;&#65038; **Development Status: Active**  
> GAMA is currently in early development. Interfaces, methods, and outputs may change.  
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
| Arabidopsis thaliana | 378 | 242822 | 243175 | 17.8104 | 24.8002 | 12.4015 | 55.0121 |
| Glycine max | 51 | 46476 | 58420 | 15.8777 | 21.4934 | 10.9754 | 48.3466 |
| Phaseolus vulgaris | 16 | 9620 | 9939 | 14.5109 | 18.3434 | 9.20432 | 42.0586 |
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
| Arabidopsis thaliana | 242822 | 88832 | 125720 | 27141 | 717 | 87 | 325 |
| Glycine max | 46476 | 25508 | 19906 | 1002 | 24 | 0 | 36 |
| Phaseolus vulgaris | 9620 | 6661 | 2756 | 181 | 4 | 0 | 18 |
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
| Arabidopsis thaliana | 5781 | transcriptomic | 1 | 6 | 11 | 18 | 5184 | 218.121 |
| Glycine max | 1298 | transcriptomic | 1 | 1 | 6 | 16 | 801 | 153.741 |
| Phaseolus vulgaris | 145 | transcriptomic | 1 | 5 | 12 | 21 | 200 | 53.2882 |
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
| Arabidopsis thaliana | genomic | leaf | 39090 | 22990.1 | 106.182 |
| Arabidopsis thaliana | genomic | shoot_meristem | 14 | 348.971 | -17.9313 |
| Arabidopsis thaliana | genomic | stem | 265 | 2877.39 | -48.7011 |
| Arabidopsis thaliana | genomic | root | 452 | 3589.32 | -52.3664 |
| Arabidopsis thaliana | genomic | root_meristem | 0 | 197.853 | -14.066 |
| Arabidopsis thaliana | genomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | genomic | flower | 3281 | 3463.79 | -3.10578 |
| Arabidopsis thaliana | genomic | fruit | 14 | 120.417 | -9.69768 |
| Arabidopsis thaliana | genomic | seed | 154 | 1183.36 | -29.9233 |
| Arabidopsis thaliana | genomic | whole | 1927 | 10948.8 | -86.2202 |
| Arabidopsis thaliana | genomic | in_vitro | 525 | 405.598 | 5.92875 |
| Arabidopsis thaliana | genomic | other | 1 | 171.586 | -13.0227 |
| Arabidopsis thaliana | genomic | mixed | 4736 | 3395.22 | 23.0103 |
| Arabidopsis thaliana | genomic | unknown | 893 | 1659.58 | -18.8173 |
| Arabidopsis thaliana | transcriptomic | leaf | 22382 | 35951.9 | -71.5674 |
| Arabidopsis thaliana | transcriptomic | shoot_meristem | 938 | 545.72 | 16.7924 |
| Arabidopsis thaliana | transcriptomic | stem | 7337 | 4499.65 | 42.2983 |
| Arabidopsis thaliana | transcriptomic | root | 9586 | 5612.96 | 53.0306 |
| Arabidopsis thaliana | transcriptomic | root_meristem | 552 | 309.401 | 13.792 |
| Arabidopsis thaliana | transcriptomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | transcriptomic | flower | 5133 | 5416.65 | -3.8541 |
| Arabidopsis thaliana | transcriptomic | fruit | 313 | 188.308 | 9.08666 |
| Arabidopsis thaliana | transcriptomic | seed | 2770 | 1850.54 | 21.3739 |
| Arabidopsis thaliana | transcriptomic | whole | 24236 | 17121.6 | 54.3707 |
| Arabidopsis thaliana | transcriptomic | in_vitro | 440 | 634.272 | -7.71388 |
| Arabidopsis thaliana | transcriptomic | other | 459 | 268.325 | 11.6402 |
| Arabidopsis thaliana | transcriptomic | mixed | 3106 | 5309.43 | -30.2396 |
| Arabidopsis thaliana | transcriptomic | unknown | 3052 | 2595.24 | 8.96611 |
| Arabidopsis thaliana | epigenomic | leaf | 5701 | 8149.88 | -27.1264 |
| Arabidopsis thaliana | epigenomic | shoot_meristem | 71 | 123.708 | -4.73893 |
| Arabidopsis thaliana | epigenomic | stem | 751 | 1020.02 | -8.42326 |
| Arabidopsis thaliana | epigenomic | root | 443 | 1272.39 | -23.2515 |
| Arabidopsis thaliana | epigenomic | root_meristem | 28 | 70.1377 | -5.03147 |
| Arabidopsis thaliana | epigenomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | epigenomic | flower | 1670 | 1227.89 | 12.6167 |
| Arabidopsis thaliana | epigenomic | fruit | 26 | 42.6873 | -2.55409 |
| Arabidopsis thaliana | epigenomic | seed | 520 | 419.496 | 4.90704 |
| Arabidopsis thaliana | epigenomic | whole | 5796 | 3881.28 | 30.734 |
| Arabidopsis thaliana | epigenomic | in_vitro | 217 | 143.782 | 6.10609 |
| Arabidopsis thaliana | epigenomic | other | 43 | 60.8263 | -2.28568 |
| Arabidopsis thaliana | epigenomic | mixed | 2075 | 1203.59 | 25.118 |
| Arabidopsis thaliana | epigenomic | unknown | 863 | 588.31 | 11.325 |
| Arabidopsis thaliana | chromatin | leaf | 18 | 51.4852 | -4.66672 |
| Arabidopsis thaliana | chromatin | shoot_meristem | 0 | 0.781502 | -0.884026 |
| Arabidopsis thaliana | chromatin | stem | 29 | 6.44376 | 8.88581 |
| Arabidopsis thaliana | chromatin | root | 4 | 8.03809 | -1.42429 |
| Arabidopsis thaliana | chromatin | root_meristem | 0 | 0.44308 | -0.665643 |
| Arabidopsis thaliana | chromatin | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | chromatin | flower | 3 | 7.75696 | -1.70798 |
| Arabidopsis thaliana | chromatin | fruit | 0 | 0.269668 | -0.519296 |
| Arabidopsis thaliana | chromatin | seed | 6 | 2.65008 | 2.05781 |
| Arabidopsis thaliana | chromatin | whole | 30 | 24.5192 | 1.10687 |
| Arabidopsis thaliana | chromatin | in_vitro | 0 | 0.908315 | -0.953056 |
| Arabidopsis thaliana | chromatin | other | 0 | 0.384258 | -0.619885 |
| Arabidopsis thaliana | chromatin | mixed | 4 | 7.60341 | -1.3068 |
| Arabidopsis thaliana | chromatin | unknown | 21 | 3.71653 | 8.96525 |
| Arabidopsis thaliana | other | leaf | 87 | 38.9497 | 7.69918 |
| Arabidopsis thaliana | other | shoot_meristem | 0 | 0.591223 | -0.768911 |
| Arabidopsis thaliana | other | stem | 0 | 4.87485 | -2.20791 |
| Arabidopsis thaliana | other | root | 0 | 6.08099 | -2.46597 |
| Arabidopsis thaliana | other | root_meristem | 0 | 0.3352 | -0.578965 |
| Arabidopsis thaliana | other | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | other | flower | 0 | 5.86831 | -2.42246 |
| Arabidopsis thaliana | other | fruit | 0 | 0.20401 | -0.451674 |
| Arabidopsis thaliana | other | seed | 0 | 2.00484 | -1.41592 |
| Arabidopsis thaliana | other | whole | 0 | 18.5493 | -4.30689 |
| Arabidopsis thaliana | other | in_vitro | 0 | 0.68716 | -0.828951 |
| Arabidopsis thaliana | other | other | 0 | 0.290699 | -0.539165 |
| Arabidopsis thaliana | other | mixed | 0 | 5.75215 | -2.39836 |
| Arabidopsis thaliana | other | unknown | 0 | 2.81163 | -1.67679 |
| Arabidopsis thaliana | mixed | leaf | 77 | 120.878 | -3.99094 |
| Arabidopsis thaliana | mixed | shoot_meristem | 0 | 1.83483 | -1.35456 |
| Arabidopsis thaliana | mixed | stem | 35 | 15.1288 | 5.10882 |
| Arabidopsis thaliana | mixed | root | 20 | 18.872 | 0.259648 |
| Arabidopsis thaliana | mixed | root_meristem | 0 | 1.04028 | -1.01994 |
| Arabidopsis thaliana | mixed | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | mixed | flower | 29 | 18.212 | 2.52791 |
| Arabidopsis thaliana | mixed | fruit | 0 | 0.633133 | -0.795697 |
| Arabidopsis thaliana | mixed | seed | 7 | 6.22193 | 0.311931 |
| Arabidopsis thaliana | mixed | whole | 52 | 57.5667 | -0.733691 |
| Arabidopsis thaliana | mixed | in_vitro | 7 | 2.13257 | 3.3331 |
| Arabidopsis thaliana | mixed | other | 0 | 0.90217 | -0.949826 |
| Arabidopsis thaliana | mixed | mixed | 24 | 17.8515 | 1.45523 |
| Arabidopsis thaliana | mixed | unknown | 19 | 8.72576 | 3.47815 |
| Arabidopsis thaliana | unknown | leaf | 40 | 91.7779 | -5.40475 |
| Arabidopsis thaliana | unknown | shoot_meristem | 0 | 1.39311 | -1.1803 |
| Arabidopsis thaliana | unknown | stem | 18 | 11.4867 | 1.92177 |
| Arabidopsis thaliana | unknown | root | 17 | 14.3288 | 0.705679 |
| Arabidopsis thaliana | unknown | root_meristem | 0 | 0.789839 | -0.888729 |
| Arabidopsis thaliana | unknown | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | unknown | flower | 38 | 13.8276 | 6.50048 |
| Arabidopsis thaliana | unknown | fruit | 0 | 0.480712 | -0.693334 |
| Arabidopsis thaliana | unknown | seed | 12 | 4.72405 | 3.34759 |
| Arabidopsis thaliana | unknown | whole | 55 | 43.7081 | 1.708 |
| Arabidopsis thaliana | unknown | in_vitro | 0 | 1.61917 | -1.27247 |
| Arabidopsis thaliana | unknown | other | 0 | 0.684981 | -0.827636 |
| Arabidopsis thaliana | unknown | mixed | 8 | 13.5539 | -1.50857 |
| Arabidopsis thaliana | unknown | unknown | 17 | 6.62512 | 4.03076 |

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

S = 2·ln(1 + SRA)  
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
- Runtime increases with species list size 
- Novel protocols may not be fully captured by the modality ontology 
- The anatomy ontology is broad but not exhaustive and will require refinement
- Results should be interpreted cautiously during early development 

---

## Licence

See the `LICENSE` files for details.

---
