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

See the **[GAMA user guide](docs/GAMA_user_guide.pdf)** for a comprehensive overview of functions and methods.

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

Data richness provides a weighted overview of sequence availability across species, combining genome assembly quality with transformed SRA and BioSample accession counts.

<details>
<summary>Data richness table</summary>

| species | Assembly | SRA | BioSample | A | S | B | score |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis thaliana | 378 | 240585 | 241449 | 17.8104 | 24.7817 | 12.3944 | 54.9864 |
| Glycine max | 51 | 45994 | 57969 | 15.8777 | 21.4726 | 10.9677 | 48.318 |
| Phaseolus vulgaris | 16 | 9614 | 9923 | 14.5109 | 18.3422 | 9.20271 | 42.0557 |
| Vigna radiata | 9 | 5539 | 5623 | 13.9703 | 17.2395 | 8.6348 | 39.8446 |

</details>

### Sequencing modality

<p align='center'>
  <img src='man/figures/Modality_plot.png' alt='GAMA modality plot' width='800'>
</p>

Sequencing modality summarises the proportional composition of experimental classes. Similar outputs are available for Assembly level and BioSample anatomy composition.

<details>
<summary>Sequencing modality table</summary>

| species | SRA | genomic | transcriptomic | epigenomic | chromatin | other | unknown |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis thaliana | 240585 | 87989 | 124647 | 26820 | 717 | 87 | 325 |
| Glycine max | 45994 | 25504 | 19468 | 962 | 24 | 0 | 36 |
| Phaseolus vulgaris | 9614 | 6661 | 2750 | 181 | 4 | 0 | 18 |
| Vigna radiata | 5539 | 4761 | 733 | 18 | 27 | 0 | 0 |

</details>

### Replication skew

<p align='center'>
  <img src='man/figures/Skew_plot.png' alt='GAMA skew plot' width='800'>
</p>

Replication skew reveals how broadly experiments are distributed across BioProjects or BioSamples.

<details>
<summary>Replication skew table</summary>

| species | BioProject | class | min | q25 | med | q75 | max | eff |
| :--- | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis thaliana | 5728 | transcriptomic | 1 | 6 | 11 | 18 | 5184 | 214.565 |
| Glycine max | 1274 | transcriptomic | 1 | 1 | 6 | 16 | 801 | 147.873 |
| Phaseolus vulgaris | 144 | transcriptomic | 1 | 4.75 | 12 | 21 | 200 | 53.0541 |
| Vigna radiata | 64 | transcriptomic | 1 | 4 | 7 | 12.75 | 65 | 29.0066 |

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
| Arabidopsis thaliana | genomic | leaf | 38385 | 22595.7 | 105.039 |
| Arabidopsis thaliana | genomic | shoot_meristem | 14 | 345.066 | -17.8223 |
| Arabidopsis thaliana | genomic | stem | 253 | 2839.48 | -48.5389 |
| Arabidopsis thaliana | genomic | root | 440 | 3528.59 | -51.9948 |
| Arabidopsis thaliana | genomic | root_meristem | 0 | 190.569 | -13.8047 |
| Arabidopsis thaliana | genomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | genomic | flower | 3281 | 3453.38 | -2.93342 |
| Arabidopsis thaliana | genomic | fruit | 14 | 117.064 | -9.52566 |
| Arabidopsis thaliana | genomic | seed | 154 | 1173.02 | -29.753 |
| Arabidopsis thaliana | genomic | whole | 1912 | 10780.8 | -85.4158 |
| Arabidopsis thaliana | genomic | in_vitro | 525 | 392.368 | 6.69578 |
| Arabidopsis thaliana | genomic | other | 1 | 170.491 | -12.9806 |
| Arabidopsis thaliana | genomic | mixed | 4736 | 3368.65 | 23.5588 |
| Arabidopsis thaliana | genomic | unknown | 893 | 1652.85 | -18.69 |
| Arabidopsis thaliana | transcriptomic | leaf | 22141 | 35488 | -70.8503 |
| Arabidopsis thaliana | transcriptomic | shoot_meristem | 929 | 541.948 | 16.6261 |
| Arabidopsis thaliana | transcriptomic | stem | 7259 | 4459.58 | 41.92 |
| Arabidopsis thaliana | transcriptomic | root | 9465 | 5541.87 | 52.6993 |
| Arabidopsis thaliana | transcriptomic | root_meristem | 532 | 299.301 | 13.4506 |
| Arabidopsis thaliana | transcriptomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | transcriptomic | flower | 5127 | 5423.75 | -4.02945 |
| Arabidopsis thaliana | transcriptomic | fruit | 304 | 183.856 | 8.8606 |
| Arabidopsis thaliana | transcriptomic | seed | 2748 | 1842.3 | 21.101 |
| Arabidopsis thaliana | transcriptomic | whole | 23988 | 16931.9 | 54.2269 |
| Arabidopsis thaliana | transcriptomic | in_vitro | 425 | 616.238 | -7.70373 |
| Arabidopsis thaliana | transcriptomic | other | 457 | 267.767 | 11.5643 |
| Arabidopsis thaliana | transcriptomic | mixed | 3052 | 5290.67 | -30.7776 |
| Arabidopsis thaliana | transcriptomic | unknown | 3056 | 2595.9 | 9.03047 |
| Arabidopsis thaliana | epigenomic | leaf | 5652 | 8013.96 | -26.3845 |
| Arabidopsis thaliana | epigenomic | shoot_meristem | 71 | 122.384 | -4.64476 |
| Arabidopsis thaliana | epigenomic | stem | 751 | 1007.07 | -8.06918 |
| Arabidopsis thaliana | epigenomic | root | 423 | 1251.48 | -23.419 |
| Arabidopsis thaliana | epigenomic | root_meristem | 28 | 67.5886 | -4.81541 |
| Arabidopsis thaliana | epigenomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | epigenomic | flower | 1670 | 1224.8 | 12.721 |
| Arabidopsis thaliana | epigenomic | fruit | 26 | 41.5187 | -2.40843 |
| Arabidopsis thaliana | epigenomic | seed | 520 | 416.032 | 5.09725 |
| Arabidopsis thaliana | epigenomic | whole | 5643 | 3823.58 | 29.4237 |
| Arabidopsis thaliana | epigenomic | in_vitro | 196 | 139.16 | 4.81832 |
| Arabidopsis thaliana | epigenomic | other | 43 | 60.4677 | -2.24633 |
| Arabidopsis thaliana | epigenomic | mixed | 2075 | 1194.75 | 25.4664 |
| Arabidopsis thaliana | epigenomic | unknown | 851 | 586.21 | 10.9364 |
| Arabidopsis thaliana | chromatin | leaf | 18 | 51.3458 | -4.65359 |
| Arabidopsis thaliana | chromatin | shoot_meristem | 0 | 0.784117 | -0.885504 |
| Arabidopsis thaliana | chromatin | stem | 29 | 6.45234 | 8.87652 |
| Arabidopsis thaliana | chromatin | root | 4 | 8.01826 | -1.41905 |
| Arabidopsis thaliana | chromatin | root_meristem | 0 | 0.433043 | -0.65806 |
| Arabidopsis thaliana | chromatin | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | chromatin | flower | 3 | 7.84736 | -1.73039 |
| Arabidopsis thaliana | chromatin | fruit | 0 | 0.266012 | -0.515764 |
| Arabidopsis thaliana | chromatin | seed | 6 | 2.66553 | 2.04237 |
| Arabidopsis thaliana | chromatin | whole | 30 | 24.4979 | 1.11165 |
| Arabidopsis thaliana | chromatin | in_vitro | 0 | 0.891605 | -0.944248 |
| Arabidopsis thaliana | chromatin | other | 0 | 0.387419 | -0.62243 |
| Arabidopsis thaliana | chromatin | mixed | 4 | 7.65481 | -1.32098 |
| Arabidopsis thaliana | chromatin | unknown | 21 | 3.75588 | 8.89786 |
| Arabidopsis thaliana | other | leaf | 87 | 38.8442 | 7.72655 |
| Arabidopsis thaliana | other | shoot_meristem | 0 | 0.593202 | -0.770196 |
| Arabidopsis thaliana | other | stem | 0 | 4.88134 | -2.20937 |
| Arabidopsis thaliana | other | root | 0 | 6.06599 | -2.46292 |
| Arabidopsis thaliana | other | root_meristem | 0 | 0.327606 | -0.572369 |
| Arabidopsis thaliana | other | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | other | flower | 0 | 5.9367 | -2.43653 |
| Arabidopsis thaliana | other | fruit | 0 | 0.201244 | -0.448602 |
| Arabidopsis thaliana | other | seed | 0 | 2.01653 | -1.42005 |
| Arabidopsis thaliana | other | whole | 0 | 18.5332 | -4.30502 |
| Arabidopsis thaliana | other | in_vitro | 0 | 0.674518 | -0.821291 |
| Arabidopsis thaliana | other | other | 0 | 0.293091 | -0.541379 |
| Arabidopsis thaliana | other | mixed | 0 | 5.79103 | -2.40646 |
| Arabidopsis thaliana | other | unknown | 0 | 2.8414 | -1.68565 |
| Arabidopsis thaliana | mixed | leaf | 76 | 119.658 | -3.9911 |
| Arabidopsis thaliana | mixed | shoot_meristem | 0 | 1.82733 | -1.35179 |
| Arabidopsis thaliana | mixed | stem | 34 | 15.0368 | 4.8903 |
| Arabidopsis thaliana | mixed | root | 20 | 18.686 | 0.303969 |
| Arabidopsis thaliana | mixed | root_meristem | 0 | 1.00918 | -1.00458 |
| Arabidopsis thaliana | mixed | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | mixed | flower | 29 | 18.2878 | 2.50496 |
| Arabidopsis thaliana | mixed | fruit | 0 | 0.619924 | -0.787353 |
| Arabidopsis thaliana | mixed | seed | 7 | 6.21185 | 0.316225 |
| Arabidopsis thaliana | mixed | whole | 52 | 57.0907 | -0.673741 |
| Arabidopsis thaliana | mixed | in_vitro | 7 | 2.07783 | 3.4147 |
| Arabidopsis thaliana | mixed | other | 0 | 0.902854 | -0.950187 |
| Arabidopsis thaliana | mixed | mixed | 24 | 17.839 | 1.45869 |
| Arabidopsis thaliana | mixed | unknown | 19 | 8.75282 | 3.46362 |
| Arabidopsis thaliana | unknown | leaf | 40 | 91.5294 | -5.3861 |
| Arabidopsis thaliana | unknown | shoot_meristem | 0 | 1.39777 | -1.18228 |
| Arabidopsis thaliana | unknown | stem | 18 | 11.502 | 1.91599 |
| Arabidopsis thaliana | unknown | root | 17 | 14.2934 | 0.715903 |
| Arabidopsis thaliana | unknown | root_meristem | 0 | 0.771946 | -0.878605 |
| Arabidopsis thaliana | unknown | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | unknown | flower | 38 | 13.9888 | 6.41985 |
| Arabidopsis thaliana | unknown | fruit | 0 | 0.474196 | -0.688619 |
| Arabidopsis thaliana | unknown | seed | 12 | 4.75161 | 3.32523 |
| Arabidopsis thaliana | unknown | whole | 55 | 43.6701 | 1.71448 |
| Arabidopsis thaliana | unknown | in_vitro | 0 | 1.58938 | -1.26071 |
| Arabidopsis thaliana | unknown | other | 0 | 0.690616 | -0.831033 |
| Arabidopsis thaliana | unknown | mixed | 8 | 13.6455 | -1.5283 |
| Arabidopsis thaliana | unknown | unknown | 17 | 6.69526 | 3.98248 |

</details>

---

## Transparency

GAMA outputs retain query provenance, including tool version, timestamp, database sources, search terms, and synonym handling. Full methodological details can be found in the **[GAMA user guide](docs/GAMA_user_guide.pdf)**.

### Data richness

The data richness score is defined as:

Score = A + S + B

Where A, S, and B are the transformed contributions of Assembly, SRA, and BioSample accession counts.

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

SRA experiments are classified into sequencing modality classes and subclasses using a curated ontology.

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
