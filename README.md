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
ASM <- extract_assembly_metadata(RESULTS, best = TRUE)
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
| Arabidopsis thaliana | 378 | 240113 | 239446 | 17.8104 | 24.7777 | 12.3861 | 54.9742 |
| Glycine max | 51 | 45881 | 57903 | 15.8777 | 21.4677 | 10.9665 | 48.3119 |
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
| Arabidopsis thaliana | 240113 | 87932 | 124293 | 26759 | 717 | 87 | 325 |
| Glycine max | 45881 | 25455 | 19404 | 962 | 24 | 0 | 36 |
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
| Arabidopsis thaliana | 5712 | transcriptomic | 1 | 6 | 11 | 18 | 5184 | 213.417 |
| Glycine max | 1270 | transcriptomic | 1 | 1 | 6 | 16 | 801 | 146.983 |
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
| Arabidopsis thaliana | genomic | leaf | 36794 | 21542.2 | 103.914 |
| Arabidopsis thaliana | genomic | shoot_meristem | 14 | 338.38 | -17.634 |
| Arabidopsis thaliana | genomic | stem | 253 | 2776.45 | -47.8906 |
| Arabidopsis thaliana | genomic | root | 440 | 3454.21 | -51.2861 |
| Arabidopsis thaliana | genomic | root_meristem | 0 | 186.877 | -13.6703 |
| Arabidopsis thaliana | genomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | genomic | flower | 3281 | 3380.13 | -1.70506 |
| Arabidopsis thaliana | genomic | fruit | 14 | 114.796 | -9.4076 |
| Arabidopsis thaliana | genomic | seed | 154 | 1150.29 | -29.3753 |
| Arabidopsis thaliana | genomic | whole | 1855 | 10543.8 | -84.618 |
| Arabidopsis thaliana | genomic | in_vitro | 525 | 383.431 | 7.22979 |
| Arabidopsis thaliana | genomic | other | 1 | 167.188 | -12.8528 |
| Arabidopsis thaliana | genomic | mixed | 4736 | 3301.38 | 24.9684 |
| Arabidopsis thaliana | genomic | unknown | 893 | 1620.82 | -18.0783 |
| Arabidopsis thaliana | transcriptomic | leaf | 21894 | 34824.4 | -69.2897 |
| Arabidopsis thaliana | transcriptomic | shoot_meristem | 929 | 547.013 | 16.3324 |
| Arabidopsis thaliana | transcriptomic | stem | 7244 | 4488.31 | 41.1328 |
| Arabidopsis thaliana | transcriptomic | root | 9447 | 5583.96 | 51.6961 |
| Arabidopsis thaliana | transcriptomic | root_meristem | 532 | 302.098 | 13.2272 |
| Arabidopsis thaliana | transcriptomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | transcriptomic | flower | 5108 | 5464.2 | -4.8187 |
| Arabidopsis thaliana | transcriptomic | fruit | 304 | 185.575 | 8.69332 |
| Arabidopsis thaliana | transcriptomic | seed | 2748 | 1859.52 | 20.6037 |
| Arabidopsis thaliana | transcriptomic | whole | 23961 | 17044.8 | 52.975 |
| Arabidopsis thaliana | transcriptomic | in_vitro | 421 | 619.841 | -7.98665 |
| Arabidopsis thaliana | transcriptomic | other | 457 | 270.27 | 11.3584 |
| Arabidopsis thaliana | transcriptomic | mixed | 3046 | 5336.89 | -31.3588 |
| Arabidopsis thaliana | transcriptomic | unknown | 3056 | 2620.16 | 8.51455 |
| Arabidopsis thaliana | epigenomic | leaf | 5645 | 7890.45 | -25.2785 |
| Arabidopsis thaliana | epigenomic | shoot_meristem | 71 | 123.941 | -4.7554 |
| Arabidopsis thaliana | epigenomic | stem | 742 | 1016.96 | -8.62206 |
| Arabidopsis thaliana | epigenomic | root | 423 | 1265.2 | -23.6776 |
| Arabidopsis thaliana | epigenomic | root_meristem | 28 | 68.4489 | -4.88904 |
| Arabidopsis thaliana | epigenomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | epigenomic | flower | 1670 | 1238.07 | 12.2756 |
| Arabidopsis thaliana | epigenomic | fruit | 26 | 42.0472 | -2.47474 |
| Arabidopsis thaliana | epigenomic | seed | 520 | 421.327 | 4.80714 |
| Arabidopsis thaliana | epigenomic | whole | 5643 | 3861.98 | 28.6591 |
| Arabidopsis thaliana | epigenomic | in_vitro | 196 | 140.442 | 4.68807 |
| Arabidopsis thaliana | epigenomic | other | 43 | 61.2373 | -2.33052 |
| Arabidopsis thaliana | epigenomic | mixed | 2075 | 1209.22 | 24.8973 |
| Arabidopsis thaliana | epigenomic | unknown | 851 | 593.672 | 10.5612 |
| Arabidopsis thaliana | chromatin | leaf | 18 | 50.5995 | -4.58288 |
| Arabidopsis thaliana | chromatin | shoot_meristem | 0 | 0.794806 | -0.891519 |
| Arabidopsis thaliana | chromatin | stem | 29 | 6.52149 | 8.80226 |
| Arabidopsis thaliana | chromatin | root | 4 | 8.11345 | -1.44412 |
| Arabidopsis thaliana | chromatin | root_meristem | 0 | 0.438946 | -0.66253 |
| Arabidopsis thaliana | chromatin | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | chromatin | flower | 3 | 7.93944 | -1.753 |
| Arabidopsis thaliana | chromatin | fruit | 0 | 0.269638 | -0.519267 |
| Arabidopsis thaliana | chromatin | seed | 6 | 2.70187 | 2.00648 |
| Arabidopsis thaliana | chromatin | whole | 30 | 24.766 | 1.05174 |
| Arabidopsis thaliana | chromatin | in_vitro | 0 | 0.900624 | -0.949012 |
| Arabidopsis thaliana | chromatin | other | 0 | 0.3927 | -0.626658 |
| Arabidopsis thaliana | chromatin | mixed | 4 | 7.75446 | -1.34825 |
| Arabidopsis thaliana | chromatin | unknown | 21 | 3.80707 | 8.81159 |
| Arabidopsis thaliana | other | leaf | 87 | 38.2796 | 7.87456 |
| Arabidopsis thaliana | other | shoot_meristem | 0 | 0.601288 | -0.775428 |
| Arabidopsis thaliana | other | stem | 0 | 4.93365 | -2.22118 |
| Arabidopsis thaliana | other | root | 0 | 6.138 | -2.4775 |
| Arabidopsis thaliana | other | root_meristem | 0 | 0.332072 | -0.576257 |
| Arabidopsis thaliana | other | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | other | flower | 0 | 6.00636 | -2.45079 |
| Arabidopsis thaliana | other | fruit | 0 | 0.203987 | -0.45165 |
| Arabidopsis thaliana | other | seed | 0 | 2.04402 | -1.42969 |
| Arabidopsis thaliana | other | whole | 0 | 18.736 | -4.32851 |
| Arabidopsis thaliana | other | in_vitro | 0 | 0.681341 | -0.825434 |
| Arabidopsis thaliana | other | other | 0 | 0.297086 | -0.545056 |
| Arabidopsis thaliana | other | mixed | 0 | 5.86641 | -2.42207 |
| Arabidopsis thaliana | other | unknown | 0 | 2.88013 | -1.6971 |
| Arabidopsis thaliana | mixed | leaf | 76 | 117.919 | -3.86027 |
| Arabidopsis thaliana | mixed | shoot_meristem | 0 | 1.85224 | -1.36097 |
| Arabidopsis thaliana | mixed | stem | 34 | 15.1979 | 4.82297 |
| Arabidopsis thaliana | mixed | root | 20 | 18.9079 | 0.251162 |
| Arabidopsis thaliana | mixed | root_meristem | 0 | 1.02294 | -1.0114 |
| Arabidopsis thaliana | mixed | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | mixed | flower | 29 | 18.5023 | 2.4405 |
| Arabidopsis thaliana | mixed | fruit | 0 | 0.628375 | -0.792701 |
| Arabidopsis thaliana | mixed | seed | 7 | 6.29653 | 0.280345 |
| Arabidopsis thaliana | mixed | whole | 52 | 57.7155 | -0.752328 |
| Arabidopsis thaliana | mixed | in_vitro | 7 | 2.09884 | 3.38305 |
| Arabidopsis thaliana | mixed | other | 0 | 0.915162 | -0.956641 |
| Arabidopsis thaliana | mixed | mixed | 24 | 18.0713 | 1.39466 |
| Arabidopsis thaliana | mixed | unknown | 19 | 8.87214 | 3.40019 |
| Arabidopsis thaliana | unknown | leaf | 40 | 90.1992 | -5.28561 |
| Arabidopsis thaliana | unknown | shoot_meristem | 0 | 1.41683 | -1.19031 |
| Arabidopsis thaliana | unknown | stem | 18 | 11.6253 | 1.86965 |
| Arabidopsis thaliana | unknown | root | 17 | 14.4631 | 0.667069 |
| Arabidopsis thaliana | unknown | root_meristem | 0 | 0.782469 | -0.884573 |
| Arabidopsis thaliana | unknown | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | unknown | flower | 38 | 14.1529 | 6.33888 |
| Arabidopsis thaliana | unknown | fruit | 0 | 0.48066 | -0.693296 |
| Arabidopsis thaliana | unknown | seed | 12 | 4.81638 | 3.27328 |
| Arabidopsis thaliana | unknown | whole | 55 | 44.148 | 1.63325 |
| Arabidopsis thaliana | unknown | in_vitro | 0 | 1.60546 | -1.26707 |
| Arabidopsis thaliana | unknown | other | 0 | 0.700031 | -0.836678 |
| Arabidopsis thaliana | unknown | mixed | 8 | 13.8232 | -1.56623 |
| Arabidopsis thaliana | unknown | unknown | 17 | 6.78652 | 3.92058 |

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
