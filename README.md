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
| Arabidopsis thaliana | 378 | 241447 | 241916 | 17.8104 | 24.7888 | 12.3963 | 54.9955 |
| Glycine max | 51 | 46348 | 58055 | 15.8777 | 21.4879 | 10.9692 | 48.3348 |
| Phaseolus vulgaris | 16 | 9614 | 9933 | 14.5109 | 18.3422 | 9.20372 | 42.0567 |
| Vigna radiata | 9 | 5539 | 5623 | 13.9703 | 17.2395 | 8.6348 | 39.8446 |

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
| Arabidopsis thaliana | 241447 | 88052 | 125349 | 26917 | 717 | 87 | 325 |
| Glycine max | 46348 | 25508 | 19782 | 998 | 24 | 0 | 36 |
| Phaseolus vulgaris | 9614 | 6661 | 2750 | 181 | 4 | 0 | 18 |
| Vigna radiata | 5539 | 4761 | 733 | 18 | 27 | 0 | 0 |

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
| Arabidopsis thaliana | 5757 | transcriptomic | 1 | 6 | 11 | 18 | 5184 | 216.856 |
| Glycine max | 1290 | transcriptomic | 1 | 1 | 6 | 16 | 801 | 151.978 |
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
| Arabidopsis thaliana | genomic | leaf | 38393 | 22589.8 | 105.145 |
| Arabidopsis thaliana | genomic | shoot_meristem | 14 | 343.821 | -17.7874 |
| Arabidopsis thaliana | genomic | stem | 265 | 2844.15 | -48.3616 |
| Arabidopsis thaliana | genomic | root | 452 | 3554.51 | -52.0383 |
| Arabidopsis thaliana | genomic | root_meristem | 0 | 189.881 | -13.7797 |
| Arabidopsis thaliana | genomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | genomic | flower | 3281 | 3442.95 | -2.76008 |
| Arabidopsis thaliana | genomic | fruit | 14 | 119.693 | -9.66077 |
| Arabidopsis thaliana | genomic | seed | 154 | 1176.25 | -29.8062 |
| Arabidopsis thaliana | genomic | whole | 1915 | 10785.2 | -85.4124 |
| Arabidopsis thaliana | genomic | in_vitro | 525 | 396.038 | 6.48028 |
| Arabidopsis thaliana | genomic | other | 1 | 169.876 | -12.9569 |
| Arabidopsis thaliana | genomic | mixed | 4736 | 3374.8 | 23.4314 |
| Arabidopsis thaliana | genomic | unknown | 893 | 1656.03 | -18.7504 |
| Arabidopsis thaliana | transcriptomic | leaf | 22343 | 35695.9 | -70.6749 |
| Arabidopsis thaliana | transcriptomic | shoot_meristem | 929 | 543.298 | 16.5475 |
| Arabidopsis thaliana | transcriptomic | stem | 7291 | 4494.26 | 41.7179 |
| Arabidopsis thaliana | transcriptomic | root | 9543 | 5616.76 | 52.3883 |
| Arabidopsis thaliana | transcriptomic | root_meristem | 532 | 300.046 | 13.3908 |
| Arabidopsis thaliana | transcriptomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | transcriptomic | flower | 5133 | 5440.48 | -4.16868 |
| Arabidopsis thaliana | transcriptomic | fruit | 313 | 189.136 | 9.00651 |
| Arabidopsis thaliana | transcriptomic | seed | 2770 | 1858.68 | 21.1383 |
| Arabidopsis thaliana | transcriptomic | whole | 24097 | 17042.6 | 54.0369 |
| Arabidopsis thaliana | transcriptomic | in_vitro | 440 | 625.811 | -7.42761 |
| Arabidopsis thaliana | transcriptomic | other | 457 | 268.434 | 11.5092 |
| Arabidopsis thaliana | transcriptomic | mixed | 3106 | 5332.79 | -30.4931 |
| Arabidopsis thaliana | transcriptomic | unknown | 3071 | 2616.83 | 8.87829 |
| Arabidopsis thaliana | epigenomic | leaf | 5665 | 8035.3 | -26.4425 |
| Arabidopsis thaliana | epigenomic | shoot_meristem | 71 | 122.299 | -4.6387 |
| Arabidopsis thaliana | epigenomic | stem | 751 | 1011.68 | -8.19569 |
| Arabidopsis thaliana | epigenomic | root | 447 | 1264.36 | -22.9867 |
| Arabidopsis thaliana | epigenomic | root_meristem | 28 | 67.5418 | -4.81138 |
| Arabidopsis thaliana | epigenomic | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | epigenomic | flower | 1670 | 1224.68 | 12.7252 |
| Arabidopsis thaliana | epigenomic | fruit | 26 | 42.5755 | -2.5403 |
| Arabidopsis thaliana | epigenomic | seed | 520 | 418.397 | 4.96719 |
| Arabidopsis thaliana | epigenomic | whole | 5659 | 3836.37 | 29.4264 |
| Arabidopsis thaliana | epigenomic | in_vitro | 196 | 140.873 | 4.64463 |
| Arabidopsis thaliana | epigenomic | other | 43 | 60.4258 | -2.24172 |
| Arabidopsis thaliana | epigenomic | mixed | 2075 | 1200.43 | 25.2419 |
| Arabidopsis thaliana | epigenomic | unknown | 863 | 589.061 | 11.2869 |
| Arabidopsis thaliana | chromatin | leaf | 18 | 51.2968 | -4.64897 |
| Arabidopsis thaliana | chromatin | shoot_meristem | 0 | 0.780747 | -0.883599 |
| Arabidopsis thaliana | chromatin | stem | 29 | 6.45849 | 8.86988 |
| Arabidopsis thaliana | chromatin | root | 4 | 8.07157 | -1.43312 |
| Arabidopsis thaliana | chromatin | root_meristem | 0 | 0.431182 | -0.656644 |
| Arabidopsis thaliana | chromatin | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | chromatin | flower | 3 | 7.81825 | -1.72319 |
| Arabidopsis thaliana | chromatin | fruit | 0 | 0.271798 | -0.521343 |
| Arabidopsis thaliana | chromatin | seed | 6 | 2.67102 | 2.03692 |
| Arabidopsis thaliana | chromatin | whole | 30 | 24.4911 | 1.11316 |
| Arabidopsis thaliana | chromatin | in_vitro | 0 | 0.899322 | -0.948326 |
| Arabidopsis thaliana | chromatin | other | 0 | 0.385754 | -0.621091 |
| Arabidopsis thaliana | chromatin | mixed | 4 | 7.66348 | -1.32337 |
| Arabidopsis thaliana | chromatin | unknown | 21 | 3.76052 | 8.88997 |
| Arabidopsis thaliana | other | leaf | 87 | 38.8071 | 7.73619 |
| Arabidopsis thaliana | other | shoot_meristem | 0 | 0.590652 | -0.768539 |
| Arabidopsis thaliana | other | stem | 0 | 4.88598 | -2.21043 |
| Arabidopsis thaliana | other | root | 0 | 6.10632 | -2.4711 |
| Arabidopsis thaliana | other | root_meristem | 0 | 0.326198 | -0.571138 |
| Arabidopsis thaliana | other | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | other | flower | 0 | 5.91467 | -2.43201 |
| Arabidopsis thaliana | other | fruit | 0 | 0.205621 | -0.453455 |
| Arabidopsis thaliana | other | seed | 0 | 2.02068 | -1.42151 |
| Arabidopsis thaliana | other | whole | 0 | 18.5281 | -4.30442 |
| Arabidopsis thaliana | other | in_vitro | 0 | 0.680356 | -0.824837 |
| Arabidopsis thaliana | other | other | 0 | 0.291831 | -0.540214 |
| Arabidopsis thaliana | other | mixed | 0 | 5.79759 | -2.40782 |
| Arabidopsis thaliana | other | unknown | 0 | 2.84492 | -1.68669 |
| Arabidopsis thaliana | mixed | leaf | 76 | 119.544 | -3.98256 |
| Arabidopsis thaliana | mixed | shoot_meristem | 0 | 1.81948 | -1.34888 |
| Arabidopsis thaliana | mixed | stem | 34 | 15.0511 | 4.88428 |
| Arabidopsis thaliana | mixed | root | 20 | 18.8103 | 0.274318 |
| Arabidopsis thaliana | mixed | root_meristem | 0 | 1.00484 | -1.00242 |
| Arabidopsis thaliana | mixed | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | mixed | flower | 29 | 18.2199 | 2.52551 |
| Arabidopsis thaliana | mixed | fruit | 0 | 0.633409 | -0.79587 |
| Arabidopsis thaliana | mixed | seed | 7 | 6.22463 | 0.31078 |
| Arabidopsis thaliana | mixed | whole | 52 | 57.075 | -0.671753 |
| Arabidopsis thaliana | mixed | in_vitro | 7 | 2.09581 | 3.38759 |
| Arabidopsis thaliana | mixed | other | 0 | 0.898974 | -0.948142 |
| Arabidopsis thaliana | mixed | mixed | 24 | 17.8592 | 1.45308 |
| Arabidopsis thaliana | mixed | unknown | 19 | 8.76365 | 3.45782 |
| Arabidopsis thaliana | unknown | leaf | 40 | 91.442 | -5.37954 |
| Arabidopsis thaliana | unknown | shoot_meristem | 0 | 1.39177 | -1.17973 |
| Arabidopsis thaliana | unknown | stem | 18 | 11.513 | 1.91185 |
| Arabidopsis thaliana | unknown | root | 17 | 14.3884 | 0.688481 |
| Arabidopsis thaliana | unknown | root_meristem | 0 | 0.768628 | -0.876714 |
| Arabidopsis thaliana | unknown | storage | 0 | 0 | 0 |
| Arabidopsis thaliana | unknown | flower | 38 | 13.9369 | 6.44569 |
| Arabidopsis thaliana | unknown | fruit | 0 | 0.48451 | -0.696068 |
| Arabidopsis thaliana | unknown | seed | 12 | 4.76138 | 3.31734 |
| Arabidopsis thaliana | unknown | whole | 55 | 43.6581 | 1.71654 |
| Arabidopsis thaliana | unknown | in_vitro | 0 | 1.60314 | -1.26615 |
| Arabidopsis thaliana | unknown | other | 0 | 0.687648 | -0.829245 |
| Arabidopsis thaliana | unknown | mixed | 8 | 13.661 | -1.53162 |
| Arabidopsis thaliana | unknown | unknown | 17 | 6.70354 | 3.97682 |

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
