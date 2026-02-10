
# GAMA  
**Genomic Availability & Metadata Analysis Tool**

> ⚠️ **Development Status: Active / Early Stage**  
> GAMA is currently in early development. Interfaces, outputs, and scoring methods may change.  
> Users are encouraged to validate results independently and report issues.

GAMA is an R-based framework for surveying publicly available sequencing data across NCBI Assembly, SRA, and BioSample. Its aim is to support feasibility assessments for in silico research on underutilised plant species.

---

## Overview

GAMA:

- Unifies NCBI database searches  
- Computes a data richness score  
- Classifies SRA accessions by experimental modality  
- Enables strategic parsing of Assembly and SRA accession metadata  
- Generates publication-ready visuals 

---

## Installation

Install the development version from GitHub using `pak`:

``` r
install.packages('pak')
pak::pak('JLewis-dev/GAMA')
```

---

## Quick-Start Example

### 1. Load package

``` r
library(GAMA)
```

### 2. Configure NCBI access (recommended)

To improve rate limits and ensure responsible use of NCBI services:

``` r
options(ENTREZ_EMAIL = 'your.email@example.com')
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
META <- summarise_sra_availability(RESULTS)
print(META)
```

### 7. Visualise SRA modality composition

``` r
plot_sra_availability(META)
```

### 8. Extract filtered Assembly accession metadata

``` r
ASM_META <- extract_assembly_metadata(RESULTS, best = TRUE)
print(ASM_META)
```

### 9. Extract filtered SRA accession metadata

``` r
SRA_META <- extract_sra_metadata(RESULTS, species = 'Vigna vexillata', class = 'genomic')
print(SRA_META)
```

### 10. Cite

``` r
citation('GAMA')
```

---

## Visualisation Tools

GAMA includes built-in plotting functions for rapid assessment.

### Data richness plots

`plot_availability()` produces stacked bar plots showing:

- Assembly contribution  
- SRA contribution  
- BioSample contribution  
- Overall data richness score  

Supports custom ranking, colour palettes, and ggplot2 theming.

---

### SRA modality plots

`plot_sra_availability()` visualises:

- Relative abundance of sequencing strategies  
- Ontology-classified experiment types  
- Cross-species comparisons  

Optional GEO overlays are available via extended plotting functions.

---

## Transparency

Queries automatically record:

- Timestamp  
- Tool version  
- Database sources  
- Search terms  

This metadata is embedded in outputs.

---

## Data Richness

The data richness score is defined as:

Score = A + S + B

Where A, S, and B are the transformed contributions of Assembly, SRA, and BioSample accession counts.

A = best + ln(1 + total − best), with assemblies weighted as:

- Complete = 10  
- Chromosome = 8  
- Scaffold = 5  
- Contig = 2  

Here, best is the maximum weight assembly (ties broken by highest N50) and total is the sum of all accession weights.

S = 2·ln(1 + SRA)  
B = ln(1 + BioSample)

This formulation prioritises high-quality assemblies while incorporating diminishing returns for extensively sampled taxa.

---

## Ontology-Driven Classification

SRA experiments are classified using an ontology derived from large-scale metadata mining and manual curation.

### Genomic

- WGS  
- Amplicon-seq  
- RAD-seq  
- Targeted-Capture  
- Clone-based  

### Transcriptomic

- RNA-seq  
- small-RNA  
- Long-read  

### Epigenomic

- Bisulfite-seq  
- ChIP-seq  
- CUT&RUN  
- CUT&Tag  
- ATAC-seq  
- DNase-seq  
- FAIRE-seq  
- MNase-seq  
- SELEX  

### Chromatin

- Hi-C  
- 3C-based  
- ChIA-PET  
- TCC  

### Other

- Other  

Fallback strategies are applied when primary metadata fields are missing or ambiguous.

---

## Intended Use

GAMA is designed for:

- Grant and project scoping  
- Identification of under-studied taxa  
- Strategic prioritisation of existing data sets  

It is particularly suited to investigations of underutilised and non-model plant species.

---

## Limitations

- Dependent on NCBI metadata quality  
- Runtime increases with species list size  
- Novel protocols may not be fully captured by the ontology  
- Results should be interpreted cautiously during early development  

---

## Licence

See the `LICENSE` file for details.

---
