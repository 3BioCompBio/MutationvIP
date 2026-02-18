# Data Files Description

## Files Overview

### 1. `clinvar_mutation_context.csv`

ClinVar mutations with context and clinical significance annotations.

**Columns:**

- `chromosome`: Chromosome identifier
- `position`:  HG38 Genomic position of the mutation
- `context`: 11-nucleotide sequence context centered on the mutation site
- `impact`: Clinical significance classification (Likely_benign, Uncertain_significance, Pathogenic)
- `type`: Functional impact type (intronic, missence, silent, splice_region)
- `ref` : ref nucleotide
- `alt`: alt nucleotide

### 2. `cmc_mutation_context.csv`

Cancer Mutation Census (CMC) mutations with context.

**Columns:**

- `chromosome`: Chromosome identifier
- `position`: HG38 Genomic position of the mutation
- `context`: 11-nucleotide sequence context centered on the mutation site
- `type`: Mutation classification (1,2,3 or Other)
- `ref` : ref nucleotide
- `alt`: alt nucleotide

### 3. `cosmic_mutation_context.csv`

COSMIC (Catalogue Of Somatic Mutations In Cancer) mutations with context.

**Columns:**

- `chromosome`: Chromosome identifier
- `position`:  HG38 Genomic position of the mutation
- `context`: 11-nucleotide sequence context centered on the mutation site
- `type`: Functional impact type (missence, intronic, intergenic)
- `cancer_type`: Tissue/organ where cancer originated (breast, pancreas, ovary, liver, etc.)
- `ref`: ref nucleotide
- `alt`: alt nucleotide

### 4. `somamutdb_mutation_context.csv`

Somamutdb mutations with context and diverse tissue/cell type origins.

**Columns:**

- `chromosome`: Chromosome identifier
- `position`:  HG38 Genomic position of the mutation
- `context`: 11-nucleotide sequence context centered on the mutation site
- `type`: Functional location type (intergenic, intronic)
- `cancer_type`: Tissue/cell type of origin (Brain, Colon, iPSC, Bone marrow, Lung, Blood)
- `ref`: ref nucleotide
- `alt`: alt nucleotide

### 5. `experimental.tsv`

Experimental data with electronic properties and mutation density measurements.

**Columns:**

- `chr`: Chromosome identifier
- `start`: Start position of the genomic window
- `end`: End position of the genomic window (200bp windows)
- `RADD` RADD score
- `AP1`: AP score of the first triplicat
- `AP2`: AP score of the second triplicat
- `AP3`: AP score of the third triplicat
- `vIP`: Mean Vertical Ionization Potential
