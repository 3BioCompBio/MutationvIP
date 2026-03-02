> All the files in this directory can be obtained from Zenodo (DOI: 10.5281/zenodo.17870586)

## COSMIC signatures

### cosmic_sbs_etiologies.tsv
Etiologies for COSMIC SBS signatures, collected from the COSMIC website. Each signature was classified by hand; justifications are provided in the supplementary materials of the paper.

### COSMIC_v3.4_SBS_GRCh38.txt
SBS mutational signatures from COSMIC v3.4 (GRCh38, 96-trinucleotide context).

---

## Mutation context datasets

Single-base substitutions were downloaded in June 2025 from COSMIC (GenomeScreensMutant v102), ClinVar, and SomaMutDB and filtered to retain only SNVs on canonical chromosomes. Each variant was annotated with snpEff v5.2e (GRCh38.mane.1.2.refseq) and classified into seven mutually exclusive functional categories using the mapping in the table below. When multiple annotations were present, the highest-priority (most functionally impactful) category was retained. The 11-nucleotide genomic context centered on each variant was extracted from GRCh38.

| snpEff consequence | Category |
|---|---|
| missense_variant | missense |
| start_lost | missense |
| stop_lost | missense |
| initiator_codon_variant | missense |
| stop_gained | nonsense |
| synonymous_variant | silent |
| start_retained_variant | silent |
| stop_retained_variant | silent |
| splice_region_variant | splice_region |
| splice_donor_variant | splice_region |
| splice_acceptor_variant | splice_region |
| 5_prime_UTR_variant | UTR |
| 3_prime_UTR_variant | UTR |
| 5_prime_UTR_premature_start_codon_gain_variant | UTR |
| intron_variant | intronic |
| intragenic_variant | — (excluded) |
| upstream_gene_variant | intergenic |
| downstream_gene_variant | intergenic |
| intergenic_region | intergenic |

### COSMIC — cosmic_mutation_context.csv
Somatic SNVs from COSMIC GenomeScreensMutant v102 (GRCh38). Cancer type (PRIMARY_SITE) was assigned by linking each mutation to its phenotype record. Columns: chromosome, position, context, type, cancer_type, ref, alt.

### ClinVar — clinvar_mutation_context.csv
SNVs from ClinVar (GRCh38, June 2025). Clinical significance (CLNSIG) was harmonised to five classes: Benign, Likely_benign, Uncertain_significance, Likely_pathogenic, Pathogenic. Columns: chromosome, position, context, type, impact, ref, alt.

### SomaMutDB — somamutdb_mutation_context.csv
Somatic SNVs from SomaMutDB (June 2025). Columns: chromosome, position, context, type, cancer_type, ref, alt.

---

## CMC

### cmc_mutation_context.csv
For each somatic substitution in the COSMIC Cancer Mutation Census (CMC, v102), the 11-nucleotide genomic context centered on the mutated position was retrieved from the GRCh38 human reference genome. Only single-nucleotide substitutions mapping to canonical chromosomes were retained, along with their significance tier (MUTATION_SIGNIFICANCE_TIER), reference allele, and alternate allele.
q

---

## Experimental

### experimental.tsv
`experimental.tsv` can be obtained from Zenodo and also be computed directly.