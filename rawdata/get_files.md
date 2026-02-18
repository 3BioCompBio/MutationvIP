## Experimental Data
Run `scripts/download_and_preprocess_experimental.sh`

- Following samples are selectioned for the AP dammages (GSE121005)
  - GSM3428223
  - GSM3428224
  - GSM3428225
- RADD file can be obtained from GSE166843


## COSMIC signature
Download signature from this page:
- Release Version: Release v3.5
- Variant Class : SBS
- Context Size: 96 
- Download Files: Cosmic_v3.4_SBS_GRCh38

Download the zip and extract `COSMIC_v3.4_SBS_GRCh38.txt` to `rawdata/` subdirectory

## CMC export
Download version 102 of cancer mutation census.
Extract cmc_export.tsv to `rawdata/cmc/cmc_export.tsv`
https://cancer.sanger.ac.uk/cosmic/download/cancer-mutation-census/v102/alldata-cmc