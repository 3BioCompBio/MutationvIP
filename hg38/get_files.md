# FASTA FILE
Has been done with.
samtools 1.16.1                                                                                                                                                                                                                                                                                                     
Using htslib 1.16
You can download samtools easily with `sudo apt-get install samtools`
```bash
cd chromosomes
rsync -avzP --exclude='chr*_*.fa.gz' --exclude='chrM.fa.gz' rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/ .
gunzip *.fa.gz
cat chr{1..22}.fa chrX.fa chrY.fa > ../genome.fasta
samtools faidx ../genome.fasta
```

# KMER count file

This script uses [KMC tools](https://github.com/refresh-bio/KMC) - It has been done with version 3.2.4

The script will automatically download KMC 3.2.4 if not present (with user confirmation).

```bash
# Usage: ./get_kmer_count.sh <kmer_size> <input_dir> <output_file>
cd ../scripts
./get_kmer_count.sh 3 ../hg38/chromosomes ../hg38/kmer_counts_3N.csv # already present
./get_kmer_count.sh 5 ../hg38/chromosomes ../hg38/kmer_counts_5N.csv # already present
./get_kmer_count.sh 7 ../hg38/chromosomes ../hg38/kmer_counts_7N.csv
./get_kmer_count.sh 9 ../hg38/chromosomes ../hg38/kmer_counts_9N.csv
./get_kmer_count.sh 11 ../hg38/chromosomes ../hg38/kmer_counts_11N.csv
rm -rf bin/ include/ kmc_tmp/
```

The script will:
1. Check/install KMC tools if needed
2. Count k-mers for each chromosome fasta file
3. Generate a CSV file with k-mer counts per chromosome

# CYTOBAND FILE
```
wget https://gist.githubusercontent.com/alexllc/59e4d221efa402de7509b7bec965cb02/raw/12ab39d867d95fd0e122a9a018c3da134247c3a8/cytoBand-human-GRCh38-hg38.txt
```