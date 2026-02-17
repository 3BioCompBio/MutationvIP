#!/bin/bash

# Usage: ./script_name.sh <kmer> <input_dir> <output_file>

# Check and download KMC if not present
check_and_install_kmc() {
    if ! command -v bin/kmc &> /dev/null && ! [ -f "./bin/kmc" ]; then
        echo "KMC not found."
        read -p "Do you want to download KMC version 3.2.4? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            echo "Downloading KMC version 3.2.4..."
            KMC_VERSION="3.2.4"
            KMC_URL="https://github.com/refresh-bio/KMC/releases/download/v${KMC_VERSION}/KMC3.2.4.linux.x64.tar.gz"

            wget -O kmc.tar.gz "$KMC_URL"
            tar -xzf kmc.tar.gz
            rm kmc.tar.gz

            chmod +x bin/kmc bin/kmc_dump bin/kmc_tools
            echo "KMC version ${KMC_VERSION} installed successfully."
        else
            echo "KMC is required to run this script. Exiting."
            exit 1
        fi
    else
        echo "KMC found."
    fi
}

# Install KMC if needed
check_and_install_kmc

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <kmer> <input_dir> <output_file>"
    exit 1
fi

# Assign arguments to variables
kmer=$1
input_dir=$2
output_file=$3

# Temporary directory for KMC output
tmp_dir="kmc_tmp"
mkdir -p "$tmp_dir"

# Extract chromosome names (without extensions) from fasta files
chromosomes=($(ls "$input_dir"/*.fa | xargs -n 1 basename | sed 's/\.fa//'))

# Create output CSV header
echo -n "context" > "$output_file"
for chr in "${chromosomes[@]}"; do
    if [ "$chr" != "genome" ]; then
        echo -n ",$chr" >> "$output_file"
    fi
done
echo "" >> "$output_file"

# Declare associative array to store k-mer counts
declare -A kmer_counts

# Process each chromosome fasta file
for chr in "${chromosomes[@]}"; do
    if [ "$chr" == "genome" ]; then
        continue
    fi
    input_file="$input_dir/${chr}.fa"
    tmp_kmer_output="${tmp_dir}/${chr}_kmer.tmp"

    # Run KMC
    echo "Processing $chr with k-mer size $kmer..."
    bin/kmc -k${kmer} -b -fm -cs989999999999 -ci1 "$input_file" tmp .
    bin/kmc_tools transform tmp dump "$tmp_kmer_output" -cx99999999999999

    # Read k-mer counts and store them in the associative array
    while IFS=$'\t' read -r kmer_val count; do
        kmer_counts["$kmer_val,$chr"]=$count
    done < "$tmp_kmer_output"

done

# Write k-mer counts to output file
all_kmers=($(printf "%s\n" "${!kmer_counts[@]}" | cut -d',' -f1 | sort -u))
for kmer_val in "${all_kmers[@]}"; do
    echo -n "$kmer_val" >> "$output_file"
    for chr in "${chromosomes[@]}"; do
        if [ "$chr" != "genome" ]; then
            count=${kmer_counts["$kmer_val,$chr"]}
            echo -n ",${count:-0}" >> "$output_file"
        fi
    done
    echo "" >> "$output_file"
done

# Cleanup
rm -rf "$tmp_dir"
echo "Process completed. Output saved in $output_file"
