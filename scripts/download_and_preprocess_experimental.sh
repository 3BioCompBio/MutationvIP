#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# download_and_preprocess_experimental.sh
# Downloads raw experimental files, preprocesses them, and merges results.
#
# Usage: ./download_and_preprocess_experimental.sh [OPTIONS]
#   -f, --force-download     Re-download even if files already exist
#   -F, --force-preprocess   Re-preprocess even if output already exists
#   -a, --force-all          Force all steps
#   -h, --help               Show this help
# ---------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAWDATA_DIR="$SCRIPT_DIR/../rawdata/experimental"

# --- URLs -------------------------------------------------------------------
declare -A AP_URLS=(
    ["GSM3428223_04_APseq_6Gy_30min_rep1_input.bw"]="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3428223&format=file&file=GSM3428223%5F04%5FAPseq%5F6Gy%5F30min%5Frep1%5Finput%2Ebw"
    ["GSM3428224_05_APseq_6Gy_30min_rep2_input.bw"]="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3428224&format=file&file=GSM3428224%5F05%5FAPseq%5F6Gy%5F30min%5Frep2%5Finput%2Ebw"
    ["GSM3428225_06_APseq_6Gy_30min_rep3_input.bw"]="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3428225&format=file&file=GSM3428225%5F06%5FAPseq%5F6Gy%5F30min%5Frep3%5Finput%2Ebw"
)
declare -A RADD_URLS=(
    ["GSE166843_Damage.bedgraph.gz"]="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE166843&format=file&file=GSE166843%5FDamage%2Ebedgraph%2Egz"
)

# --- Flags ------------------------------------------------------------------
FORCE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--force)    FORCE=true ;;
        -h|--help)
            sed -n '4,10p' "$0" | sed 's/^# \?//'
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

# ---------------------------------------------------------------------------
# STEP 1 — Download
# ---------------------------------------------------------------------------
echo "=== Step 1: Download ==="
mkdir -p "$RAWDATA_DIR"

download_file() {
    local filename="$1"
    local url="$2"
    local dest="$RAWDATA_DIR/$filename"

    if [[ -f "$dest" ]] && [[ "$FORCE" == false ]]; then
        echo "  [skip] $filename already exists"
        return
    fi

    echo "  Downloading $filename ..."
    curl -L --fail --progress-bar -o "$dest" "$url"
    echo "  Done: $dest"
}

for filename in "${!AP_URLS[@]}"; do
    download_file "$filename" "${AP_URLS[$filename]}"
done
for filename in "${!RADD_URLS[@]}"; do
    download_file "$filename" "${RADD_URLS[$filename]}"
done

# ---------------------------------------------------------------------------
# STEP 2 — Preprocess: liftover AP files from hg19 → hg38
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 2: Preprocess (liftover AP files hg19 → hg38) ==="

CHAIN_URL="https://raw.githubusercontent.com/liguowang/CrossMap/master/chain_files/human/GRCh37_to_GRCh38.chain.gz"
CHAIN_FILE="$RAWDATA_DIR/GRCh37_to_GRCh38.chain.gz"
FAI="$SCRIPT_DIR/../hg38/genome.fasta.fai"

# --- Dependency checks & installation ---------------------------------------
install_pip_package() {
    local pkg="$1"
    if python -c "import $pkg" &>/dev/null; then return; fi
    echo "  Python package '$pkg' not found."
    read -r -p "  Install it now via 'pip install $pkg'? [y/N] " reply
    if [[ "$reply" =~ ^[Yy]$ ]]; then
        pip install "$pkg"
    else
        echo "  ERROR: '$pkg' is required. Aborting."; exit 1
    fi
}

install_pip_package pyBigWig
install_pip_package pyliftover

# Download chain file if needed
if [[ ! -f "$CHAIN_FILE" ]]; then
    echo "  Downloading chain file..."
    curl -L --fail --progress-bar -o "$CHAIN_FILE" "$CHAIN_URL"
fi

liftover_bw() {
    local input="$1"
    local output="${input%.bw}_hg38.bw"

    if [[ -f "$output" ]] && [[ "$FORCE" == false ]]; then
        echo "  [skip] $(basename "$output") already exists"
        return
    fi

    echo "  Lifting over $(basename "$input") ..."
    python "$SCRIPT_DIR/liftover_bw.py" "$CHAIN_FILE" "$input" "$output"
    echo "  Done: $(basename "$output")"
}

for filename in "${!AP_URLS[@]}"; do
    liftover_bw "$RAWDATA_DIR/$filename"
done

# ---------------------------------------------------------------------------
# STEP 3 — Merge AP replicates with RADD template
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 3: Merge ==="

TEMPLATE="$RAWDATA_DIR/GSE166843_Damage.bedgraph.gz"
OUTPUT="$SCRIPT_DIR/../data/experimental/AP_RADD_merged.tsv"
mkdir -p "$(dirname "$OUTPUT")"

[[ -f "$TEMPLATE" ]] || { echo "  ERROR: template $TEMPLATE not found — re-run step 1"; exit 1; }

if [[ -f "$OUTPUT" ]] && [[ "$FORCE" == false ]]; then
    echo "  [skip] $(basename "$OUTPUT") already exists"
else
    # Collect hg38-lifted AP files in a deterministic order
    BW_FILES=()
    for filename in "${!AP_URLS[@]}"; do
        BW_FILES+=("$RAWDATA_DIR/${filename%.bw}_hg38.bw")
    done

    echo "  Running merge_big_wig.py ..."
    python "$SCRIPT_DIR/merge_big_wig.py" "$TEMPLATE" "$OUTPUT" "${BW_FILES[@]}"
    echo "  Done: $OUTPUT"
fi