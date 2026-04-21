#!/usr/bin/env bash
#
# Download the GRCh38 reference genome (no-alt + hs38d1 decoys).
# Required for:
#   - Running CRAM fixture integration tests (cargo test)
#   - Regenerating CRAM test fixtures (build_test_fixtures)
#   - Browser testing with CRAM files
#
# The reference is ~3 GB and takes a few minutes to download.
# Files are placed in ref/ which is gitignored.
#
# Usage:
#   ./scripts/download-reference.sh
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REF_DIR="${REPO_ROOT}/ref"
FASTA="${REF_DIR}/GRCh38_no_alt_plus_hs38d1.fna"
FAI="${FASTA}.fai"

NCBI_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"

mkdir -p "${REF_DIR}"

if [ -f "${FASTA}" ]; then
    echo "Reference already exists: ${FASTA}"
    echo "Size: $(du -h "${FASTA}" | cut -f1)"

    if [ ! -f "${FAI}" ]; then
        echo "Index missing. Creating .fai..."
        samtools faidx "${FASTA}"
        echo "Done: ${FAI}"
    fi
    exit 0
fi

echo "Downloading GRCh38 reference genome (~900 MB compressed, ~3 GB uncompressed)..."
echo "Source: NCBI"
echo ""

if command -v curl &>/dev/null; then
    curl -L --progress-bar -o "${FASTA}.gz" "${NCBI_URL}"
elif command -v wget &>/dev/null; then
    wget --show-progress -O "${FASTA}.gz" "${NCBI_URL}"
else
    echo "ERROR: curl or wget required" >&2
    exit 1
fi

echo ""
echo "Decompressing..."
gunzip "${FASTA}.gz"

echo "Creating FASTA index (.fai)..."
if command -v samtools &>/dev/null; then
    samtools faidx "${FASTA}"
    echo "Index created: ${FAI}"
else
    echo "WARNING: samtools not found. Install samtools to create .fai index."
    echo "  brew install samtools   # macOS"
    echo "  apt install samtools    # Ubuntu/Debian"
fi

echo ""
echo "Done. Reference genome: ${FASTA}"
echo "Size: $(du -h "${FASTA}" | cut -f1)"
