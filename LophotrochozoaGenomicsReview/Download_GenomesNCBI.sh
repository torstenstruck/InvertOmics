#!/bin/bash

# (c) Torsten Hugo Struck @05.06.2025, with the assistance of ChatGPT 

# Exit on error, undefined variable check (treat as error), pipeline failure detection
set -euo pipefail

### Usage of pipeline ###
# sh Annotation_v1.sh [Taxon1] [Taxon2] [Taxon3] [and so on]
#
# 1) Taxa names: Names of taxa for which genomes shall be checked and downloaded
#
# example of usage: sh Download_GenomesNCBI.sh Gastrotricha Platyhelminthes Annelida

###### Activate the required conda environment #######
# >>> mamba initialize >>>
# !! Contents within this block are managed by 'mamba shell init' !!
export MAMBA_EXE='/storage/conda/torsths/miniforge3/bin/mamba';
export MAMBA_ROOT_PREFIX='/storage/conda/torsths/miniforge3';
__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias mamba="$MAMBA_EXE"  # Fallback on help from mamba activate
fi
unset __mamba_setup
# <<< mamba initialize <<<

mamba activate base
mamba activate ncbi_datasets

# Set the taxon name for the root folder
for TAXON in "$@"; do

# Create a directory for the taxon if it doesn't exist
mkdir -p "$TAXON"
cd "$TAXON" || exit 1

# Create the annotated genomes list and log file
touch Annotated_Genomes.tsv
touch logfile.txt

# Log everything with echo into logfile.txt instead of STDOUT
log() {
  echo "$@" >> logfile.txt
}

# Start logging
log "=== Processing taxon: $TAXON ==="
log "Started at: $(date)"

# Query genome metadata from NCBI Datasets for the taxon
# Restrict results to: GenBank source, latest version, chromosome/complete assemblies
# Exclude MAGs
# Output as JSON lines → convert to TSV with only organism name and accession fields
datasets summary genome taxon "$TAXON" \
  --assembly-version latest \
  --assembly-level chromosome,complete \
  --mag exclude \
  --report genome \
  --as-json-lines | \
dataformat tsv genome \
  --fields organism-name,accession > "${TAXON}_accessions.tsv"
  
# Input TSV file
INPUT_FILE="${TAXON}_accessions.tsv"
# Output TSV file
OUTPUT_FILE="${TAXON}_accessions_red.tsv"

# Check if the table has only header (i.e., no genome data)
  if [ "$(wc -l < "$INPUT_FILE")" -le 1 ]; then
    log "No chromosome-level or complete genomes available for $TAXON at NCBI"
    cd ..
    continue
  fi


# Use AWK to process the file
# If it is the header line (first line), just copy it to the output and move on
# Extract the accession string (e.g., GCF_000699445.3 or GCA_000699445.3)
# Split accession into parts: e.g., GCF_000699445.3 becomes ["GCF", "000699445", "3"]
# Extract base accession (without prefix or version), e.g., 000699445
# Extract version as a number, e.g., 3
# Check if the accession is GCF (preferred)
# Use the base accession number as the key for tracking
# Logic to decide whether to update the best entry for this base accession:
    # - If this is the first time seeing this key
    # - OR if this is a GCF and the stored one is not
    # - OR if both are GCA, keep the one with the higher version number
        # Save the best line seen so far for this base accession
        # Update stored version and GCF status
# After processing all lines, output the best lines to the output file

awk -F '\t' '
NR == 1 {
    print > "'"$OUTPUT_FILE"'"
    next
}
{
    acc = $2
    split(acc, parts, /[_.]/)
    base = parts[2]
    version = parts[3]
    is_gcf = ($2 ~ /^GCF_/)
    key = base
    if (!(key in best_version) || \
        (is_gcf && !best_is_gcf[key]) || \
        (!best_is_gcf[key] && version + 0 > best_version[key])) {
        best_line[key] = $0
        best_version[key] = version
        best_is_gcf[key] = is_gcf
    }
}
END {
    for (k in best_line)
        print best_line[k] >> "'"$OUTPUT_FILE"'"
}
' "$INPUT_FILE"

# Input TSV file: organism name <TAB> accession
INPUT_TABLE="${TAXON}_accessions_red.tsv"

# Initialize counters
  TOTAL_GENOMES=$(($(wc -l < "$INPUT_TABLE") - 1))  # Minus header
  VALID_GENOMES=0

# Read the TSV, skipping the header
{
  read  # Skip header
  while IFS=$'\t' read -r ORGANISM ACCESSION; do
    SAFE_NAME=$(echo "$ORGANISM" | tr ' ' '_')
    SPECIES_FOLDER="${SAFE_NAME}__${ACCESSION}"

    if [ -d "$SPECIES_FOLDER" ] && [ -d "$SPECIES_FOLDER/$ACCESSION" ]; then
        log "Skipping as already downloaded: $SPECIES_FOLDER"
        continue
    fi

    log "Downloading: $ORGANISM ($ACCESSION)"
    mkdir -p "$SPECIES_FOLDER"

    # Download genome directly into the species folder
    datasets download genome accession "$ACCESSION" --no-progressbar \
        --include gff3,rna,cds,protein,genome,seq-report \
        --filename "$SPECIES_FOLDER/ncbi_dataset.zip"

    unzip -q "$SPECIES_FOLDER/ncbi_dataset.zip" -d "$SPECIES_FOLDER"

    # Flatten directory structure if needed
    if [ -d "$SPECIES_FOLDER/ncbi_dataset/data" ]; then
        mv "$SPECIES_FOLDER"/ncbi_dataset/data/* "$SPECIES_FOLDER"/
        rm -r "$SPECIES_FOLDER/ncbi_dataset"
    fi

    ACCESSION_DIR="$SPECIES_FOLDER/$ACCESSION"
    GFF_FILE="$ACCESSION_DIR/genomic.gff"
    PROTEIN_FILE="$ACCESSION_DIR/protein.faa"

    if [[ -s "$GFF_FILE" && -s "$PROTEIN_FILE" ]]; then
        log "Genome annotated for $ORGANISM ($ACCESSION)"
        echo -e "$ORGANISM\t$ACCESSION" >> Annotated_Genomes.tsv
    else
        log "Genome not annotated for $ORGANISM ($ACCESSION)"
    fi

  done
} < "$INPUT_TABLE"

VALID_GENOMES=$(wc -l < "Annotated_Genomes.tsv") 

log "=== Summary for $TAXON ==="
log "Total genomes listed: $TOTAL_GENOMES"
log "Genomes with annotations: $VALID_GENOMES"
log "Finished at: $(date)"
cd ..

done