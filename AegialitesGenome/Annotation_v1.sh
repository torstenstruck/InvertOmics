#!/bin/bash

# (c) Torsten Hugo Struck @13.05.2025, based on snakemake from EBP-Nor at 08.05.2025 (https://github.com/ebp-nor/GenomeAnnotation) with the assistance of ChatGPT 

# Exit on error, undefined variable check (treat as error), pipeline failure detection
set -euo pipefail 

### Usage of pipeline ###
# sh Annotation_v1.sh [Prefix] [Softmasked assembly] [Proteins from model organism] [Proteins from Uniprot] [Proteins from OrthoDB] [Additional BUSCO Lineage] [Minimum AA size] [Uniprot database] [Maximum number of threads] > FilteredDatasets_QualityAssessment_Pipeline.out
#
# 1) Prefix: Name to be used a prefix to the final files generated
# 2) Softmasked assembly: Path to and name of multi-line fasta file of the softmasked assembly to be annotated (if single-line, the programm will generate a multi-line one)
# 3) Proteins from model organism: Path to and name of fasta file with proteins from model organism
# 4) Proteins from Uniprot: Path to and name of fasta file with proteins from uniprot
# 5) Proteins from OrthoDB: Path to and name of fasta file with proteins from OrthoDB for a specific taxonomic unit
# 6) Additional BUSCO Lineage: Name of one additional BUSCO lineage (metazoa_odb12 is done anyway) to be used for BUSCO analysis, for example, lophotrochozoa_odb12, mollusca_odb12, arthropoda_odb12 or coleoptera_odb12
# 7) Minimum AA size: Minimum ORF size that gene models must pass
# 8) Uniprot database: Path to and name of database from uniprot
# 9) Maximum number of threads: The maximum number of threads that should be used (e.g., 64) IMPORTANT NOTE: It must be an even number!!!
#
# example of usage: sh Annotation_v1.sh Aeg55_final_asm_anno \
#										 /storage/InvertOmics/03_RepeatMasking/Arthropoda/Aegialites_55802/128cpus/Aeg55_final_pri_asm.fasta.masked \
#										 /storage/InvertOmics/04_Annotation/Arthropoda/Resources/Tribolium_castaneum/ncbi_dataset/data/GCF_031307605.1/protein.faa \
#										 /storage/InvertOmics/00_Scripts/01_Databases/Funannotate_DBs/uniprot_sprot.fasta \
#										 /storage/InvertOmics/00_Scripts/01_Databases/OrthoDB_v12/Arthropoda.fa \
#										 coleoptera_odb12 50 \
#										 /storage/InvertOmics/00_Scripts/01_Databases/Funannotate_DBs/uniprot.dmnd \
#										 64 > Aegialites_55802_Annotation.out
#
# IMPORTANT NOTE: The parameters have to be provided in this order and a value for each parameter has to be provided!!!
# IMPORTANT NOTE: If you do changes to the script ALWAYS save the file under a different name. Do NOT modify this master script.

####### Activate the required conda environment #######
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
mamba activate annotation

####### Configuration of variables and folders #######
# Variables to be set via the command line
echo "Setting variables ..." > annotation.log
# Name used for the identifying the analyses
PREFIX=$1 
echo "Variable PREFIX is ${PREFIX}." >> annotation.log
# Softmasked assembly
ASSEMBLY=$2 
echo "Variable ASSEMBLY is ${ASSEMBLY}." >> annotation.log
# Proteins from a model species (zebrafish for fish, human for mammals, etc.)
MODEL_PROTEINS=$3 
echo "Variable MODEL_PROTEINS is ${MODEL_PROTEINS}." >> annotation.log
# SwissProt/UniProtKB database as protein sequences in fasta file
UNIPROT_PROTEINS=$4 
echo "Variable UNIPROT_PROTEINS is ${UNIPROT_PROTEINS}." >> annotation.log
# OrthoDB v12 clade of proteins (precompiled or individually compiled)
ORTHODB_PROTEINS=$5 
echo "Variable ORTHODB_PROTEINS is ${ORTHODB_PROTEINS}." >> annotation.log
# Provide one additional BUSCO database (metazoa_odb12 is done anyway), for example, lophotrochozoa_odb12, mollusca_odb12, arthropoda_odb12 or coleoptera_odb12
BUSCO_LINEAGE=$6
echo "Variable BUSCO_LINEAGE is ${BUSCO_LINEAGE}." >> annotation.log
# Minimum ORF size that gene models must pass (e.g., min_aa_size = 50).
MIN_AA_SIZE=$7 
echo "Variable MIN_AA_SIZE is ${MIN_AA_SIZE}." >> annotation.log
# SwissProt/UniProtKB database 
UNIPROT_DB=$8
echo "Variable UNIPROT_DB is ${UNIPROT_DB}." >> annotation.log
# Number of threads to be used
THREADS=$9 
echo "Variable THREADS is ${THREADS}." >> annotation.log
# Setting home directory
HOME_DIR="$PWD"
echo "Variable HOME_DIR is ${HOME_DIR}." >> annotation.log

# Make working directories
echo "Making directories ..." >> annotation.log
rm -rf galba
mkdir -p proteins miniprot galba

####### Step 1: Preprocessing #######
# Check for soft-masking (lowercase bases)
echo "Checking assembly file ..." >> annotation.log
if ! grep -q '[acgtn]' "$ASSEMBLY"; then
    echo "ERROR: Assembly must be soft-masked (contain lowercase bases for repeats)."
    exit 1
fi

# Count number of header lines and sequence lines
HEADER_COUNT=$(grep -c '^>' "$ASSEMBLY")
SEQ_LINE_COUNT=$(grep -v '^>' "$ASSEMBLY" | wc -l)

# If the fasta file is single-lined, break the sequence lines in units of 80 bp while leaving the headers unchangded
# If it is multi-lined the original is copied to ensure that all assemblies have the same name
if [[ "$SEQ_LINE_COUNT" -eq "$HEADER_COUNT" ]]; then
    echo "Detected single-line FASTA. Folding to 80 characters per line..."
    awk '
        /^>/ {print; next}
        {
            while (length > 0) {
                print substr($0, 1, 80)
                $0 = substr($0, 81)
            }
        }
    ' "$ASSEMBLY" > ${PREFIX}.fold.fa
else
    echo "Detected multi-line FASTA. Copying as is..."
    cp "$ASSEMBLY" ${PREFIX}.fold.fa
fi

# Preparing the input data for the different model searches
echo "Preparing input data ..." >> annotation.log
cp "$MODEL_PROTEINS" proteins/model.fasta
ln -s "$UNIPROT_PROTEINS" proteins/uniprot.fasta
ln -s "$UNIPROT_DB" proteins/uniprot.dmnd
ln -s "$ORTHODB_PROTEINS" proteins/orthodb.fasta

####### Step 2: Miniprot alignments and extracting protein and mRNA sequences as well as preparing GFF3 files for EvidenceModeler #######
for protein_set in model uniprot orthodb; do
	echo "Running miniprot on ${protein_set} ..." >> annotation.log
    miniprot -It${THREADS} --gff ${PREFIX}.fold.fa proteins/${protein_set}.fasta > miniprot/${protein_set}_mp_aln.gff 2> miniprot/miniprot_${protein_set}.err
    if [[ "$protein_set" != "orthodb" ]]; then
        agat_sp_extract_sequences.pl --gff miniprot/${protein_set}_mp_aln.gff -f ${PREFIX}.fold.fa -p -o miniprot/${protein_set}_mp.proteins.fa
        agat_sp_extract_sequences.pl --gff miniprot/${protein_set}_mp_aln.gff -f ${PREFIX}.fold.fa -t exon --merge -o miniprot/${protein_set}_mp.mrna.fa
    else
        touch miniprot/orthodb_mp.proteins.fa
        touch miniprot/orthodb_mp.mrna.fa
    fi
    python "$EVM_HOME"/EvmUtils/misc/miniprot_GFF_2_EVM_GFF3.py miniprot/${protein_set}_mp_aln.gff > miniprot/${protein_set}_mp4evm.gff3
	if [[ ! -s miniprot/${protein_set}_mp4evm.gff3 ]]; then
    echo "ERROR: ${protein_set}_mp4evm.gff3 was not created." >> annotation.log
    exit 1
fi

done

####### Step 3: Run GALBA for ab initio protein-model prediction #######
echo "Running galba ..." >> annotation.log
singularity exec -B "$HOME_DIR":/data /storage/conda/torsths/singularities/galba/galba.sif galba.pl --version > galba.version
singularity exec -B "$HOME_DIR":/data /storage/conda/torsths/singularities/galba/galba.sif cp -rf /opt/Augustus/config /data/galba
singularity exec -B "$HOME_DIR":/data /storage/conda/torsths/singularities/galba/galba.sif galba.pl --species="$PREFIX" --threads="$THREADS" --genome=/data/${PREFIX}.fold.fa --prot_seq=/data/proteins/model.fasta --workingdir=/data/galba --augustus_args="--stopCodonExcludedFromCDS=False" --gff3 --AUGUSTUS_CONFIG_PATH=/data/galba/config
# Wait/check for GTF output
echo "Waiting for galba.gtf to be ready..." >> annotation.log
for i in {1..10}; do
    if [[ -s galba/galba.gtf ]]; then
        echo "galba.gtf is ready." >> annotation.log
        break
    fi
    echo "${i} - Still waiting..." >> annotation.log
    sleep 60
done

if [[ ! -s galba/galba.gtf ]]; then
    echo "ERROR: galba.gtf not found or empty after waiting." >> annotation.log
    exit 1
fi

# Extracting protein sequences as well as preparing GFF3 file for EvidenceModeler
echo "Extracting galba sequences ..." >> annotation.log
agat_sp_extract_sequences.pl --gff galba/galba.gtf -f ${PREFIX}.fold.fa -t cds -p -o galba/galba.proteins.fa

if [[ ! -s galba/galba.proteins.fa ]]; then
    echo "ERROR: galba.proteins.fa was not created." >> annotation.log
    exit 1
fi

####### Step 4: EvidenceModeler (EVM) to generate gene models from all predictions above #######
# Make working directories
echo "Making directories ..." >> annotation.log
rm -rf filter
mkdir -p prediction evm filter ipr uniprot functional

# Convert Galba gtf to evm.gff3 file 
echo "Converting galba files ..." >> annotation.log
"$EVM_HOME"/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl galba/galba.gtf > galba/galba.evm.gff3

if [[ ! -s galba/galba.evm.gff3 ]]; then
    echo "ERROR: evm/galba.evm.gff3 was not created." >> annotation.log
    exit 1
fi

# Combine all gff3 files from the four different approaches
echo "Preparing files for EVM ..." >> annotation.log
cat galba/galba.evm.gff3 miniprot/model_mp4evm.gff3 miniprot/uniprot_mp4evm.gff3 miniprot/orthodb_mp4evm.gff3 | awk '/gene/ {printf "\n"} 1' > prediction/gene_predictions.gff3

if [[ ! -s prediction/gene_predictions.gff3 ]]; then
    echo "ERROR: gene_predictions.gff3 was not created." >> annotation.log
    exit 1
fi

printf "ABINITIO_PREDICTION\tAugustus\t1\n" > prediction/weights.evm.txt 
printf "OTHER_PREDICTION\tminiprot\t4\n" >> prediction/weights.evm.txt

if [[ ! -s prediction/weights.evm.txt ]]; then
    echo "ERROR: weights.evm.txt was not created." >> annotation.log
    exit 1
fi

echo "Running EVM ..." >> annotation.log
/storage/conda/torsths/miniforge3/envs/annotation/lib/python3.9/site-packages/funannotate/aux_scripts/funannotate-runEVM.py -w prediction/weights.evm.txt -c "$THREADS" -d evm -g prediction/gene_predictions.gff3 -f ${PREFIX}.fold.fa -l evm.log -m 10 -i 1500 -o evm/evm.gff3 --EVM_HOME "$EVM_HOME"

if [[ ! -s evm/evm.gff3 ]]; then
    echo "ERROR: evm.gff3 was not created." >> annotation.log
    exit 1
fi

# Extract protein and mRNA sequences
echo "Converting EVM files ..." >> annotation.log
agat_sp_extract_sequences.pl --gff evm/evm.gff3 -f ${PREFIX}.fold.fa -t cds -p -o evm/evm.proteins.fa
agat_sp_extract_sequences.pl --gff evm/evm.gff3 -f ${PREFIX}.fold.fa -t exon --merge -o evm/evm.mrna.fa

if [[ ! -s evm/evm.proteins.fa ]]; then
    echo "ERROR: evm.proteins.fa was not created." >> annotation.log
    exit 1
fi

####### Step 5: Filter out low-quality gene perdictions #######
# Detect repeat-like proteins using DIAMOND
echo "Filtering protein detection ..." >> annotation.log
diamond blastp --sensitive --query evm/evm.proteins.fa --db /storage/InvertOmics/00_Scripts/01_Databases/Funannotate_DBs/repeats.dmnd --evalue 1e-10 --max-target-seqs 1 --out filter/repeats.tsv --outfmt 6 --threads "$THREADS"

if [[ ! -s filter/repeats.tsv ]]; then
    echo "ERROR: repeats.tsv was not created." >> annotation.log
    exit 1
fi

# Add matching protein IDs into a kill list
cut -f 1 filter/repeats.tsv | sort -u > filter/kill.list

# Detect and add proteins with "X" (e.g., translation across gaps) to kill list
echo "Detect genes with X ..." >> annotation.log
awk '
    /^>/ {
        if (seq ~ /X/) print header
        header=$0; seq=""
        next
    }
    {
        seq = seq $0
    }
    END {
        if (seq ~ /X/) print header
    }
' evm/evm.proteins.fa | cut -d " " -f1 | sed 's/^>//' >> filter/kill.list 2>> kill.err

if [[ ! -s filter/kill.list ]]; then
    echo "ERROR: kill.list was not created." >> annotation.log
    exit 1
fi

# Remove GFF entries corresponding to "kill list"
echo "Removing low quality genes ..." >> annotation.log
agat_sp_filter_feature_from_kill_list.pl --gff evm/evm.gff3 --kill_list filter/kill.list -o filter/removed_repeats.gff

if [[ ! -s filter/removed_repeats.gff ]]; then
    echo "ERROR: removed_repeats.gff was not created." >> annotation.log
    exit 1
fi

# Filter out short ORFs by size
agat_sp_filter_by_ORF_size.pl --gff filter/removed_repeats.gff -s "$MIN_AA_SIZE" -o filter/filtered_genes.gff

if [[ ! -s filter/filtered_genes_sup${MIN_AA_SIZE}.gff ]]; then
    echo "ERROR: filtered_genes_sup${MIN_AA_SIZE}.gff was not created." >> annotation.log
    exit 1
fi

# Extract cleaned protein and mRNA sequences
agat_sp_extract_sequences.pl --gff filter/filtered_genes_sup${MIN_AA_SIZE}.gff -f ${PREFIX}.fold.fa -t cds -p -o filter/filtered.proteins.fa
agat_sp_extract_sequences.pl --gff filter/filtered_genes_sup${MIN_AA_SIZE}.gff -f ${PREFIX}.fold.fa -t exon --merge -o filter/filtered.mrna.fa

if [[ ! -s filter/filtered.proteins.fa ]]; then
    echo "ERROR: filtered.proteins.fa was not created." >> annotation.log
    exit 1
fi

####### Step 6: InterProScan to annotate filtered protein sequences with known protein domains, GO terms, and family classifications #######
mamba deactivate
mamba deactivate

echo "Running InterProScan ..." >> annotation.log
module --quiet purge
module load InterProScan_data/5.62-94.0-foss-2023a
# Clean the input FASTA as stop codons might cause parsing problems
cat filter/filtered.proteins.fa | sed 's/*//g' > ipr/input.proteins.fa

if [[ ! -s ipr/input.proteins.fa ]]; then
    echo "ERROR: input.proteins.fa was not created." >> annotation.log
    exit 1
fi

# Run the interproscan run with goterms, three output formats and specifying a temporary directory 
cd ipr
interproscan.sh --cpu "$THREADS" -dp -i input.proteins.fa  -f XML,GFF3,TSV -goterms -iprlookup -b ipr 1> interproscan.out 2> interproscan.err 
cd ..

if [[ ! -s ipr/ipr.tsv ]]; then
    echo "ERROR: ipr.tsv was not created." >> annotation.log
    exit 1
fi
module --quiet purge

####### Step 7: UniProt BLAST to functionally annotate filtered proteins by aligning them to known proteins in the UniProt database #######
mamba activate base
mamba activate annotation

echo "Running UniProt Blast ..." >> annotation.log
diamond blastp --query filter/filtered.proteins.fa --db proteins/uniprot.dmnd --out uniprot/diamond.blastp.out --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 --threads "$THREADS"

if [[ ! -s uniprot/diamond.blastp.out ]]; then
    echo "ERROR: diamond.blastp.out was not created." >> annotation.log
    exit 1
fi

####### Step 8: Merge functional annotations #######
# Merge annotations from InterProScan and UniProt Blast using a perl script
echo "Merging functional annotation ..." >> annotation.log
agat_sp_manage_functional_annotation.pl --gff filter/filtered_genes_sup${MIN_AA_SIZE}.gff -i ipr/ipr.tsv -b uniprot/diamond.blastp.out --ID FUNC -o functional/fun --clean_name -db proteins/uniprot.fasta

if [[ ! -s functional/fun/filtered_genes_sup${MIN_AA_SIZE}.gff ]]; then
    echo "ERROR: filtered_genes_sup${MIN_AA_SIZE}.gff was not created." >> annotation.log
    exit 1
fi

cp functional/fun/filtered_genes_sup${MIN_AA_SIZE}.gff functional/${PREFIX}.gff

if [[ ! -s functional/${PREFIX}.gff ]]; then
    echo "ERROR: ${PREFIX}.gff was not created." >> annotation.log
    exit 1
fi

# Generate summary stats on the funcationally annotated genes
echo "Generate summary statistics ..." >> annotation.log
agat_sp_statistics.pl --gff functional/${PREFIX}.gff -o functional/gff_stats

# Extract cleaned protein and mRNA sequences
echo "Extract final sequences ..." >> annotation.log
agat_sp_extract_sequences.pl --gff functional/${PREFIX}.gff -f ${PREFIX}.fold.fa -t cds -p -o functional/${PREFIX}.proteins.fa
agat_sp_extract_sequences.pl --gff functional/${PREFIX}.gff -f ${PREFIX}.fold.fa -t exon --merge -o functional/${PREFIX}.mrna.fa

if [[ ! -s functional/${PREFIX}.proteins.fa ]]; then
    echo "ERROR: ${PREFIX}.proteins.fa was not created." >> annotation.log
    exit 1
fi

mamba deactivate
mamba deactivate

####### Step 9: BUSCO analysis of the filtered genome against the metazoa orthodb database #######
# Determine BUSCO score for metazoa_odb12
module --quiet purge
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4 
export AUGUSTUS_CONFIG_PATH=/storage/ConfigFiles/BRAKER2/config/
export BUSCO_CONFIG_FILE=/storage/ConfigFiles/BUSCO/config.ini

echo "Running BUSCO ..." >> annotation.log
cp functional/${PREFIX}.proteins.fa functional/${PREFIX}.mrna.fa .

busco -i ${PREFIX}.proteins.fa -l metazoa_odb10 -o BUSCO_metazoa_odb10 -m proteins > Busco_metazoa_odb10.log 

# Determine BUSCO score for another selected database
busco -i ${PREFIX}.proteins.fa -l "$BUSCO_LINEAGE" -o BUSCO_${BUSCO_LINEAGE} -m proteins > Busco_${BUSCO_LINEAGE}.log 

module --quiet purge

####### Final clean up #######
rm -rf ipr/temp
