#!/bin/sh

# (c) Torsten Hugo Struck @04.05.2022 

set -o errexit

### Usage of pipeline ###
# sh AssemblyHifiAsm_QualityAssessment_v1.7.sh [Name of fastq file] [Name of barcode file] [e value] [Maximum number of hits] [Maximum number of threads] [RAM for Merquery analysis] [Genome size] > AssemblyHifiAsm_QualityAssessment_Pipeline.out
#
# 1) Name of fastq file: name of the fastq or fasta to be used for the assembly
# 2) Name of the barcode: name of the barcode file to find contigs with the barcode sequences (e.g., mitochondrial genes), sequence must be amino acids
# 3) e value: the e value to be used for the different blast searches (e.g., 1e-20)
# 4) Maximum number of hits: the maximum number of hits to be saved in the blast searches (e.g., 1)
# 5) Maximum number of threads: the maximum number of threads that should be used (e.g., 64) IMPORTANT NOTE: It must be an even number!!!
# 6) RAM: RAM needed for Merquery analysis in T or G (e.g., 1T)
# 7) Genome size: the estimated or known genome size in K, M or G (e.g., 926M)
#
# example of usage: sh AssemblyHifiAsm_QualityAssessment_v1.7.sh seqdata.fastq Query_Polychaeta_ProtMito.fasta 1e-20 1 64 1T 926M > AssemblyHifiAsm_QualityAssessment_Pipeline.out
#
# IMPORTANT NOTE: The parameters have to be provided in this order and a value for each parameter has to be provided!!!
# IMPORTANT NOTE: If you do changes to the script ALWAYS save the file under a different name. Do NOT modify this master script.
# IMPORTANT NOTE: HifiAsm saves the primary assembly as a gfa file, but the pipeline will also generate a fasta file of this named "${FASTQ}_pri_asm.contigs.fasta"; the two haplotigs are saved as "${FASTQ}_hap1|hap2_asm.contigs.fasta"
# IMPORTANT NOTE: Recovery points have still to be implemented.

### Version history ###
# v1.1: fixing minor bugs
# v1.2: changed blobtools to BlobToolKit v3.1.0; haploidy analyses included using HapPy; additional purging of duplicates included using purge_dups; QC on both purged and unpurged assemblies except for species detection; fixing minor bugs
# v1.3: separated purge_dups into its own pipeline; fixing minor bugs
# v1.4: changed BLAST+ 2.11.0 to 2.13.0
# v1.5: Changed to Blobtools 3.1.4; fixed an issue with the wait commands
# v1.6: Upgraded to HifiAsm0.18; generate fasta file of alternative assembly; added Merquery to the pipeline; fixed some minor text issues
# v1.7: Commented out HapPy analyses due to re-occuring error messages (it needs to be check thoroughly checked); fixed a bug concerning BUSCO scores and Blobtool analyses

### Variables ###
FASTQ=$1
BARCODE=$2
EVALUE=$3
MAXHITS=$4
THREADS=$5
REDTHREADS=$((${THREADS}/2))
MEMORY=$6
GENOMESIZE=$7
SECONDS=0

# print out settings for the pipeline to the log
echo -e "The pipeline was called as follows:" > AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "sh AssemblyHifiAsm_QualityAssessment_v1.5.sh ${FASTQ} ${BARCODE} ${EVALUE} ${MAXHITS} ${THREADS} ${MEMORY} ${GENOMESIZE}" >> AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "Name of fastq file: ${FASTQ}" >> AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "Name of barcode query file: ${BARCODE}" >> AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "e value: ${EVALUE}" >> AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "Maximum number of saved blast hits: ${MAXHITS}" >> AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "Maximum number of threads: ${THREADS}" >> AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "Maximum size of RAM: ${MEMORY}" >> AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "Genome size: ${GENOMESIZE}" >> AssemblyHifiAsm_QualityAssessment_Pipeline.log
DATE=$(date)
echo -e "\nPipeline was started: $DATE" >> AssemblyHifiAsm_QualityAssessment_Pipeline.log 

### Different subscripts used during the script ###

function RetrieveTaxonomyDatabase () {
# (c) Frédéric Mahé & modified by Torsten Hugo Struck @12.05.2021 adjusting to new databases in NCBI

## Download NCBI's taxonomic data and GI (GenBank ID) taxonomic
## assignation.

## Variables
local NCBI="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
local TAXDUMP="taxdump.tar.gz"
local TAXIDPATH="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
local TAXID="nucl_gb.accession2taxid.gz"
local NAMES="names.dmp"
local NODES="nodes.dmp"
local DMP=$(echo {citations,division,gencode,merged,delnodes}.dmp)
local USELESS_FILES="${TAXDUMP} ${DMP} gc.prt readme.txt"

## Download taxdump
rm -rf ${USELESS_FILES} "${NODES}" "${NAMES}"
wget "${NCBI}${TAXDUMP}" && \
    tar zxvf "${TAXDUMP}" && \
    rm -rf ${USELESS_FILES}

## Limit search space to scientific names
grep "scientific name" "${NAMES}" > "${NAMES/.dmp/_reduced.dmp}" && \
    rm -f "${NAMES}" && \
    mv "${NAMES/.dmp/_reduced.dmp}" "${NAMES}"

## Download gi_taxid_nucl
rm -f "${TAXID/.gz/}*"
wget "${TAXIDPATH}${TAXID}" && \
    gunzip "${TAXID}"
}

function RetrieveTaxonPath () {
# (c) Frédéric Mahé & modified by Torsten Hugo Struck @12.05.2021 adjusting to new databases in NCBI

local NAMES="names.dmp"
local NODES="nodes.dmp"
local GI_TO_TAXID="nucl_gb.accession2taxid"
local TAXONOMY=""
local GI="${1}"

# Obtain the name corresponding to a taxid or the taxid of the parent taxa
get_name_or_taxid()
{
    grep --max-count=1 "^${1}"$'\t' "${2}" | cut -f"${3}"
}

# Get the taxid corresponding to the GI number
TAXID=$(get_name_or_taxid "${GI}" "${GI_TO_TAXID}" "3")

# Loop until you reach the root of the taxonomy (i.e. taxid = 1)
while [[ "${TAXID}" -gt 1 ]] ; do
    # Obtain the scientific name corresponding to a taxid
    NAME=$(get_name_or_taxid "${TAXID}" "${NAMES}" "3")
    # Obtain the parent taxa taxid
    PARENT=$(get_name_or_taxid "${TAXID}" "${NODES}" "3")
    # Build the taxonomy path
    TAXONOMY="${NAME};${TAXONOMY}"
    TAXID="${PARENT}"
done

echo -e "${GI}\t${TAXONOMY}"
}

function SpeciesTaxaDetectionPipeline () {
# Set variables with the relevant paths
local ASSEMBLY=$1
echo "Settings of the species detection pipeline:	Assembly file: ${ASSEMBLY}	Query file: ${BARCODE}	e value: ${EVALUE}	max. hits: ${MAXHITS}" > Settings_Detection_Pipeline_${FASTQ}_${ASSEMBLY}.txt

#Generate blast-searchable database
makeblastdb -in ${ASSEMBLY} -parse_seqids -dbtype nucl

#Blast the barcode sequence(s) against the assembly library
tblastn -evalue ${EVALUE} -outfmt "6 qseqid sseqid pident score evalue" -num_threads ${THREADS} -out ${ASSEMBLY}_Species_BLASThits.txt -db ${ASSEMBLY} -query ${BARCODE}


#Retrieve the corresponding sequences to a hit and write it out to new fasta files
cut -f2 < ${ASSEMBLY}_Species_BLASThits.txt | sort | uniq | while read ID
do
echo ${ID}
sed -n "/^>${ID}/,/>/p" < ${ASSEMBLY} | sed '$d' >> ${ASSEMBLY}_Sequences_of_best_hit.fas
done

#Blast the retrieved sequences against NCBI nr database remotely (version 2.6 is needed for that)
blastn -max_target_seqs ${MAXHITS} -outfmt "6 sacc qseqid sseqid length pident score evalue" -out ${ASSEMBLY}_nr_BLAST_Matches.txt -db nt -num_threads ${THREADS}  -query ${ASSEMBLY}_Sequences_of_best_hit.fas

#Determine the taxonomic position of the hit
echo "Retrieving taxa from blast results"
cut -f1 < ${ASSEMBLY}_nr_BLAST_Matches.txt | sort | uniq | while read GI
do
RetrieveTaxonPath ${GI} >> ${ASSEMBLY}_Taxa_found.txt
done
}

function Map_reads_to_assembly () {
# Set the local variables
local ASSEMBLY=$1

module load minimap2/2.17-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
minimap2 -ax map-pb -t 16  ${ASSEMBLY} ../${FASTQ} | samtools view -buS - | samtools sort - -o MappedReads_${ASSEMBLY}.sorted.bam
samtools index MappedReads_${ASSEMBLY}.sorted.bam 
samtools stats MappedReads_${ASSEMBLY}.sorted.bam > MappedReads_${ASSEMBLY}.sorted.bam.stats
grep ^COV MappedReads_${ASSEMBLY}.sorted.bam.stats | cut -f 2- > MappedReads_${ASSEMBLY}.sorted.bam.stats.coverage
module load gnuplot/5.2.8-GCCcore-8.3.0
plot-bamstats -p MappedReads_${ASSEMBLY}.sorted.bam.sorted.bam.stats.vis MappedReads_${ASSEMBLY}.sorted.bam.stats
echo -e "Mapping ended at time $(date +%T)."
}

function BlobtoolAnalyses () {
# Set the local variables
local ASSEMBLY=$1

module load BLAST+/2.13.0-gompi-2020a 

# Generate Blast hit table for Blobtool analyses and map the reads to the assembly
echo -e "Blast search and mapping started at $(date +%T)"
blastn -task megablast -query ${ASSEMBLY} -db nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs ${MAXHITS} -max_hsps 1 -num_threads ${REDTHREADS} -evalue ${EVALUE} -out ${ASSEMBLY}.blobplot.blast.out &
local PID1=$!
Map_reads_to_assembly ${ASSEMBLY} &
local PID2=$!
wait $PID1
echo -e "Blast search ended at time $(date +%T)."
wait $PID2

### Run actual blobtool analyses ###
module load Blobtools/3.1.4

#Create the BlobtoolsResults folder containing the assembly
blobtools create --fasta ${ASSEMBLY} BlobtoolsResults_Primary
#Add blast hit information to the BlobtoolsResults folder
blobtools add --hits ${ASSEMBLY}.blobplot.blast.out --taxrule bestsumorder --taxdump ./taxdump BlobtoolsResults_Primary
#Add mapping coverage to the BlobtoolsResults folder
blobtools add --cov MappedReads_${ASSEMBLY}.sorted.bam BlobtoolsResults_Primary
#Add BUSCO scores to the BlobtoolsResults folder
blobtools add --busco ./BUSCO_Metazoa_Primary/run_metazoa_odb10/full_table.tsv BlobtoolsResults_Primary
}

function HapPyAnalyses () {
# Set the local variables
local ASSEMBLY=$1

module load HapPy/0.2.1-foss-2020a

mkdir -p HaploidyResults_Primary
cd HaploidyResults_Primary

#HapPy commands to determine haploidy
happy coverage --outdir=./HapPy_coverage/ ../MappedReads_${ASSEMBLY}.sorted.bam
happy estimate --size $GENOMESIZE --plot --outstats HapPy_estimate_${ASSEMBLY}.txt ./HapPy_coverage/MappedReads_${ASSEMBLY}.sorted.bam.hist
happy autoest --size=$GENOMESIZE --plot --outstats HapPy_autoest_${ASSEMBLY}.txt ./HapPy_coverage/MappedReads_${ASSEMBLY}.sorted.bam.hist

cd ..
}

### Run assembly using HifiAsm ###
echo -e "\nRunning the assembly using HifiAsm with following command line: hifiasm -o ${FASTQ}_asm -t ${THREADS} ../${FASTQ}" >> AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\nStep 1: Running HifiAsm assembly"

module --quiet purge  # Reset the modules to the system default
module load hifiasm/0.18.2-GCCcore-10.3.0

# set HifiAsm settings
mkdir -p HifiAsm_assembly
cd HifiAsm_assembly
hifiasm -o ${FASTQ}_asm -t ${THREADS} ../${FASTQ}
awk '/^S/{print ">"$2;print $3}' ${FASTQ}_asm.bp.p_ctg.gfa > ${FASTQ}_pri_asm.contigs.fasta
ASM_UNPUR="${FASTQ}_pri_asm.contigs.fasta"
awk '/^S/{print ">"$2;print $3}' ${FASTQ}_asm.bp.hap1.p_ctg.gfa > ${FASTQ}_hap1_asm.contigs.fasta
ASM_HAP1="${FASTQ}_hap1_asm.contigs.fasta"
awk '/^S/{print ">"$2;print $3}' ${FASTQ}_asm.bp.hap2.p_ctg.gfa > ${FASTQ}_hap2_asm.contigs.fasta
ASM_HAP2="${FASTQ}_hap2_asm.contigs.fasta"

### Quality assessment using Quast ###
DURA_ASS_SEC=$SECONDS
DURA_ASS=$(($DURA_ASS_SEC/3600))
DATE=$(date)
echo -e "Duration of HifiAsm analysis in hours: ${DURA_ASS}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\nRunning QUAST analyses on the primary assembly and the two haplotigs" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 2: Running QUAST analyses on the primary assembly"

module --quiet purge
module load QUAST/5.0.2-foss-2020a-Python-3.8.2 

# run quast on primary assmebly
quast.py -o Quast_Primary ${ASM_UNPUR} 
quast.py -o Quast_Haplotig1 ${ASM_HAP1} 
quast.py -o Quast_Haplotig2 ${ASM_HAP2} 

### Quality assessment using Quast and Merquery ###
echo -e "\nRunning Merqury analyses on the primary and alternative assembly" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 3: Running Merqury analyses on the primary and alternative assembly"

module --quiet purge
module load merqury/1.3

#build meryl database for the merquery analysis
meryl k=21 threads=${THREADS} memory=100g count output ${FASTQ}.meryl ../${FASTQ}

#define path to merqury srcipts
export MERQURY=/opt/software/custom/software/merqury-1.3/

#run merqury analysis
sh /opt/software/custom/software/merqury-1.3/merqury.sh ${FASTQ}.meryl ${ASM_UNPUR} UnpurgedPri > UnpurgedPri_merqury.out 2> UnpurgedPri_merqury.err
sh /opt/software/custom/software/merqury-1.3/merqury.sh ${FASTQ}.meryl ${ASM_HAP1} ${ASM_HAP2} UnpurgedHaps > UnpurgedHaps_merqury.out 2> UnpurgedHaps_merqury.err

### Species identifications using provided barcode sequences ###
DURA_QUAST_SEC=$SECONDS
DURA_QUAST=$((($DURA_QUAST_SEC-$DURA_ASS_SEC)/3600))
DATE=$(date)
echo -e "Duration of QUAST and Merqury analyses in hours: ${DURA_QUAST}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\nRunning barcoding detection on the primary assembly" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "The following settings have been used for the BLAST analyses:" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "makeblastdb -in ${ASM_UNPUR} -parse_seqids -dbtype nucl" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "tblastn -evalue ${EVALUE} -outfmt \"6 qseqid sseqid pident score evalue\" -num_threads ${THREADS} -out ${ASM_UNPUR}_Species_BLASThits.txt -db ${ASM_UNPUR} -query ${BARCODE}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "blastn -max_target_seqs ${MAXHITS} -outfmt \"6 sacc qseqid sseqid length pident score evalue\" -out ${ASM_UNPUR}_nr_BLAST_Matches.txt -db nt -num_threads ${THREADS}  -query ${ASM_UNPUR}_Sequences_of_best_hit.fas" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 4: Running barcoding detection"
module --quiet purge
module load BLAST+/2.13.0-gompi-2020a 
mv ../${BARCODE} .

# run the species detection script retrieving all contigs and unassembled reads matching the provided query sequences as well as the corresponding taxon name
RetrieveTaxonomyDatabase
perl /storage/software/easybuild/software/BLAST+/2.13.0-gompi-2020a/bin/update_blastdb.pl --passive --decompress nt
SpeciesTaxaDetectionPipeline ${ASM_UNPUR}

### BUSCO check on primary assembly and haplotig #1 ###
DURA_MIT_SEC=$SECONDS
DURA_MIT=$((($DURA_MIT_SEC-$DURA_QUAST_SEC)/3600))
DATE=$(date)
echo -e "Duration of barcoding detection in hours: ${DURA_MIT}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 5: Running the BUSCO analysis"
module --quiet purge
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4 
export AUGUSTUS_CONFIG_PATH=/storage/ConfigFiles/BRAKER2/config/
export BUSCO_CONFIG_FILE=/storage/ConfigFiles/BUSCO/config.ini

# set the BUSCO settings
echo -e "\nRunning BUSCO analyses on the primary assembly using the following command line: busco -i ${ASM_UNPUR} -l metazoa_odb10 --out BUSCO_Metazoa_Primary -m geno -c ${THREADS} -f > Busco_${ASM_UNPUR}.log" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
busco -i ${ASM_UNPUR} -l metazoa_odb10 --out BUSCO_Metazoa_Primary -m geno -c ${THREADS} -f > Busco_${ASM_UNPUR}.log
echo -e "\nRunning BUSCO analyses on the haplotig #1 using the following command line: busco -i ${ASM_HAP1} -l metazoa_odb10 --out BUSCO_Metazoa_Haplotig1 -m geno -c ${THREADS} -f > Busco_${ASM_BAP1}.log" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
busco -i ${ASM_HAP1} -l metazoa_odb10 --out BUSCO_Metazoa_Haplotig1 -m geno -c ${THREADS} -f > Busco_${ASM_HAP1}.log

### Run blobtool analyses ###
### load necessary for modules for the preparation steps ###
DURA_BUSCO_SEC=$SECONDS
DURA_BUSCO=$((($DURA_BUSCO_SEC-$DURA_MIT_SEC)/3600))
DATE=$(date)
echo -e "Duration of BUSCO analysis in hours: ${DURA_BUSCO}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\nRunning the Blobtools analysis on the primary assembly using the following command lines:" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "blastn -task megablast -query ${ASM_UNPUR} -db nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs ${MAXHITS} -max_hsps 1 -num_threads ${REDTHREADS} -evalue ${EVALUE} -out ${ASM_UNPUR}.blobplot.blast.out" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "blobtools create --fasta ${ASM_UNPUR} BlobtoolsResults_${ASM_UNPUR}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "blobtools add --hits ${ASM_UNPUR}.blobplot.blast.out --taxrule bestsumorder --taxdump ./taxdump BlobtoolsResults_${ASM_UNPUR}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "blobtools add --cov MappedReads_${ASM_UNPUR}.sorted.bam BlobtoolsResults_${ASM_UNPUR}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "blobtools add --busco ./BUSCO_Metazoa_${ASM_UNPUR}/run_metazoa_odb10/full_table.tsv BlobtoolsResults_${ASM_UNPUR}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "To display the data in a webserver do the following:" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "blobtools view --remote ./BlobtoolsResults_${ASM_UNPUR}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "follow the instructions at https://github.com/blobtoolkit/tutorials/blob/main/2021-04-28_blobtoolkit_commandline_mollusc.md for \"Run a viewer for a BlobDir\", start at \"Run the blobtools view command:\"" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 6: Running the Blobtools analysis"

#Fetch the NCBI Taxdump
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;

module --quiet purge
BlobtoolAnalyses ${ASM_UNPUR}

#Background reading and help:
#https://blobtoolkit.genomehubs.org/blobtools2/
#https://www.youtube.com/channel/UC6QimXygwjaX2CGZ8cgjNxw
#https://github.com/blobtoolkit/tutorials/blob/main/2021-04-28_blobtoolkit_commandline_mollusc.md

#To display the data in a webserver do the following
#blobtools view --remote ./BlobtoolsResults_${ASM_UNPUR} or ./BlobtoolsResults_${ASM_PUR} 
#follow the instructions in the last link above for "Run a viewer for a BlobDir", start at "Run the blobtools view command:"

### Run HapPy analyses ###
# the bam files from the BlobToolKit are also needed here
DURA_BLOB_SEC=$SECONDS
DURA_BLOB=$((($DURA_BLOB_SEC-$DURA_BUSCO_SEC)/3600))
DATE=$(date)
echo -e "Duration of the Blobtools analysis in hours: ${DURA_BLOB}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
#echo -e "\nRunning the HapPy analysis on the primary assembly using the following command lines:" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
#echo -e "happy coverage --outdir=./HapPy_coverage/ MappedReads_${ASM_UNPUR}.sorted.bam" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
#echo -e "happy estimate --size ${GENOMESIZE} --plot --outstats HapPy_estimate_${ASM_UNPUR}.txt ./HapPy_coverage/MappedReads_${ASM_UNPUR}.sorted.bam.hist" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
#echo -e "happy autoest --size=${GENOMESIZE} --plot --outstats HapPy_autoest_${ASM_UNPUR}.txt ./HapPy_coverage/MappedReads_${ASM_UNPUR}.sorted.bam.hist" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
#echo -e "\n${DATE}"
#echo -e "\nStep 7: Running the HapPy analysis"

#module --quiet purge
#HapPyAnalyses ${ASM_UNPUR}


### Cleaning up folder ###
#DURA_HAPPY_SEC=$SECONDS
#DURA_HAPPY=$((($DURA_HAPPY_SEC-$DURA_BLOB_SEC)/3600))
#DATE=$(date)
#echo -e "Duration of the HapPy analysis in hours: ${DURA_HAPPY}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\nCleaning up the data and moving results into appropriate folders" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 8: Cleaning up the data and moving results"
rm -rf nt.*
rm -rf *.dmp nodesDB.txt nucl_gb.accession2taxid
rm -rf busco_downloads/
rm -rf busco_data
rm -rf taxdump/
rm -rf taxdb.*

### Sorting results ###
mkdir -p MitochondrialResults
mv *_nr_BLAST_Matches.txt *_Sequences_of_best_hit.fas *_Species_BLASThits.txt *_Taxa_found.txt ${BARCODE} Settings_Detection_Pipeline_* MitochondrialResults/
mkdir -p BAM_results
mv MappedReads_* BAM_results/
mv Busco_${ASM_UNPUR}.log BUSCO_Metazoa_Primary/
mv Busco_${ASM_HAP1}.log BUSCO_Metazoa_Haplotig1/
mkdir Merquery_results
mv UnpurgedPri* UnpurgedHaps* Merquery_results/

DURA_ALL=$(($SECONDS/3600))
DATE=$(date)
echo -e "Duration of all analyses in hours: ${DURA_ALL}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\nThe assembly and quality assessment pipeline ended at ${DATE}" >> ../AssemblyHifiAsm_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nAnalyses finished."
