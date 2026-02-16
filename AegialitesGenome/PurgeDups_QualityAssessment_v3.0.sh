#!/bin/sh

# (c) Torsten Hugo Struck @04.05.2022 

set -o errexit

### Usage of pipeline ###
# sh PurgeDups_QualityAssessment_v3.0.sh [Name of fastq file] [Name of primary assembly file] [Name of haplotig1 file] [Name of haplotig2 file] [Name to meryl database] [e value] [Maximum number of hits] [Maximum number of threads] >> PurgeDups_QualityAssessment_Pipeline.out
#
# 1) Name of fastq file: name of the fastq or fasta to be used for the original assembly (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 2) Name of primary assembly file: name of the fasta file to be used for the purging of duplicates from the primary assembly (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 3) Name of haplotig1 file: name of the fasta file to be used for the purging of duplicates from haplotig1 (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 4) Name of haplotig2 file: name of the fasta file to be used for the purging of duplicates from haplotig2 (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 5) Name of meryl database: name of the meryl database (.meryl) to be used for the merqury analyses (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 6) e value: the e value to be used for the different blast searches (e.g., 1e-20)
# 7) Maximum number of hits: the maximum number of hits to be saved in the blast searches (e.g., 1)
# 8) Maximum number of threads: the maximum number of threads that should be used (e.g., 64) IMPORTANT NOTE: It must be an even number!!!
#
# example of usage: sh PurgeDups_QualityAssessment_v3.0.sh ../reads.fastq pri_assembly.fasta hap1_assembly.fasta hap2_assembly.fasta reads.fastq.meryl 1e-20 1 64 >> PurgeDups_QualityAssessment_Pipeline.out
#
# IMPORTANT NOTE: The parameters have to be provided in this order and a value for each parameter has to be provided!!!
# IMPORTANT NOTE: If you do changes to the script ALWAYS save the file under a different name. Do NOT modify this master script.
# IMPORTANT NOTE: Recovery points have still to be implemented.

### Version history ###
# v1.1: Changed from Blast+ 2.11.0 to 2.13.0
# v1.2: Changed to Blobtools 3.1.4; fixed an issue with the wait commands
# v1.3: Added Purge_dups on haplotigs; added Merqury analysis
# v1.4: Commented out HapPy analyses due to re-occuring error messages (it needs to be check thoroughly checked); fixed a bug concerning BUSCO scores and Blobtool analyses
# v2.0: Took out HapPy analyses, fixed some minor bugs in reporting
# v2.1: Instead of downloading a local copy each time for the BLAST search, a local copy is now download at regular intervals and the BLAST search is done against this copy
# v3.0: Added Blobtool analyses also for the haplotig assemblies, added BUSCO analysis for the second haplotig analysis, removed all code relating to the HapPy analyses, fixed some minor issues with dates and reorder moving files to folders


### Variables ###
FASTQ=$1
ASS=$2
HAP1=$3
HAP2=$4
MERYL=$5
EVALUE=$6
MAXHITS=$7
THREADS=$8
REDTHREADS=$((${THREADS}/2))
SECONDS=0

# print out settings for the pipeline to the log
echo -e "The pipeline was called as follows:" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "sh PurgeDups_QualityAssessment_v2.1.sh ${FASTQ} ${ASS} ${HAP1} ${HAP2} ${MERYL} ${EVALUE} ${MAXHITS} ${THREADS}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Name of the fastq file: ${FASTQ}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Name of the fasta file of the primary assembly: ${ASS}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Name of the fasta file of the haplotig1: ${HAP1}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Name of the fasta file of the haplotig2: ${HAP2}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Name of the meryl database: ${MERYL}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "e value of blast search: ${EVALUE}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Maximum number of hits in blast searcg: ${MAXHITS}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Maximum number of threads: ${THREADS}" >> PurgeDups_QualityAssessment_Pipeline.log
DATE=$(date)
echo -e "\nPipeline was started: $DATE" >> PurgeDups_QualityAssessment_Pipeline.log 

### Different subscripts used during the script ###
function Map_reads_to_assembly_PAF () {
local ASSEMBLY=$1
minimap2 -xmap-pb -t 16 $ASSEMBLY ${FASTQ} | gzip -c - > ${FASTQ}.paf.gz
pbcstat ${FASTQ}.paf.gz
calcuts PB.stat > cutoffs 2> calcults.log
echo -e "Mapping ended at time $(date +%T).\n"
}

function SelfMap_assembly_PAF () {
local ASSEMBLY=$1
split_fa $ASSEMBLY > ${ASSEMBLY}.split
minimap2 -xasm5 -DP -t 16 ${ASSEMBLY}.split ${ASSEMBLY}.split | gzip -c - > ${ASSEMBLY}.split.self.paf.gz
echo -e "Self mapping ended at time $(date +%T).\n"
}

function Purge_dups_func () {
local ASSEMBLY=$1

module --quiet purge  # Reset the modules to the system default
module load Purge_dups/1.0.1

#Run minimap2 to align pacbio data and generate paf files, then calculate read depth histogram and base-level read depth.
echo -e "Self mapping and mapping started at $(date +%T)"
Map_reads_to_assembly_PAF $ASSEMBLY &
local PID1=$!
#Split an assembly and do a self-self alignment.
SelfMap_assembly_PAF $ASSEMBLY &
local PID2=$!
wait $PID1
wait $PID2


#Purge haplotigs and overlaps.
purge_dups -2 -T cutoffs -c PB.base.cov $ASSEMBLY.split.self.paf.gz > dups.bed 2> purge_dups.log

#Get purged primary and haplotig sequences from draft assembly.
get_seqs -p $ASSEMBLY dups.bed $ASSEMBLY
}

function Map_reads_to_assembly () {
# Set the local variables
local ASSEMBLY=$1

module load minimap2/2.17-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
minimap2 -ax map-pb -t 16  ${ASSEMBLY} ${FASTQ} | samtools view -buS - | samtools sort - -o MappedReads_${ASSEMBLY}.sorted.bam
samtools index MappedReads_${ASSEMBLY}.sorted.bam 
samtools stats MappedReads_${ASSEMBLY}.sorted.bam > MappedReads_${ASSEMBLY}.sorted.bam.stats
grep ^COV MappedReads_${ASSEMBLY}.sorted.bam.stats | cut -f 2- > MappedReads_${ASSEMBLY}.sorted.bam.stats.coverage
module load gnuplot/5.2.8-GCCcore-8.3.0
plot-bamstats -p MappedReads_${ASSEMBLY}.sorted.bam.sorted.bam.stats.vis MappedReads_${ASSEMBLY}.sorted.bam.stats
echo -e "Mapping ended at time $(date).\n"
}

function BlobtoolAnalyses () {
# Set the local variables
local ASSEMBLY=$1

module load BLAST+/2.13.0-gompi-2020a 

# Generate Blast hit table for Blobtool analyses and map the reads to the assembly
echo -e "Blast search and mapping started at $(date)"
blastn -task megablast -query ${ASSEMBLY} -db /storage/InvertOmics/00_Scripts/NCBI_nt/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs ${MAXHITS} -max_hsps 1 -num_threads ${REDTHREADS} -evalue ${EVALUE} -out ${ASSEMBLY}.blobplot.blast.out &
local PID1=$!
Map_reads_to_assembly ${ASSEMBLY} &
local PID2=$!
wait $PID1
echo -e "Blast search ended at time $(date).\n"
wait $PID2

### Run actual blobtool analyses ###
module load Blobtools/3.1.4

#Create the BlobtoolsResults folder containing the assembly
blobtools create --fasta ${ASSEMBLY} BlobtoolsResults_Purged
#Add blast hit information to the BlobtoolsResults folder
blobtools add --hits ${ASSEMBLY}.blobplot.blast.out --taxrule bestsumorder --taxdump ./taxdump BlobtoolsResults_Purged
#Add mapping coverage to the BlobtoolsResults folder
blobtools add --cov MappedReads_${ASSEMBLY}.sorted.bam BlobtoolsResults_Purged
#Add BUSCO scores to the BlobtoolsResults folder
blobtools add --busco ./BUSCO_Metazoa_Purged/run_metazoa_odb10/full_table.tsv BlobtoolsResults_Purged
}


### Purge duplicates ###
echo -e "\nPurging duplicates" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\nStep 1: Purging duplicates"

#Purging duplicates on the primary assembly
Purge_dups_func ${ASS}
mv ${ASS}.purged.fa ${ASS}.purged.fasta
ASM_PUR="${ASS}.purged.fasta"
mkdir -p PurgeDups_Primary
mv PB.* purge_dups.log calcults.log cutoffs dups.bed PurgeDups_Primary/

#Purging duplicates on the haplotig1
Purge_dups_func ${HAP1} 
mv ${HAP1}.purged.fa ${HAP1}.purged.fasta
HAP1_PUR="${HAP1}.purged.fasta"
mkdir -p PurgeDups_Haplotig1
mv PB.* purge_dups.log calcults.log cutoffs dups.bed PurgeDups_Haplotig1/

#Purging duplicates on the haplotig2
Purge_dups_func ${HAP2} 
mv ${HAP2}.purged.fa ${HAP2}.purged.fasta
HAP2_PUR="${HAP2}.purged.fasta"
mkdir -p PurgeDups_Haplotig2
mv PB.* purge_dups.log calcults.log cutoffs dups.bed PurgeDups_Haplotig2/

### Quality assessment using Quast on purged primary and haplotig data ###
DURA_PURGE_SEC=$SECONDS
DURA_PURGE=$(($DURA_PURGE_SEC/3600))
DATE=$(date)
echo -e "Duration of purging duplicates in hours: ${DURA_PURGE}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\nRunning QUAST analyses on the purged files" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 2: Running QUAST analyses on the purged data"

module --quiet purge
module load QUAST/5.0.2-foss-2020a-Python-3.8.2 

# run quast quast on all assemblies
quast.py -o Quast_Purged_Primary ${ASM_PUR}
quast.py -o Quast_Purged_Haplotig1 ${HAP1_PUR}
quast.py -o Quast_Purged_Haplotig2 ${HAP2_PUR}

### Quality assessment using Merquery on purged primary and haplotig assembly###
echo -e "\nRunning Merqury analyses on the purged on the primary and haplotig assemblies" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 3: Running Merqury analyses on the purged primary and alternative assembly"

module --quiet purge
module load merqury/1.3

#define path to merqury srcipts
export MERQURY=/opt/software/custom/software/merqury-1.3/

#run merqury analysis
sh /opt/software/custom/software/merqury-1.3/merqury.sh ${MERYL} ${ASM_PUR} PurgedPri > PurgedPri_merqury.out 2> PurgedPri_merqury.err
sh /opt/software/custom/software/merqury-1.3/merqury.sh ${MERYL} ${HAP1_PUR} ${HAP2_PUR} PurgedHaps > PurgedHaps_merqury.out 2> PurgedHaps_merqury.err
mkdir Merquery_results_Purged
mv PurgedPri* PurgedHaps* Merquery_results_Purged/

### BUSCO and Blobtool analyses ###
### load necessary for modules for the preparation steps ###

#Fetch the NCBI Taxdump
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;

### BUSCO and Blobtool analyses on the purged primary assembly ###
### BUSCO check on purged primary assembly###
DURA_QUAST_SEC=$SECONDS
DURA_QUAST=$((($DURA_QUAST_SEC-$DURA_PURGE_SEC)/3600))
DATE=$(date)
echo -e "Duration of QUAST and Merqury analyses in hours: ${DURA_QUAST}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 4a: Running the BUSCO analysis on the purged primary assembly"
module --quiet purge
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4 
export AUGUSTUS_CONFIG_PATH=/storage/ConfigFiles/BRAKER2/config/
export BUSCO_CONFIG_FILE=/storage/ConfigFiles/BUSCO/config.ini

# set the BUSCO settings
echo -e "\nRunning BUSCO analyses on the purged primary assembly using the following command line: busco -i ${ASM_PUR} -l metazoa_odb10 --out BUSCO_Metazoa_Purged -m geno -c ${THREADS} -f > Busco_${ASM_PUR}.log" >> PurgeDups_QualityAssessment_Pipeline.log
busco -i ${ASM_PUR} -l metazoa_odb10 --out BUSCO_Metazoa_Purged -m geno -c ${THREADS} -f > Busco_${ASM_PUR}.log

### Run blobtool analyses on the purged primary assembly###
DURA_BUSCO_SEC=$SECONDS
DURA_BUSCO=$((($DURA_BUSCO_SEC-$DURA_QUAST_SEC)/3600))
DATE=$(date)
echo -e "Duration of BUSCO analysis on the purged primary assembly in hours: ${DURA_BUSCO}" >> PurgeDups_QualityAssessment_Pipeline.log

echo -e "\nRunning the Blobtools analysis on the purged primary assembly using the following command lines:" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blastn -task megablast -query ${ASM_PUR} -db /storage/InvertOmics/00_Scripts/NCBI_nt/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs ${MAXHITS} -max_hsps 1 -num_threads ${REDTHREADS} -evalue ${EVALUE} -out ${ASM_PUR}.blobplot.blast.out" >> PurgeDups_QualityAssessment_Pipeline.log

echo -e "\nblobtools create --fasta ${ASM_PUR} BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --hits ${ASM_PUR}.blobplot.blast.out --taxrule bestsumorder --taxdump ./taxdump BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --cov MappedReads_${ASM_PUR}.sorted.bam BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --busco ./BUSCO_Metazoa_Purged/run_metazoa_odb10/full_table.tsv BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log

echo -e "\n${DATE}"
echo -e "\nStep 5a: Running the Blobtools analysis on the purged primary assembly"
module --quiet purge
BlobtoolAnalyses ${ASM_PUR}

mv BlobtoolsResults_Purged/ BlobtoolsResults_Purged_Primary/
mv BUSCO_Metazoa_Purged/ BUSCO_Metazoa_Purged_Primary/
mv Busco_${ASM_PUR}.log BUSCO_Metazoa_Purged_Primary/
rm -rf busco_downloads/
rm -rf busco_data
mkdir -p BAM_results_Purged_Primary
mv MappedReads_* BAM_results_Purged_Primary/

echo -e "\n\n********** To display the data in a webserver do the following: **********" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools view --remote ./BlobtoolsResults_Purged_Primary" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "follow the instructions at https://github.com/blobtoolkit/tutorials/blob/main/2021-04-28_blobtoolkit_commandline_mollusc.md for \"Run a viewer for a BlobDir\", start at \"Run the blobtools view command:\"" >> PurgeDups_QualityAssessment_Pipeline.log

### BUSCO and Blobtool analyses on the purged haplotig1 assembly ###
### BUSCO check on purged haplotig1 assembly###
DURA_BLOB_SEC=$SECONDS
DURA_BLOB=$((($DURA_BLOB_SEC-$DURA_BUSCO_SEC)/3600))
DATE=$(date)
echo -e "Duration of of the Blobtools analysis on the purged primary assembly in hours: ${DURA_QUAST}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 4b: Running the BUSCO analysis on the purged haplotig1 assembly"
module --quiet purge
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4 
export AUGUSTUS_CONFIG_PATH=/storage/ConfigFiles/BRAKER2/config/
export BUSCO_CONFIG_FILE=/storage/ConfigFiles/BUSCO/config.ini

# set the BUSCO settings
echo -e "\nRunning BUSCO analyses on the purged haplotig1 using the following command line: busco -i ${HAP1_PUR} -l metazoa_odb10 --out BUSCO_Metazoa_Purged -m geno -c ${THREADS} -f > Busco_${HAP1_PUR}.log" >> PurgeDups_QualityAssessment_Pipeline.log
busco -i ${HAP1_PUR} -l metazoa_odb10 --out BUSCO_Metazoa_Purged -m geno -c ${THREADS} -f > Busco_${HAP1_PUR}.log

### Run blobtool analyses on the purged haplotig1 assembly###
DURA_BUSCO_SEC=$SECONDS
DURA_BUSCO=$((($DURA_BUSCO_SEC-$DURA_BLOB_SEC)/3600))
DATE=$(date)
echo -e "Duration of BUSCO analysis on the purged haplotig1 assembly in hours: ${DURA_BUSCO}" >> PurgeDups_QualityAssessment_Pipeline.log

echo -e "\nRunning the Blobtools analysis on the purged haplotig1 assembly using the following command lines:" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blastn -task megablast -query ${HAP1_PUR} -db /storage/InvertOmics/00_Scripts/NCBI_nt/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs ${MAXHITS} -max_hsps 1 -num_threads ${REDTHREADS} -evalue ${EVALUE} -out ${HAP1_PUR}.blobplot.blast.out" >> PurgeDups_QualityAssessment_Pipeline.log

echo -e "\nblobtools create --fasta ${HAP1_PUR} BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --hits ${HAP1_PUR}.blobplot.blast.out --taxrule bestsumorder --taxdump ./taxdump BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --cov MappedReads_${HAP1_PUR}.sorted.bam BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --busco ./BUSCO_Metazoa_Purged/run_metazoa_odb10/full_table.tsv BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log

echo -e "\n${DATE}"
echo -e "\nStep 5b: Running the Blobtools analysis on the purged haplotig1 assembly"
module --quiet purge
BlobtoolAnalyses ${HAP1_PUR}

mv BlobtoolsResults_Purged/ BlobtoolsResults_Purged_Haplotig1/
mv BUSCO_Metazoa_Purged/ BUSCO_Metazoa_Purged_Haplotig1/
mv Busco_${HAP1_PUR}.log BUSCO_Metazoa_Purged_Haplotig1/
rm -rf busco_downloads/
rm -rf busco_data
mkdir -p BAM_results_Purged_Haplotig1
mv MappedReads_* BAM_results_Purged_Haplotig1/

echo -e "\n\n********** To display the data in a webserver do the following: **********" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools view --remote ./BlobtoolsResults_Purged_Haplotig1" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "follow the instructions at https://github.com/blobtoolkit/tutorials/blob/main/2021-04-28_blobtoolkit_commandline_mollusc.md for \"Run a viewer for a BlobDir\", start at \"Run the blobtools view command:\"" >> PurgeDups_QualityAssessment_Pipeline.log

### BUSCO and Blobtool analyses on the purged haplotig2 assembly ###
### BUSCO check on purged haplotig2 assembly###
DURA_BLOB_SEC=$SECONDS
DURA_BLOB=$((($DURA_BLOB_SEC-$DURA_BUSCO_SEC)/3600))
DATE=$(date)
echo -e "Duration of of the Blobtools analysis on the purged haplotig1 assembly in hours: ${DURA_QUAST}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 4c: Running the BUSCO analysis on the purged haplotig2 assembly"
module --quiet purge
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4 
export AUGUSTUS_CONFIG_PATH=/storage/ConfigFiles/BRAKER2/config/
export BUSCO_CONFIG_FILE=/storage/ConfigFiles/BUSCO/config.ini

# set the BUSCO settings
echo -e "\nRunning BUSCO analyses on the purged haplotig2 using the following command line: busco -i ${HAP2_PUR} -l metazoa_odb10 --out BUSCO_Metazoa_Purged -m geno -c ${THREADS} -f > Busco_${HAP2_PUR}.log" >> PurgeDups_QualityAssessment_Pipeline.log
busco -i ${HAP2_PUR} -l metazoa_odb10 --out BUSCO_Metazoa_Purged -m geno -c ${THREADS} -f > Busco_${HAP2_PUR}.log

### Run blobtool analyses on the purged haplotig2 assembly###
DURA_BUSCO_SEC=$SECONDS
DURA_BUSCO=$((($DURA_BUSCO_SEC-$DURA_BLOB_SEC)/3600))
DATE=$(date)
echo -e "Duration of BUSCO analysis on the purged haplotig2 assembly in hours: ${DURA_BUSCO}" >> PurgeDups_QualityAssessment_Pipeline.log

echo -e "\nRunning the Blobtools analysis on the purged haplotig2 assembly using the following command lines:" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blastn -task megablast -query ${HAP2_PUR} -db /storage/InvertOmics/00_Scripts/NCBI_nt/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs ${MAXHITS} -max_hsps 1 -num_threads ${REDTHREADS} -evalue ${EVALUE} -out ${HAP2_PUR}.blobplot.blast.out" >> PurgeDups_QualityAssessment_Pipeline.log

echo -e "\nblobtools create --fasta ${HAP2_PUR} BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --hits ${HAP2_PUR}.blobplot.blast.out --taxrule bestsumorder --taxdump ./taxdump BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --cov MappedReads_${HAP2_PUR}.sorted.bam BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --busco ./BUSCO_Metazoa_Purged/run_metazoa_odb10/full_table.tsv BlobtoolsResults_Purged" >> PurgeDups_QualityAssessment_Pipeline.log

echo -e "\n${DATE}"
echo -e "\nStep 5c: Running the Blobtools analysis on the purged haplotig2 assembly"
module --quiet purge
BlobtoolAnalyses ${HAP2_PUR}

mv BlobtoolsResults_Purged/ BlobtoolsResults_Purged_Haplotig2/
mv BUSCO_Metazoa_Purged/ BUSCO_Metazoa_Purged_Haplotig2/
mv Busco_${HAP2_PUR}.log BUSCO_Metazoa_Purged_Haplotig2/
rm -rf busco_downloads/
rm -rf busco_data
mkdir -p BAM_results_Purged_Haplotig2
mv MappedReads_* BAM_results_Purged_Haplotig2/

echo -e "\n\n********** To display the data in a webserver do the following: **********" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools view --remote ./BlobtoolsResults_Purged_Haplotig2" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "follow the instructions at https://github.com/blobtoolkit/tutorials/blob/main/2021-04-28_blobtoolkit_commandline_mollusc.md for \"Run a viewer for a BlobDir\", start at \"Run the blobtools view command:\"" >> PurgeDups_QualityAssessment_Pipeline.log

#Background reading and help:
#https://blobtoolkit.genomehubs.org/blobtools2/
#https://www.youtube.com/channel/UC6QimXygwjaX2CGZ8cgjNxw
#https://github.com/blobtoolkit/tutorials/blob/main/2021-04-28_blobtoolkit_commandline_mollusc.md

#To display the data in a webserver do the following
#blobtools view --remote ./BlobtoolsResults_Primary_Purged 
#follow the instructions in the last link above for "Run a viewer for a BlobDir", start at "Run the blobtools view command:"

### Cleaning up folder ###
DURA_BLOB_SEC=$SECONDS
DURA_BLOB=$((($DURA_BLOB_SEC-$DURA_BUSCO_SEC)/3600))
DATE=$(date)
echo -e "Duration of the Blobtools analysis on the purged haplotig1 assembly in hours: ${DURA_BLOB}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\nCleaning up the folder" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 6: Cleaning up the folder"
rm -rf *.dmp nodesDB.txt nucl_gb.accession2taxid
rm -rf taxdump/

DURA_ALL=$(($SECONDS/3600))
DATE=$(date)
echo -e "Duration of all analyses in hours: ${DURA_ALL}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\nThe purging and quality assessment pipeline ended at ${DATE}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nAnalyses finished.\n\n"


