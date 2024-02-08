#!/bin/sh

# (c) Torsten Hugo Struck @04.05.2022 

set -o errexit

### Usage of pipeline ###
# sh PurgeDups_QualityAssessment_v1.4.sh [Name of fastq file] [Name of primary assembly file] [Name of haplotig1 file] [Name of haplotig2 file] [Name to meryl database] [e value] [Maximum number of hits] [Maximum number of threads] [Genome size] >> PurgeDups_QualityAssessment_Pipeline.out
#
# 1) Name of fastq file: name of the fastq or fasta to be used for the assembly (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 2) Name of primary assembly file: name of the fasta file to be used for the purging of duplicates from the primary assembly (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 3) Name of haplotig1 file: name of the fasta file to be used for the purging of duplicates from haplotig1 (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 4) Name of haplotig2 file: name of the fasta file to be used for the purging of duplicates from haplotig2 (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 5) Name of meryl database: name of the meryl database (.meryl) to be used for the merqury analyses (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 6) e value: the e value to be used for the different blast searches (e.g., 1e-20)
# 7) Maximum number of hits: the maximum number of hits to be saved in the blast searches (e.g., 1)
# 8) Maximum number of threads: the maximum number of threads that should be used (e.g., 64) IMPORTANT NOTE: It must be an even number!!!
# 9) Genome size: the estimated or known genome size in K, M or G (e.g., 926M)
#
# example of usage: sh PurgeDups_QualityAssessment_v1.4.sh ../reads.fastq pri_assembly.fasta hap1_assembly.fasta hap2_assembly.fasta reads.fastq.meryl 1e-20 1 64 926M >> PurgeDups_QualityAssessment_Pipeline.out
#
# IMPORTANT NOTE: The parameters have to be provided in this order and a value for each parameter has to be provided!!!
# IMPORTANT NOTE: If you do changes to the script ALWAYS save the file under a different name. Do NOT modify this master script.
# IMPORTANT NOTE: Recovery points have still to be implemented.

### Version history ###
# v1.1: Changed from Blast+ 2.11.0 to 2.13.0
# v1.2: Changed to Blobtools 3.1.4; fixed an issue with the wait commands
# v1.3: Added Purge_dups on haplotigs; added Merqury analysis
# v1.4: Commented out HapPy analyses due to re-occuring error messages (it needs to be check thoroughly checked); fixed a bug concerning BUSCO scores and Blobtool analyses

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
GENOMESIZE=$9
SECONDS=0

# print out settings for the pipeline to the log
echo -e "The pipeline was called as follows:" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "sh PurgeDups_QualityAssessment_v1.3.sh ${FASTQ} ${ASS} ${HAP1} ${HAP2} ${MERYL} ${EVALUE} ${MAXHITS} ${THREADS} ${GENOMESIZE}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Name of the fastq file: ${FASTQ}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Name of the fasta file of the primary assembly: ${ASS}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Name of the fasta file of the haplotig1: ${HAP1}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Name of the fasta file of the haplotig2: ${HAP2}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Name of the meryl database: ${MERYL}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "e value of blast search: ${EVALUE}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Maximum number of hits in blast searcg: ${MAXHITS}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Maximum number of threads: ${THREADS}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "Genome size: ${GENOMESIZE}" >> PurgeDups_QualityAssessment_Pipeline.log
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
echo -e "Mapping ended at time $(date +%T).\n"
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
echo -e "Blast search ended at time $(date +%T).\n"
wait $PID2

### Run actual blobtool analyses ###
module load Blobtools/3.1.4

#Create the BlobtoolsResults folder containing the assembly
blobtools create --fasta ${ASSEMBLY} BlobtoolsResults_Primary_Purged
#Add blast hit information to the BlobtoolsResults folder
blobtools add --hits ${ASSEMBLY}.blobplot.blast.out --taxrule bestsumorder --taxdump ./taxdump BlobtoolsResults_Primary_Purged
#Add mapping coverage to the BlobtoolsResults folder
blobtools add --cov MappedReads_${ASSEMBLY}.sorted.bam BlobtoolsResults_Primary_Purged
#Add BUSCO scores to the BlobtoolsResults folder
blobtools add --busco ./BUSCO_Metazoa_Primary_Purged/run_metazoa_odb10/full_table.tsv BlobtoolsResults_Primary_Purged
}

function HapPyAnalyses () {
# Set the local variables
local ASSEMBLY=$1

module load HapPy/0.2.1-foss-2020a

mkdir -p HaploidyResults_Primary_Purged
cd HaploidyResults_Primary_Purged

#HapPy commands to determine haploidy
happy coverage --outdir=./HapPy_coverage/ ../MappedReads_${ASSEMBLY}.sorted.bam
happy estimate --size $GENOMESIZE --plot --outstats HapPy_estimate_${ASSEMBLY}.txt ./HapPy_coverage/MappedReads_${ASSEMBLY}.sorted.bam.hist
happy autoest --size=$GENOMESIZE --plot --outstats HapPy_autoest_${ASSEMBLY}.txt ./HapPy_coverage/MappedReads_${ASSEMBLY}.sorted.bam.hist

cd ..
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

### Quality assessment using Quast on purged primary data ###
DURA_PURGE_SEC=$SECONDS
DURA_PURGE=$(($DURA_PURGE_SEC/3600))
DATE=$(date)
echo -e "Duration of purging duplicates in hours: ${DURA_PURGE}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\nRunning QUAST analyses on the purged files" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 2: Running QUAST analyses on the purged data"
module --quiet purge
module load QUAST/5.0.2-foss-2020a-Python-3.8.2 

# run quast on both data
quast.py -o Quast_Primary_Purged ${ASM_PUR}
quast.py -o Quast_Haplotig1_Purged ${HAP1_PUR}
quast.py -o Quast_Haplotig2_Purged ${HAP2_PUR}

### Quality assessment using Merquery on purged primary and alternative assembly###
echo -e "\nRunning Merqury analyses on the purged primary assembly and haplotigs" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 3: Running Merqury analyses on the purged primary and alternative assembly"

module --quiet purge
module load merqury/1.3

#define path to merqury srcipts
export MERQURY=/opt/software/custom/software/merqury-1.3/

#run merqury analysis
sh /opt/software/custom/software/merqury-1.3/merqury.sh ${MERYL} ${ASM_PUR} PurgedPri > PurgedPri_merqury.out 2> PurgedPri_merqury.err
sh /opt/software/custom/software/merqury-1.3/merqury.sh ${MERYL} ${HAP1_PUR} ${HAP2_PUR} PurgedHaps > PurgedHaps_merqury.out 2> PurgedHaps_merqury.err

### BUSCO check on assembled data ###
DURA_QUAST_SEC=$SECONDS
DURA_QUAST=$((($DURA_QUAST_SEC-$DURA_PURGE_SEC)/3600))
DATE=$(date)
echo -e "Duration of QUAST and Merqury analyses in hours: ${DURA_QUAST}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 3: Running the BUSCO analysis"
module --quiet purge
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4 
export AUGUSTUS_CONFIG_PATH=/storage/ConfigFiles/BRAKER2/config/
export BUSCO_CONFIG_FILE=/storage/ConfigFiles/BUSCO/config.ini

# set the BUSCO settings
echo -e "\nRunning BUSCO analyses on the purged primary assembly using the following command line: busco -i ${ASM_PUR} -l metazoa_odb10 --out BUSCO_Metazoa_${ASM_PUR} -m geno -c ${THREADS} -f > Busco_${ASM_PUR}.log" >> PurgeDups_QualityAssessment_Pipeline.log
busco -i ${ASM_PUR} -l metazoa_odb10 --out BUSCO_Metazoa_Primary_Purged -m geno -c ${THREADS} -f > Busco_${ASM_PUR}.log
echo -e "\nRunning BUSCO analyses on the purged haplotig1 using the following command line: busco -i ${HAP1_PUR} -l metazoa_odb10 --out BUSCO_Metazoa_Haplotig1_Purged -m geno -c ${THREADS} -f > Busco_${HAP1_PUR}.log" >> PurgeDups_QualityAssessment_Pipeline.log
busco -i ${HAP1_PUR} -l metazoa_odb10 --out BUSCO_Metazoa_Haplotig1_Purged -m geno -c ${THREADS} -f > Busco_${HAP1_PUR}.log

### Run blobtool analyses ###
### load necessary for modules for the preparation steps ###
DURA_BUSCO_SEC=$SECONDS
DURA_BUSCO=$((($DURA_BUSCO_SEC-$DURA_QUAST_SEC)/3600))
DATE=$(date)
echo -e "Duration of BUSCO analysis in hours: ${DURA_BUSCO}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\nRunning the Blobtools analysis on the purged primary assembly using the following command lines:" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blastn -task megablast -query ${ASM_PUR} -db nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs ${MAXHITS} -max_hsps 1 -num_threads ${REDTHREADS} -evalue ${EVALUE} -out ${ASM_PUR}.blobplot.blast.out" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools create --fasta ${ASM_PUR} BlobtoolsResults_${ASM_PUR}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --hits ${ASM_PUR}.blobplot.blast.out --taxrule bestsumorder --taxdump ./taxdump BlobtoolsResults_${ASM_PUR}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --cov MappedReads_${ASM_PUR}.sorted.bam BlobtoolsResults_${ASM_PUR}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools add --busco ./BUSCO_Metazoa_${ASM_PUR}/run_metazoa_odb10/full_table.tsv BlobtoolsResults_${ASM_PUR}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "To display the data in a webserver do the following:" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "blobtools view --remote ./BlobtoolsResults_${ASM_PUR}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "follow the instructions at https://github.com/blobtoolkit/tutorials/blob/main/2021-04-28_blobtoolkit_commandline_mollusc.md for \"Run a viewer for a BlobDir\", start at \"Run the blobtools view command:\"" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 4: Running the Blobtools analysis"

#Fetch the NCBI Taxdump and blast-ready database
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
perl /storage/software/easybuild/software/BLAST+/2.13.0-gompi-2020a/bin/update_blastdb.pl --passive --decompress nt

module --quiet purge
BlobtoolAnalyses ${ASM_PUR}

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
echo -e "Duration of the Blobtools analysis in hours: ${DURA_BLOB}" >> PurgeDups_QualityAssessment_Pipeline.log
#echo -e "\nRunning the HapPy analysis on the purged contigs using the following command lines:" >> PurgeDups_QualityAssessment_Pipeline.log
#echo -e "happy coverage --outdir=./HapPy_coverage/ MappedReads_${ASM_PUR}.sorted.bam" >> PurgeDups_QualityAssessment_Pipeline.log
#echo -e "happy estimate --size ${GENOMESIZE} --plot --outstats HapPy_estimate_${ASM_PUR}.txt ./HapPy_coverage/MappedReads_${ASM_PUR}.sorted.bam.hist" >> PurgeDups_QualityAssessment_Pipeline.log
#echo -e "happy autoest --size=${GENOMESIZE} --plot --outstats HapPy_autoest_${ASM_PUR}.txt ./HapPy_coverage/MappedReads_${ASM_PUR}.sorted.bam.hist" >> PurgeDups_QualityAssessment_Pipeline.log
#echo -e "\n${DATE}"
#echo -e "\nStep 5: Running the HapPy analysis"

#module --quiet purge
#HapPyAnalyses ${ASM_PUR} 

### Cleaning up folder ###
#DURA_HAPPY_SEC=$SECONDS
#DURA_HAPPY=$((($DURA_HAPPY_SEC-$DURA_BLOB_SEC)/3600))
#DATE=$(date)
#echo -e "Duration of the HapPy analysis in hours: ${DURA_HAPPY}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\nCleaning up the data and moving results into appropriate folders" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 6: Cleaning up the data and moving results"
rm -rf nt.*
rm -rf *.dmp nodesDB.txt nucl_gb.accession2taxid
rm -rf busco_downloads/
rm -rf busco_data
rm -rf taxdump/
rm -rf taxdb.*

### Sorting results ###
mkdir -p BAM_results_Purged
mv MappedReads_* BAM_results_Purged/
mv Busco_${ASM_PUR}.log BUSCO_Metazoa_Primary_Purged/
mv Busco_${HAP1_PUR}.log BUSCO_Metazoa_Haplotig1_Purged/
mkdir Merquery_results_Purged
mv PurgedPri* PurgedHaps* Merquery_results_Purged/

DURA_ALL=$(($SECONDS/3600))
DATE=$(date)
echo -e "Duration of all analyses in hours: ${DURA_ALL}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\nThe assembly and quality assessment pipeline ended at ${DATE}" >> PurgeDups_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nAnalyses finished.\n\n"
