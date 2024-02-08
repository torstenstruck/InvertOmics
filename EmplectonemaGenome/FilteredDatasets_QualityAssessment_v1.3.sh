#!/bin/sh

# (c) Torsten Hugo Struck @29.04.2022 

set -o errexit

### Usage of pipeline ###
# sh FilteredDatasets_QualityAssessment_v1.1.sh [Path and name to fastq file] [Name of assembly file] [Name to meryl database] [e value] [Maximum number of hits] [Maximum number of threads] [Genome size] > FilteredDatasets_QualityAssessment_Pipeline.out
#
# 1) Path and name of fastq file: name of the fastq or fasta to be used for the assembly
# 2) Name of assembly file: name of the fasta file to be used for the purging of duplicates
# 3) Name of meryl database: name of the meryl database (.meryl) to be used for the merqury analyses (assuming that it is in the same folder from which the analysis is started; if this is not the case then provide the path before the name)
# 4) e value: the e value to be used for the different blast searches (e.g., 1e-20)
# 5) Maximum number of hits: the maximum number of hits to be saved in the blast searches (e.g., 1)
# 6) Maximum number of threads: the maximum number of threads that should be used (e.g., 64) IMPORTANT NOTE: It must be an even number!!!
# 7) Genome size: the estimated or known genome size in K, M or G (e.g., 926M)
#
# example of usage: sh FilteredDatasets_QualityAssessment_v1.1.sh ../reads.fastq assembly.fasta reads.fastq.meryl 1e-20 1 64 926M > FilteredDatasets_QualityAssessment_Pipeline.out
#
# IMPORTANT NOTE: The parameters have to be provided in this order and a value for each parameter has to be provided!!!
# IMPORTANT NOTE: If you do changes to the script ALWAYS save the file under a different name. Do NOT modify this master script.
# IMPORTANT NOTE: Recovery points have still to be implemented.

### Version history ###
# v1.1: Changed from Blast+ 2.11.0 to 2.13.0; changed to Blobtools 3.1.4
# v1.2: Commented out HapPy analyses due to re-occuring error messages (it needs to be check thoroughly checked)
# v1.3: added Merqury analysis

### Variables ###
FASTQ=$1
ASS=$2
MERYL=$3
EVALUE=$4
MAXHITS=$5
THREADS=$6
REDTHREADS=$((${THREADS}/2))
GENOMESIZE=$7
SECONDS=0

# print out settings for the pipeline to the log
echo -e "The pipeline was called as follows:" > FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "sh FilteredDatasets_QualityAssessment_v1.0.sh ${FASTQ} ${ASS} ${THREADS} ${GENOMESIZE}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "Name and path of the fastq file: ${FASTQ}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "Name of the fasta file: ${ASS}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "Maximum number of threads: ${THREADS}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "Genome size: ${GENOMESIZE}" >> FilteredDatasets_QualityAssessment_Pipeline.log
DATE=$(date)
echo -e "\nPipeline was started: $DATE" >> FilteredDatasets_QualityAssessment_Pipeline.log 

### Different subscripts used during the script ###
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
blobtools create --fasta ${ASSEMBLY} BlobtoolsResults_${ASSEMBLY}
#Add blast hit information to the BlobtoolsResults folder
blobtools add --hits ${ASSEMBLY}.blobplot.blast.out --taxrule bestsumorder --taxdump ./taxdump BlobtoolsResults_${ASSEMBLY}
#Add mapping coverage to the BlobtoolsResults folder
blobtools add --cov MappedReads_${ASSEMBLY}.sorted.bam BlobtoolsResults_${ASSEMBLY}
#Add BUSCO scores to the BlobtoolsResults folder
blobtools add --busco ./BUSCO_Metazoa_${ASSEMBLY}/run_metazoa_odb10/full_table.tsv BlobtoolsResults_${ASSEMBLY}
}

function HapPyAnalyses () {
# Set the local variables
local ASSEMBLY=$1

module load HapPy/0.2.1-foss-2020a

mkdir -p HaploidyResults_${ASSEMBLY}
cd HaploidyResults_${ASSEMBLY}

#HapPy commands to determine haploidy
happy coverage --outdir=./HapPy_coverage/ ../MappedReads_${ASSEMBLY}.sorted.bam
happy estimate --size $GENOMESIZE --plot --outstats HapPy_estimate_${ASSEMBLY}.txt ./HapPy_coverage/MappedReads_${ASSEMBLY}.sorted.bam.hist
happy autoest --size=$GENOMESIZE --plot --outstats HapPy_autoest_${ASSEMBLY}.txt ./HapPy_coverage/MappedReads_${ASSEMBLY}.sorted.bam.hist

cd ..
}

### Quality assessment using Quast on the filtered data ###
echo -e "\nRunning QUAST analyses on the filtered contigs" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 1: Running QUAST analyses on the filtered data"
module --quiet purge
module load QUAST/5.0.2-foss-2020a-Python-3.8.2 

# run quast on both data
quast.py -o HifiAsm_Quast_AssemblyQuality_${ASS} ${ASS}

### Quality assessment using Merquery on purged primary and alternative assembly###
echo -e "\nRunning Merqury analyses on the filtered contigs" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 2: Running Merqury analyses on the filtered contigs"

module --quiet purge
module load merqury/1.3

#define path to merqury srcipts
export MERQURY=/opt/software/custom/software/merqury-1.3/

#run merqury analysis
sh /opt/software/custom/software/merqury-1.3/merqury.sh ${MERYL} ${ASS} Filtered > Filtered_merqury.out 2> Filtered_merqury.err

### BUSCO check on assembled data ###
DURA_QUAST_SEC=$SECONDS
DURA_QUAST=$(($DURA_QUAST_SEC/3600))
DATE=$(date)
echo -e "Duration of QUAST and Merqury analyses in hours: ${DURA_QUAST}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\nRunning BUSCO analyses on the filtered contigs using the following command line: busco -i ${ASS} -l metazoa_odb10 --out BUSCO_Metazoa_${ASS} -m geno -c ${THREADS} -f > Busco_${ASS}.log" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 3: Running the BUSCO analysis"
module --quiet purge
module load BUSCO/4.1.4-foss-2019b-Python-3.7.4 
export AUGUSTUS_CONFIG_PATH=/storage/ConfigFiles/BRAKER2/config/
export BUSCO_CONFIG_FILE=/storage/ConfigFiles/BUSCO/config.ini

# set the BUSCO settings
busco -i ${ASS} -l metazoa_odb10 --out BUSCO_Metazoa_${ASS} -m geno -c ${THREADS} -f > Busco_${ASS}.log

### Run blobtool analyses ###
### load necessary for modules for the preparation steps ###
DURA_BUSCO_SEC=$SECONDS
DURA_BUSCO=$((($DURA_BUSCO_SEC-$DURA_QUAST_SEC)/3600))
DATE=$(date)
echo -e "Duration of BUSCO analysis in hours: ${DURA_BUSCO}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\nRunning the Blobtools analysis on the filtered contigs using the following command lines:" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "blastn -task megablast -query ${ASS} -db nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs ${MAXHITS} -max_hsps 1 -num_threads ${REDTHREADS} -evalue ${EVALUE} -out ${ASS}.blobplot.blast.out" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "blobtools create --fasta ${ASS} BlobtoolsResults_${ASS}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "blobtools add --hits ${ASS}.blobplot.blast.out --taxrule bestsumorder --taxdump ./taxdump BlobtoolsResults_${ASS}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "blobtools add --cov MappedReads_${ASS}.sorted.bam BlobtoolsResults_${ASS}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "blobtools add --busco ./BUSCO_Metazoa_${ASS}/run_metazoa_odb10/full_table.tsv BlobtoolsResults_${ASS}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "To display the data in a webserver do the following:" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "blobtools view --remote ./BlobtoolsResults_${ASS}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "follow the instructions at https://github.com/blobtoolkit/tutorials/blob/main/2021-04-28_blobtoolkit_commandline_mollusc.md for \"Run a viewer for a BlobDir\", start at \"Run the blobtools view command:\"" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 4: Running the Blobtools analysis"

#Fetch the NCBI Taxdump and blast-ready database
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
perl /storage/software/easybuild/software/BLAST+/2.13.0-gompi-2020a/bin/update_blastdb.pl --passive --decompress nt

module --quiet purge
BlobtoolAnalyses ${ASS}

#Background reading and help:
#https://blobtoolkit.genomehubs.org/blobtools2/
#https://www.youtube.com/channel/UC6QimXygwjaX2CGZ8cgjNxw
#https://github.com/blobtoolkit/tutorials/blob/main/2021-04-28_blobtoolkit_commandline_mollusc.md

#To display the data in a webserver do the following
#blobtools view --remote ./BlobtoolsResults_${ASM_UNPUR} or ./BlobtoolsResults_${ASM_PUR} 
#follow the instructions in the last link above for "Run a viewer for a BlobDir", start at "Run the blobtools view command:"

### Run HapPy analyses ###
DURA_BLOB_SEC=$SECONDS
DURA_BLOB=$((($DURA_BLOB_SEC-$DURA_BUSCO_SEC)/3600))
DATE=$(date)
echo -e "Duration of the BlobTools analysis in hours: ${DURA_BLOB_SEC}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\nRunning the HapPy analysis on the purged contigs using the following command lines:" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "happy coverage --outdir=./HapPy_coverage/ MappedReads_${ASS}.sorted.bam" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "happy estimate --size ${GENOMESIZE} --plot --outstats HapPy_estimate_${ASS}.txt ./HapPy_coverage/MappedReads_${ASS}.sorted.bam.hist" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "happy autoest --size=${GENOMESIZE} --plot --outstats HapPy_autoest_${ASS}.txt ./HapPy_coverage/MappedReads_${ASS}.sorted.bam.hist" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 5: Running the HapPy analysis; UPDATE: Due to problems with module HapPy analyses are commented out at the moment"

module --quiet purge
#HapPyAnalyses ${ASS} 

### Cleaning up folder ###
DURA_HAPPY_SEC=$SECONDS
DURA_HAPPY=$((($DURA_HAPPY_SEC-$DURA_BUSCO_SEC)/3600))
DATE=$(date)
echo -e "Duration of the HapPy analysis in hours: ${DURA_HAPPY}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\nCleaning up the data and moving results into appropriate folders" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nStep 6: Cleaning up the data and moving results"
rm -rf busco_downloads/
rm -rf busco_data

### Sorting results ###
mkdir -p BAM_results__${ASS}
mv MappedReads_* BAM_results__${ASS}/
mv Busco_${ASS}.log BUSCO_Metazoa_${ASS}/

DURA_ALL=$(($SECONDS/3600))
DATE=$(date)
echo -e "Duration of all analyses in hours: ${DURA_ALL}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\nThe assembly and quality assessment pipeline ended at ${DATE}" >> FilteredDatasets_QualityAssessment_Pipeline.log
echo -e "\n${DATE}"
echo -e "\nAnalyses finished.\n\n"
