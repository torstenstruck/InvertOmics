#!/bin/sh

# (c) Torsten Hugo Struck @14.11.2021, modified @22.04.2022 

set -o errexit

### Usage of pipeline ###
# sh GenomeSizeEstimation_Smudgeplot_v1.1.sh [Name of fastq file] [KMER SIZE] [THREADS] [MEMORY] > GenomeSizeEstimation_Smudgeplot.out
#
# 1) Name of fastq file: name of the fastq or fasta with reads to be assessed 
# 2) kmer size used for the analysis
# 3) Number of threads used for the analysis
# 4) Memory to be allocated in Gb
# 
# example of usage: sh GenomeSizeEstimation_Smudgeplot_v1.1.sh seqdata.fastq 21 12 5G > GenomeSizeEstimation_Smudgeplot.out
#
# IMPORTANT NOTE: The parameters have to be provided in this order and a value for each parameter has to be provided!!!
# IMPORTANT NOTE: If you do changes to the script ALWAYS save the file under a different name. Do NOT modify this master script.

### Version history ###
# v1.1: adding a higher bound for jellyfish (1 million) to included ultrahigh frequency kmers for counterbalancing that genomescope2 might otherwise underrepresent the repetitiveness of the genome (see last FAQ in https://github.com/tbenavi1/genomescope2.0)

### Variables ###
FASTQ=$1
KMER=$2
THREADS=$3
MEMORY=$4

echo -e "The genoem size estemation was called as follows:" > GenomeSizeEstimation_Smudgeplot.log
echo -e "sh GenomeSizeEstimation_Smudgeplot_v1.1.sh ${FASTQ} ${KMER} ${THREADS} ${MEMORY} " >> GenomeSizeEstimation_Smudgeplot.log

#Generate the kmer distribution by counting the kmers and generating a histogram to be used in GenomeScope2 and Smudgeplot
module --quiet purge  # Reset the modules to the system default
module load Jellyfish/2.3.0-GCC-8.3.0
jellyfish count -t $THREADS -C -m $KMER -s $MEMORY -o ${FASTQ}_${KMER}mer_out $FASTQ
jellyfish histo -o ${FASTQ}_${KMER}mer_out.histo ${FASTQ}_${KMER}mer_out
jellyfish histo -h 1000000 -o ${FASTQ}_${KMER}mer_out.high_1M.histo ${FASTQ}_${KMER}mer_out

#Estimate the genome at http://genomescope.org/genomescope2.0 using the histo-file; this can be done while Smudgeplot is running

#Smudgeplot commands
module load smudgeplot/0.2.4
 
#Extract genomic kmers using reasonable coverage thresholds
L=$(smudgeplot.py cutoff ${FASTQ}_${KMER}mer_out.histo L)
U=$(smudgeplot.py cutoff ${FASTQ}_${KMER}mer_out.histo U)
echo $L $U > Borders_${KMER}mer_smudgeplot.txt # these need to be sane values like 30 800 or so

#Extract kmers in the coverage range from L to U using Jellyfish dump and pipe them directly to smudgeplot.py hetkmers to compute the set of kmer pairs.
jellyfish dump -c -L $L -U $U ${FASTQ}_${KMER}mer_out | smudgeplot.py hetkmers -o ${FASTQ}_${KMER}mer_smudgeplot_pairs

#Generate the smudgeplot using the coverages of the kmer pairs
smudgeplot.py plot ${FASTQ}_${KMER}mer_smudgeplot_pairs_coverages.tsv -o ${FASTQ}_${KMER}mer_smudgeplot_pairs_plot

### Cleaning up folder ###
mkdir -p JellyfishResults
mv ${FASTQ}_${KMER}mer_out.histo ${FASTQ}_${KMER}mer_out.high_1M.histo ${FASTQ}_${KMER}mer_out JellyfishResults

mkdir -p Smudgeplot_Results
mv *_${KMER}mer_smudgeplot* Smudgeplot_Results/ 

