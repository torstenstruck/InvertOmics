#!/bin/sh

#SBATCH --job-name=Ga52_MitoHifi_r
#SBATCH --output=Ga52_MitoHifi_r_slurm-%j.txt
#SBATCH --time=14-00:00:00
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=1T
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=t.h.struck@nhm.uio.no

set -o errexit

# Put your commands here

singularity exec -H $PWD:/home -B $PWD:/data/ /storage/EBPNor_pipeline/MitoHifi/MitoHifi.sif mitohifi.py -r /data/Ga52allHifireads.fastq  -f /data/Ga52_COX1start.fas -g /data/Ga52_COX1start_annotation.gb -p 99 -a animal -t 32 -o 5

## -r is reads.fasta
## -f is reference.fasta
## -g is reference.genbank
## -p is percentage
## -a is animal, fungi or plant
## -t is number of threads for hifiasm, blast, minimap2, samtools etc
## -o is organism code from NCBI table, which is 2 for vertebrates, and 1 for standard. Cannot tell if there is one for plants
