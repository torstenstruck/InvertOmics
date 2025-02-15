#!/bin/sh

#SBATCH --job-name=Ga88_MitoGenomes
#SBATCH --output=Ga88_MitoGenomes_slurm-%j.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=t.h.struck@nhm.uio.no

set -o errexit

# Put your commands here

sh RetrieveSpecificSequences.sh blastn Ga88_allHifireads.fastq_pri_asm.contigs.fasta Ga52_COX1start.fas