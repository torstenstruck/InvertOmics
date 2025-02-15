#!/bin/sh

#SBATCH --job-name=Ga82_HifiAsmMito
#SBATCH --output=Ga82_HifiAsmMito_slurm-%j.txt
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=t.h.struck@nhm.uio.no

set -o errexit

# Put your commands here
module load hifiasm/0.18.2-GCCcore-10.3.0

hifiasm --primary -t 32 -f 0 --min-hist-cnt 2 -o Ga82_allHifireads2.mitohits_cleaned_assembled Ga82_allHifireads2.mitohits_cleaned.fas

awk '/^S/{print ">"$2;print $3}' Ga82_allHifireads2.mitohits_cleaned_assembled.p_ctg.gfa > Ga82_allHifireads2.mitohits_cleaned_assembled.contigs.fasta
awk '/^S/{print ">"$2;print $3}' Ga82_allHifireads2.mitohits_cleaned_assembled.a_ctg.gfa >> Ga82_allHifireads2.mitohits_cleaned_assembled.contigs.fasta
