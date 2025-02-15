#!/bin/sh

#SBATCH --job-name=Ga52_PurgePalin_MitoGenome
#SBATCH --output=Ga52_PurgePalin_MitoGenome_slurm-%j.txt
#SBATCH --time=7-0:00:00
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=t.h.struck@nhm.uio.no

set -o errexit

# Put your commands here
eval "$(/storage/EBPNor_pipeline/miniconda3/bin/conda shell.bash hook)"
conda activate purgepalindrom
python /storage/EBPNor_pipeline/opt/Pacasus/pacasus.py --filetype1 fasta -o Ga52allHifireads.mitohits_cleaned.fas --device_type=CPU --platform_name=Intel --framework=opencl Ga52allHifireads.mitohits.fas > PurgingPalindromeSequences.out
conda deactivate