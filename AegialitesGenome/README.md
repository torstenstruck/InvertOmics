# Genome assembly and annotation of Aegialites sp. (Mannerheim 1853) (Insecta: Coleoptera: Salpingidae). 
This folder contains the scripts used for the analyses conducted as part of the analyes of the genomes of Aegialites. More details can be found in the paper itself.

## Analytical procedures

1. GenomeSizeEstimation_Smudgeplot_v1.1.sh: Bash script used to determine the kmer distribution using Jellyfish and pliody level using Smudgeplot.
2. AssemblyHifiAsm_QualityAssessment_v3.0.sh: Bash script used to assemble the genome using hifiasm, to quality check the first assemblies (primary, haplotig1 and haplotig2) using QUAST, merqury, BUSCO and BlobToolKit and to retrieve the mitochondrial genomes and identify the sequenced species.
3. FilteredDatasets_QualityAssessment_v3.0.sh: Bash script used to to quality check the filtered assemblies using QUAST, merqury, BUSCO and BlobToolKit.
4. 
