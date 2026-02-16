# Genome assembly and annotation of Aegialites sp. (Mannerheim 1853) (Insecta: Coleoptera: Salpingidae). 
This folder contains the scripts used for the analyses conducted as part of the analyes of the genomes of Aegialites. More details can be found in the paper itself.

## Analytical procedures

1. GenomeSizeEstimation_Smudgeplot_v1.1.sh: Bash script used to determine the kmer distribution using Jellyfish of the PacBio reads and ploidy level using Smudgeplot.
2. AssemblyHifiAsm_QualityAssessment_v3.0.sh: Bash script used to assemble the genome using hifiasm, to quality check the first assemblies (primary, haplotig1 and haplotig2) using QUAST, merqury, BUSCO and BlobToolKit and to retrieve the mitochondrial genomes and identify the sequenced species.
3. FilteredDatasets_QualityAssessment_v3.0.sh: Bash script used to quality check the filtered assemblies using QUAST, merqury, BUSCO and BlobToolKit.
4. PurgeDups_QualityAssessment_v3.0.sh: Bash script used to purge duplicates from the filtered assemblies using Purge_Dups and to quality check the filtered and purged assemblies using QUAST, merqury, BUSCO and BlobToolKit.
5. RepeatModeling_v1.0.sh: Bash script used to model repeat elements in the primary assembly.
6. RepeatMasking_v1.1.sh: Bash script used to annotate and mask repeat elements in the primary assembly using simple repeats, Dfam databases and species-specific elements and to generate some summary statistics and graphics.
7. Annotation_v1.sh: Bash script used to structurally and functionally annotate the genome using miniprot, galba, evidencemodeler, Augustus and InterProScan. A quality check is conducted with BUSCO.
