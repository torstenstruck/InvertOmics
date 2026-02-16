# Genome assembly and annotation of Aegialites sp. (Mannerheim 1853) (Insecta: Coleoptera: Salpingidae). 
This folder contains the scripts used for the analyses conducted as part of the analyes of the genomes of Aegialites. More details can be found in the paper itself.

## Analytical procedures

1. Download_GenomesNCBI.sh: Bash script used to download public genomes from NCBI at at least the chromosome-level, including checks for redundant versions and if the genomes had alread been downloaded before.
2. RetrieveInformationPublicGenomes_v1.sh: Bash script to retrieve different parameters from the meta-data associated with the downloaded public genomes.
3. DataExplorationCompiled.R: R script for statistical analyses of the parameters of public chromosome-level genomes.
4. DataExploration_Scaffolds.R: R script for statistical analyses of the parameters of public scaffold-level genomes.
5. GenomeSizeEstimation_Smudgeplot_v1.1.sh: Bash script to conduct jellyfish and smudgeplot analyses.
6. CountReadlength.sh: Bash script to count the read length of each sequence read.
7. CountReadlength_v2.r: R script to statistically analyse the read length.
8. AssemblyHifiAsm_QualityAssessment_v2.1.sh: Bash script for the first assembly and quality check of the assembly.
9. FilteredDatasets_QualityAssessment_v2.1.sh: Bash script for the quality check of the assembly after filtering with Blobtools.
10. PurgeDups_QualityAssessment_v2.1.sh: Bash script for the purging of duplicated contigs from the filtered assembly and quality check of the filtered and purged assembly.
11. CorrelationAnalysesAssemblies.R: R script for statistical analyses of the parameters of preliminary genome assemblies of the InvertOmics project.
