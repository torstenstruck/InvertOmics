# Genome assembly of the Emplectonema genome 
This folder contains the scripts used for the genome assembly without scaffolding for the paper on the reference genome of Emplectonema gracile (Nemertea). More details can be found in the paper itself.

## Analytical procedures

1. The script "GenomeSizeEstimation_Smudgeplot_v1.1.sh" generates kmer-distrubitions using JellyFish and conducts a SmudgePlot analysis.  
2. The script "AssemblyHifiAsm_QualityAssessment_v1.7.sh" conducts a HifiAsm assembly and then conducts several quality checks (i.e., QUAST, Merqury, Blast search for barcoding sequences and species identification, BUSCO, BlobTools).
3. The script "PurgeDups_QualityAssessment_v1.4.sh" purges duplicates from the assembly and then conducts several quality checks (i.e., QUAST, Merqury, BUSCO, BlobTools).
4.   The script "FilteredDatasets_QualityAssessment_v1.3.sh" conducts several quality checks (i.e., QUAST, Merqury, BUSCO, BlobTools) on a provided assembly (i.e., after decontamination and after scaffolding).
