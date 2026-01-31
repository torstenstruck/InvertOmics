#!/bin/bash

# (c) Torsten Hugo Struck @31.07.2025, with the assistance of ChatGPT 

# This script will retrieve information from the json-files associated the public genomes of NCBI. The specific information (if available) gathered are:
# accession 
# buscoLineage
# BUSCO_percentage of complete genes
# BUSCO_percentage of duplicated genes
# BUSCO_percentage of missing genes
# ANNOTATION_number of non-coding genes
# ANNOTATION_number of protein-coding genes genes
# ANNOTATION_number of pseudogenes
# ANNOTATION_total number of genes
# ASSEMBLY_assembly level
# ASSEMBLY_contig L50
# ASSEMBLY_contig N50
# ASSEMBLY_gc percentage
# ASSEMBLY_genome coverage
# ASSEMBLY_number of contigs
# ASSEMBLY_number of scaffolds
# ASSEMBLY_scaffold L50
# ASSEMBLY_scaffold N50
# ASSEMBLY_total number of chromosomes
# ASSEMBLY_total sequence length
# CHROMOSOME_Maximal GC%
# CHROMOSOME_Minimal GC%
# CHROMOSOME_Mean GC%
# CHROMOSOME_Standard deviation GC%
# CHROMOSOME_Median GC%
# CHROMOSOME_Maximum length
# CHROMOSOME_Minimum length
# CHROMOSOME_Mean length
# CHROMOSOME_Standard deviation length
# CHROMOSOME_Median length 

# Exit on error, undefined variable check (treat as error), pipeline failure detection
set -euo pipefail

### Usage of pipeline ###
# sh RetrieveInformationPublicGenomes_v1.sh [Folder_name]
#
# 1) Folder_name: Name of folder that contains the files assembly_data_report.jsonl and $subfolder/$GCA_GCF_FOLDER_NAME/sequence_report.jsonl 
#                 For example, Bivalvia is a possible folder name. At the next level, it needs to have the folder with species and accession number (e.g., Mytilus_edulis__GCF_963676685.1). Within it has to have the structure as it was downloaded from NCBI.
#
# example of usage: sh RetrieveInformationPublicGenomes_v1.sh Bivalvia

# Check if a base directory was provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <base_directory>"
    exit 1
fi

# Define the base directory from the command-line argument
BASE_DIR="$1"

# Check if the provided base directory exists
if [ ! -d "$BASE_DIR" ]; then
    echo "Directory not found: $BASE_DIR"
    exit 1
fi

# Output file for the compiled results
OUTPUT_FILE="$BASE_DIR/${BASE_DIR}_assembly_data_summary.tsv"

# Check if the output file already exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "Error: Output file already exists: $OUTPUT_FILE"
    exit 1
fi

# Write the header to the output file
echo -e "taxon\tspecies_accession\taccession\tBUSCO_buscoLineage\tBUSCO_complete\tBUSCO_duplicated\tBUSCO_missing\tANNOTATION_nonCoding\tANNOTATION_proteinCoding\tANNOTATION_pseudogene\tANNOTATION_total\tASSEMBLY_assemblyLevel\tASSEMBLY_contigL50\tASSEMBLY_contigN50\tASSEMBLY_gcPercent\tASSEMBLY_genomeCoverage\tASSEMBLY_numberOfContigs\tASSEMBLY_numberOfScaffolds\tASSEMBLY_scaffoldL50\tASSEMBLY_scaffoldN50\tASSEMBLY_totalNumberOfChromosomes\tASSEMBLY_totalSequenceLength\tCHR_Max_GC%\tCHR_Min_GC%\tCHR_Mean_GC%\tCHR_StDev_GC%\tCHR_Median_GC%\tCHR_MaxLength\tCHR_MinLength\tCHR_MeanLength\tCHR_StDevLength\tCHR_MedianLength" > "$OUTPUT_FILE"

# Function to calculate median
calculate_median() {
    local -n array=$1
    local count=${#array[@]}
    if [ $count -eq 0 ]; then
        echo "-"
        return
    fi

    # Sort the array
    sorted=($(printf "%s\n" "${array[@]}" | sort -n))
    if (( count % 2 == 0 )); then
        mid=$(( count / 2 ))
        median=$(echo "scale=2; (${sorted[mid - 1]} + ${sorted[mid]}) / 2" | bc)
    else
        mid=$(( (count + 1) / 2 ))
        median=${sorted[mid - 1]}
    fi
    echo "$median"
}

# Function to calculate standard deviation
calculate_stddev() {
    local -n array=$1
    local count=${#array[@]}
    if [ $count -eq 0 ]; then
        echo "-"
        return
    fi

    # Calculate the mean
    sum=0
    for val in "${array[@]}"; do
        sum=$(echo "$sum + $val" | bc)
    done
    mean=$(echo "scale=2; $sum / $count" | bc)

    # Calculate the variance
    variance_sum=0
    for val in "${array[@]}"; do
        variance_sum=$(echo "$variance_sum + (($val - $mean)^2)" | bc)
    done
    variance=$(echo "scale=2; $variance_sum / $count" | bc)
    stddev=$(echo "scale=2; sqrt($variance)" | bc)
    echo "$stddev"
}

# Iterate over each folder inside the base directory
for subfolder in "$BASE_DIR"/*; do
    if [ -d "$subfolder" ]; then
        # Extract the folder name for output
        FOLDER_NAME=$(basename "$subfolder")
#        echo $FOLDER_NAME
		
        # Construct the paths to the JSONL files
        ASSEMBLY_DATA_FILE="$subfolder/assembly_data_report.jsonl"
        echo $ASSEMBLY_DATA_FILE

         if [ -f "$ASSEMBLY_DATA_FILE" ]; then
            # Use jq to extract the required fields
            ACCESSION=$(jq -r '.accession' "$ASSEMBLY_DATA_FILE")
            BUSCO_LINEAGE=$(jq -r '.annotationInfo.busco.buscoLineage' "$ASSEMBLY_DATA_FILE")
            BUSCO_COMPLETE=$(jq -r '.annotationInfo.busco.complete' "$ASSEMBLY_DATA_FILE")
            BUSCO_DUPLICATED=$(jq -r '.annotationInfo.busco.duplicated' "$ASSEMBLY_DATA_FILE")
            BUSCO_MISSING=$(jq -r '.annotationInfo.busco.missing' "$ASSEMBLY_DATA_FILE")
            NON_CODING=$(jq -r '.annotationInfo.stats.geneCounts.nonCoding' "$ASSEMBLY_DATA_FILE")
            PROTEIN_CODING=$(jq -r '.annotationInfo.stats.geneCounts.proteinCoding' "$ASSEMBLY_DATA_FILE")
            PSEUDOGENE=$(jq -r '.annotationInfo.stats.geneCounts.pseudogene' "$ASSEMBLY_DATA_FILE")
            TOTAL=$(jq -r '.annotationInfo.stats.geneCounts.total' "$ASSEMBLY_DATA_FILE")
            ASSEMBLY_LEVEL=$(jq -r '.assemblyInfo.assemblyLevel' "$ASSEMBLY_DATA_FILE")
            CONTIG_L50=$(jq -r '.assemblyStats.contigL50' "$ASSEMBLY_DATA_FILE")
            CONTIG_N50=$(jq -r '.assemblyStats.contigN50' "$ASSEMBLY_DATA_FILE")
            GC_PERCENT=$(jq -r '.assemblyStats.gcPercent' "$ASSEMBLY_DATA_FILE")
            GENOME_COVERAGE=$(jq -r '.assemblyStats.genomeCoverage' "$ASSEMBLY_DATA_FILE")
            NUMBER_OF_CONTIGS=$(jq -r '.assemblyStats.numberOfContigs' "$ASSEMBLY_DATA_FILE")
            NUMBER_OF_SCAFFOLDS=$(jq -r '.assemblyStats.numberOfScaffolds' "$ASSEMBLY_DATA_FILE")
            SCAFFOLD_L50=$(jq -r '.assemblyStats.scaffoldL50' "$ASSEMBLY_DATA_FILE")
            SCAFFOLD_N50=$(jq -r '.assemblyStats.scaffoldN50' "$ASSEMBLY_DATA_FILE")
            TOTAL_NUMBER_OF_CHROMOSOMES=$(jq -r '.assemblyStats.totalNumberOfChromosomes' "$ASSEMBLY_DATA_FILE")
            TOTAL_SEQUENCE_LENGTH=$(jq -r '.assemblyStats.totalSequenceLength' "$ASSEMBLY_DATA_FILE")

        else
            echo "ERROR: File not found: $ASSEMBLY_DATA_FILE"
            exit 1
       fi

       # Initialize a variable for the GCA/GCF folder
        GCA_GCF_FOLDER_NAME=""

        # Loop through the directories in the current subfolder to find the GCA/GCF folder
        for gca_gcf_folder in "$subfolder"/*; do
            if [ -d "$gca_gcf_folder" ] && [[ $(basename "$gca_gcf_folder") == GCA_* || $(basename "$gca_gcf_folder") == GCF_* ]]; then
                GCA_GCF_FOLDER_NAME=$(basename "$gca_gcf_folder")
                break  # Stop after finding the first matching GCA/GCF folder
            fi
        done

        # Construct the path to the sequence report file if a GCA/GCF folder is found
        SEQUENCE_REPORT_FILE="$subfolder/$GCA_GCF_FOLDER_NAME/sequence_report.jsonl"
#		echo $GCA_GCF_FOLDER_NAME
#		echo $SEQUENCE_REPORT_FILE

        if [[ -n $GCA_GCF_FOLDER_NAME && -f "$SEQUENCE_REPORT_FILE" ]]; then
            echo "Processing sequence report file: $SEQUENCE_REPORT_FILE"

            # Initialize arrays for statistics
            declare -a gc_percent_values=()
            declare -a length_values=()

            # Read sequence report and gather required data
            while IFS= read -r line; do
                chr_name=$(echo "$line" | jq -r '.chrName')
                gc_count=$(echo "$line" | jq -r '.gcPercent')
                length=$(echo "$line" | jq -r '.length')
                role=$(echo "$line" | jq -r '.role')

                # Store values based on chrName
                if [[ "$role" == "assembled-molecule" && "$chr_name" != "OsHV1_1" && ! "$chr_name" =~ ^[mM] ]]; then
                    gc_percent_values+=("$gc_count")
                    length_values+=("$length")
                else
                    # Handle non-chromosome values separately (if needed)
                    continue
                fi
            done < "$SEQUENCE_REPORT_FILE"

            # Calculate statistics for gc_percent
            if [ ${#gc_percent_values[@]} -gt 0 ]; then
                max_gc=$(printf "%s\n" "${gc_percent_values[@]}" | sort -n | tail -1)
                min_gc=$(printf "%s\n" "${gc_percent_values[@]}" | sort -n | head -1)

                # Calculate mean
                sum_gc=0
                for val in "${gc_percent_values[@]}"; do
                    sum_gc=$(echo "$sum_gc + $val" | bc)
                done
                mean_gc=$(echo "scale=2; $sum_gc / ${#gc_percent_values[@]}" | bc)
				median_gc=$(calculate_median gc_percent_values)
                median_gc=$(echo "scale=2; $median_gc / 1" | bc)                
                stddev_gc=$(calculate_stddev gc_percent_values)            
                stddev_gc=$(echo "scale=2; $stddev_gc / 1" | bc)                
				else
                max_gc="-"; min_gc="-"; mean_gc="-"; median_gc="-"; stddev_gc="-"
            fi

            # Calculate statistics for length
            if [ ${#length_values[@]} -gt 0 ]; then
                max_length=$(printf "%s\n" "${length_values[@]}" | sort -n | tail -1)
                min_length=$(printf "%s\n" "${length_values[@]}" | sort -n | head -1)

                # Calculate mean
                sum_length=0
                for val in "${length_values[@]}"; do
                    sum_length=$(echo "$sum_length + $val" | bc)
                done
                mean_length=$(echo "scale=0; $sum_length / ${#length_values[@]}" | bc)
                median_length=$(calculate_median length_values)
                median_length=$(echo "scale=0; $median_length / 1" | bc)                
				stddev_length=$(calculate_stddev length_values)
                stddev_length=$(echo "scale=0; $stddev_length / 1" | bc)                
            else
                max_length="-"; min_length="-"; mean_length="-"; median_length="-"; stddev_length="-"
            fi

        else
            echo "ERROR: File not found: $SEQUENCE_REPORT_FILE"
            exit 1
        fi
		
        # Compile the results and write to the output TSV file
        echo -e "$BASE_DIR\t$FOLDER_NAME\t$ACCESSION\t$BUSCO_LINEAGE\t$BUSCO_COMPLETE\t$BUSCO_DUPLICATED\t$BUSCO_MISSING\t$NON_CODING\t$PROTEIN_CODING\t$PSEUDOGENE\t$TOTAL\t$ASSEMBLY_LEVEL\t$CONTIG_L50\t$CONTIG_N50\t$GC_PERCENT\t$GENOME_COVERAGE\t$NUMBER_OF_CONTIGS\t$NUMBER_OF_SCAFFOLDS\t$SCAFFOLD_L50\t$SCAFFOLD_N50\t$TOTAL_NUMBER_OF_CHROMOSOMES\t$TOTAL_SEQUENCE_LENGTH\t${max_gc}\t${min_gc}\t${mean_gc}\t${stddev_gc}\t${median_gc}\t${max_length}\t${min_length}\t${mean_length}\t${stddev_length}\t${median_length}" >> "$OUTPUT_FILE"
		
    fi
done

echo "Data retrieval completed. Results stored in individual files."
