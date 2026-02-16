#!/bin/sh

# (c) Torsten Hugo Struck @23.04.2025, based on a script by Daren Card (i.e., https://darencard.net/blog/2022-07-09-genome-repeat-annotation/) and with the aid of ChatGPT

set -o errexit

### Usage of script ###
#This script does a de nove repeat modeling of a given genome
#Usage: sh RepeatMasking_v1.0.sh [Name of database] [Name of fasta file of assembled genome] [Number of threads] [Name of taxon for Dfam database] [Genome size]
#
# 1) Name of database: Name of the database generated with RepeatModeling, must be in the same folder from which the script is started
# 2) Name of fasta file of assembled genome: Name of the fasta file with the genome that shall be used for repeat modelling (needs to end on .fasta)
# 3) Number of threads: Number of threads used for the analyses
# 4) Name of taxon for Dfam database: Specify the species or taxon of the input sequence. The name must be a valid NCBI Taxonomy Database name and be contained in the Dfam database. Some examples are:
#	human
#	mouse
#	rattus
#	"ciona savignyi"
#	arabidopsis
#	Metazoa
#Other commonly used species: mammal, carnivore, rodentia, rat, cow, pig, cat, dog, chicken, fugu, danio, "ciona intestinalis",  drosophila, anopheles, worm, diatoaea, artiodactyl, arabidopsis, rice, wheat, and maize
#NOTE: The suggestion is to use a higher level such as Arthropoda, Metazoa, or Eukaryota to be more inclusive and exploiting a broader set of curated libraries for pur non-model organisms. 
# 5) Genome size: Exact genome size in bp of the assembled genome used
#
#Example of usage: sh RepeatModeling_v1.0.sh Aegialites_55802 Aeg55_final_pri_asm.fasta 16 Metazoa 100825933
#
# IMPORTANT NOTE: The parameters have to be provided in this order and a value for each parameter has to be provided!!!
# IMPORTANT NOTE: If you do changes to the script ALWAYS save the file under a different name. Do NOT modify this master script.
#
#VERSION HISTORY:
#v1.1: Added analyses of the interspersed repeat landscape

#load the RepeatModeler module
module --quiet purge
module load RepeatModeler/2.0.6-foss-2023a

#add an additional path for perl modules to the PERLLIB path
#export PERL5LIB=$PERL5LIB:/storage/InvertOmics/00_Scripts/02_PerlModules

#setting the variables
database=$1 
assembly=$2
redassembly=${assembly//".fasta"/}
threads=$3
num_pa=$((${threads}/4)) #this is necessary to account for the RMBlast search requiring 4 CPUs per independent run
taxon=$4
genome_size=$5

#generate results directories
mkdir -p MaskOut
cd MaskOut
mkdir -p 01_simple_out 02_taxon_out 03_specific_out
cd ..

#first masking round
#annotate/mask simple repeats 
RepeatMasker -pa $num_pa -a -e rmblast -dir MaskOut/01_simple_out -noint -xsmall -gff $assembly > MaskOut/01_simple_out/01_simplemask.log
#rename outputs
rename fasta simple_mask MaskOut/01_simple_out/$redassembly*
rename .masked .masked.fasta MaskOut/01_simple_out/$redassembly*

#second masking round
#annotate/mask elements sourced from Dfam database using output from 1st round of RepeatMasker and restricted to the specified taxon
RepeatMasker -pa $num_pa -a -e rmblast -dir MaskOut/02_taxon_out -nolow -species $taxon -gff MaskOut/01_simple_out/$redassembly.simple_mask.masked.fasta > MaskOut/02_taxon_out/02_taxonmask.log
#rename outputs
rename simple_mask.masked.fasta taxon_mask MaskOut/02_taxon_out/$redassembly*
rename .masked .masked.fasta MaskOut/02_taxon_out/$redassembly*

#third masking round
#annotate/mask known elements sourced from species-specific de novo repeat library using output froom 2nd round of RepeatMasker
RepeatMasker -pa $num_pa -a -e rmblast -dir MaskOut/03_specific_out -nolow -lib $database-families.fa -gff MaskOut/02_taxon_out/$redassembly.taxon_mask.masked.fasta > MaskOut/03_specific_out/03_specificmask.log
#rename outputs
rename taxon_mask.masked.fasta specific_mask MaskOut/03_specific_out/$redassembly*
rename .masked .masked.fasta MaskOut/03_specific_out/$redassembly*

##create directory for full results
cd MaskOut
mkdir -p 04_full_out

#combine full RepeatMasker result files - .cat.gz
cat 01_simple_out/$redassembly.simple_mask.cat.gz 02_taxon_out/$redassembly.taxon_mask.cat.gz 03_specific_out/$redassembly.specific_mask.cat.gz > 04_full_out/$redassembly.full_mask.cat.gz

#combine RepeatMasker tabular files for all repeats - .out
cat 01_simple_out/$redassembly.simple_mask.out > 04_full_out/$redassembly.full_mask.out
tail -n +4 02_taxon_out/$redassembly.taxon_mask.out >> 04_full_out/$redassembly.full_mask.out
tail -n +4 03_specific_out/$redassembly.specific_mask.out >> 04_full_out/$redassembly.full_mask.out

#copy RepeatMasker tabular files for simple repeats - .out
cat 01_simple_out/$redassembly.simple_mask.out > 04_full_out/$redassembly.simple_mask.out

#combine RepeatMasker tabular files for complex, interspersed repeats - .out
cat 02_taxon_out/$redassembly.taxon_mask.out > 04_full_out/$redassembly.complex_mask.out
tail -n +4 03_specific_out/$redassembly.specific_mask.out >> 04_full_out/$redassembly.complex_mask.out

#combine RepeatMasker repeat alignments for all repeats - .align
cat 01_simple_out/$redassembly.simple_mask.align 02_taxon_out/$redassembly.taxon_mask.align 03_specific_out/$redassembly.specific_mask.align > 04_full_out/$redassembly.full_mask.align

cd ..

#resummarize repeat compositions from combined analysis of all RepeatMasker rounds
ProcessRepeats -a -species $taxon -maskSource $assembly -xsmall -gff MaskOut/04_full_out/$redassembly.full_mask.cat.gz > MaskOut/04_full_out/04_fullmask.log

cd MaskOut/04_full_out/

#generate .divsum file using RepeatMasker tools
echo "Running calcDivergenceFromAlign.pl..."
perl /storage/software/easybuild/software/RepeatMasker/4.1.7-p1-foss-2023a/util/calcDivergenceFromAlign.pl -s $redassembly.full_mask.divsum $redassembly.full_mask.align

#create the landscape using RepeatMasker tools
echo "Generating repeat landscape data..."
perl /storage/software/easybuild/software/RepeatMasker/4.1.7-p1-foss-2023a/util/createRepeatLandscape.pl -div $redassembly.full_mask.divsum -g $genome_size > $redassembly.full_mask_landscape.html
perl /storage/software/easybuild/software/RepeatMasker/4.1.7-p1-foss-2023a/util/createRepeatLandscape.pl -div $redassembly.full_mask.divsum -g $genome_size -j > $redassembly.full_mask_landscape.js

#parse out interspersed landscape table information from java output to generate more tailored displays
#initialize an empty array for column names and for rows data
declare -a column_names
declare -a rows_data

#parse out genome proportion vs. divergence
#parse column names
while IFS= read -r line; do
    if [[ "$line" == *"data.addColumn("* ]]; then
        # Extract column names
        # Remove 'data.addColumn('string', ' or 'data.addColumn('number', ' and ');'
        name=$(echo "$line" | sed -E "s/data\.addColumn\('string', '([^']+)'\);/\1/; s/data\.addColumn\('number', '([^']+)'\);/\1/")
        column_names+=("$name")
    fi
done < $redassembly.full_mask_landscape.js

#parse row data
reading_rows=false
while IFS= read -r line; do
    if [[ "$line" == *"data.addRows("* ]]; then
        reading_rows=true
        continue
    fi

    if [[ "$reading_rows" == true ]]; then
        if [[ "$line" == *"var pieData"* ]]; then
            # End of rows
            break
        fi
        
        # Remove unwanted characters and split by space
        cleaned_line=$(echo "$line" | sed -E 's/[^0-9. ]//g' | tr -s ' ' '\t')
        rows_data+=("$cleaned_line")
    fi
done < $redassembly.full_mask_landscape.js

#print the column names in a TSV format at the beginning of the output file
{
    # Print the header (column names)
    (IFS=$'\t'; echo "${column_names[*]}")

    # Print each row data
    for row in "${rows_data[@]}"; do
        echo -e "$row"
    done
} > $redassembly.full_mask_landscape.tsv

echo "Data has been successfully extracted to $redassembly.full_mask_landscape.tsv"

#substitute the multiple spaces separating the column names by tab character
sed -i 's/  */\t/g' $redassembly.full_mask_landscape.tsv
sed -i 's/\t\+/\t/g' $redassembly.full_mask_landscape.tsv

#remove empty lines from the file
sed -i '/^[ \t]*$/d' $redassembly.full_mask_landscape.tsv

#remove the beginning tab at each line
sed -i 's/^\t//' $redassembly.full_mask_landscape.tsv

#parse out repeat amounts
echo "Class,bp" > $redassembly.full_mask_amounts.tmp

#parse row data
reading_rows2=false
while IFS= read -r line; do
    if [[ "$line" == *"pieData.addRows"* ]]; then
        reading_rows2=true
        continue
    fi

    if [[ "$reading_rows2" == true ]]; then
        if [[ "$line" == *"var pieOptions"* ]]; then
            # End of rows
            break
        fi
        
        # Remove unwanted characters and split by space
        cleaned_line=$(echo "$line" | sed "s/[][]//g; s/'//g")
        echo -e "$cleaned_line" >> $redassembly.full_mask_amounts.tmp
    fi
done < $redassembly.full_mask_landscape.js

echo "Data has been successfully extracted to $redassembly.full_mask_amounts.csv"

#remove empty lines from the file
tr -d ' \t' < $redassembly.full_mask_amounts.tmp > $redassembly.full_mask_amounts.csv
sed -i 's/,$//' $redassembly.full_mask_amounts.csv
sed -i 's/);//' $redassembly.full_mask_amounts.csv
rm $redassembly.full_mask_amounts.tmp

#generate summary graphics using R
module load R/4.2.1-foss-2022a

#define the input file and output file names
INPUT_FILE=$redassembly.full_mask_landscape.tsv 

#create an R script to generate histograms
cat <<EOF > generate_histograms.R
#load required library
library(ggplot2)
library(reshape2)

#load the data
data <- read.table("$INPUT_FILE", sep="\t", header=TRUE)

#generate histograms for each column against "Divergence"
columns <- names(data)[-1]  # Exclude the first column "Divergence"
for (col in columns) {
    p <- ggplot(data, aes_string(x="Divergence", y=col)) +
        geom_histogram(stat="identity", position="identity", alpha=0.6) +
        labs(title=paste("Interspersed repeat landscape of ", col), x="Kimura subsitution level (CpG adjusted)", y="Percentage of genome") +
        theme_minimal()
    
    #save as PDF
    ggsave(filename=paste0("Landscape_", col, ".pdf"), plot=p, width=8, height=6)
    
    #save as PNG
    ggsave(filename=paste0("Landscape_", col, ".png"), plot=p, width=8, height=6)
}

#melt the data to long format
  melted_data <- melt(data, id.vars="Divergence")

#create a stacked bar chart using actual values
  stacked_histogram <- ggplot(melted_data, aes(x=Divergence, y=value, fill=variable)) +
  geom_bar(stat="identity", position="stack", alpha=0.6) +
  labs(title="Interspersed repeat landscape", x="Kimura subsitution level (CpG adjusted)", y="Percentage of genome") +
  theme_minimal()

#save stacked histogram as PDF
ggsave(filename="Landscape_all_repeats.pdf", plot=stacked_histogram, width=8, height=6)

#save stacked histogram as PNG
ggsave(filename="Landscape_all_repeats.png", plot=stacked_histogram, width=8, height=6)

EOF

#run the R script
Rscript generate_histograms.R

#clean up
rm generate_histograms.R

echo "Bar plots generated and saved as PDF and PNG."

#define the input file and output file names
INPUT_PIE=$redassembly.full_mask_amounts.csv 

#create an R script to generate histograms
cat <<EOF > generate_pie.R
# Load necessary libraries
library(ggplot2)

# Read the CSV file
data <- read.csv("$INPUT_PIE")

# Calculate percentage
data\$percentage <- data\$bp / sum(data\$bp) * 100

# Format labels: e.g., "A (25.00%)"
data\$ClassLabel <- paste0(data\$Class, " (", formatC(data\$percentage, format = "f", digits = 2), "%)")

# Generate pie chart using percentages and new legend labels
p <- ggplot(data, aes(x = "", y = bp, fill = ClassLabel)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "right") +
  ggtitle("Genomic fraction")

# Save as PDF
ggsave("piechart.pdf", plot = p, width = 6, height = 6)

# Save as PNG
ggsave("piechart.png", plot = p, width = 6, height = 6, dpi = 300)
EOF

#run the R script
Rscript generate_pie.R

#clean up
rm generate_pie.R

echo "Pie chart generated and saved as PDF and PNG."

cd ../../