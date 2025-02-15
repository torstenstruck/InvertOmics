#!/bin/sh

# (c) Torsten Hugo Struck @12.07.2021 

set -o errexit

### Usage of script ###
# sh RetrieveSpecificSequences.sh [blast programme] [database] [query]
#
# 1) blast programme: name of the blast programme to be used
# 2) database: the name of the database against which the query is blasted
# 3) query: the name of th file with the query sequences 
#
# example of usage: sh RetrieveSpecificSequences.sh blastn N61_asm.contigs.fasta Lruber_hox_genes_fmp.fasta
#
# IMPORTANT NOTE: The parameters have to be provided in this order and a value for each parameter has to be provided!!!
# The script assumes that the contig file has already been made a BLAST-database (e.g., by having used the Assembly_QualityAssessment.sh script).
# If this is not the case, please use the following command: makeblastdb -in [contig file name] -parse_seqids -dbtype nucl
# Remember to load the blast module first using: module load BLAST+/2.11.0-gompi-2020a

### Variables ###
BLAST=$1
DATABASE=$2
QUERY=$3

### Different subscripts used during the script ###
function RetrieveTaxonomyDatabase () {
# (c) Frédéric Mahé & modified by Torsten Hugo Struck @12.05.2021 adjusting to new databases in NCBI

## Download NCBI's taxonomic data and GI (GenBank ID) taxonomic
## assignation.

## Variables
local NCBI="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
local TAXDUMP="taxdump.tar.gz"
local TAXIDPATH="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
local TAXID="nucl_gb.accession2taxid.gz"
local NAMES="names.dmp"
local NODES="nodes.dmp"
local DMP=$(echo {citations,division,gencode,merged,delnodes}.dmp)
local USELESS_FILES="${TAXDUMP} ${DMP} gc.prt readme.txt"

## Download taxdump
rm -rf ${USELESS_FILES} "${NODES}" "${NAMES}"
wget "${NCBI}${TAXDUMP}" && \
    tar zxvf "${TAXDUMP}" && \
    rm -rf ${USELESS_FILES}

## Limit search space to scientific names
grep "scientific name" "${NAMES}" > "${NAMES/.dmp/_reduced.dmp}" && \
    rm -f "${NAMES}" && \
    mv "${NAMES/.dmp/_reduced.dmp}" "${NAMES}"

## Download gi_taxid_nucl
rm -f "${TAXID/.gz/}*"
wget "${TAXIDPATH}${TAXID}" && \
    gunzip "${TAXID}"
}

function RetrieveTaxonPath () {
# (c) Frédéric Mahé & modified by Torsten Hugo Struck @12.05.2021 adjusting to new databases in NCBI

local NAMES="names.dmp"
local NODES="nodes.dmp"
local GI_TO_TAXID="nucl_gb.accession2taxid"
local TAXONOMY=""
local GI="${1}"

# Obtain the name corresponding to a taxid or the taxid of the parent taxa
get_name_or_taxid()
{
    grep --max-count=1 "^${1}"$'\t' "${2}" | cut -f"${3}"
}

# Get the taxid corresponding to the GI number
TAXID=$(get_name_or_taxid "${GI}" "${GI_TO_TAXID}" "3")

# Loop until you reach the root of the taxonomy (i.e. taxid = 1)
while [[ "${TAXID}" -gt 1 ]] ; do
    # Obtain the scientific name corresponding to a taxid
    NAME=$(get_name_or_taxid "${TAXID}" "${NAMES}" "3")
    # Obtain the parent taxa taxid
    PARENT=$(get_name_or_taxid "${TAXID}" "${NODES}" "3")
    # Build the taxonomy path
    TAXONOMY="${NAME};${TAXONOMY}"
    TAXID="${PARENT}"
done

echo -e "${GI}\t${TAXONOMY}"
}

module --quiet purge
module load BLAST+/2.13.0-gompi-2020a 

# run the species detection script retrieving all contigs and unassembled reads matching mitochondrial query sequences as well as the corresponding taxon name
echo ">" >> $DATABASE
RetrieveTaxonomyDatabase
perl /storage/software/easybuild/software/BLAST+/2.11.0-gompi-2020a/bin/update_blastdb.pl --passive --decompress nt

$BLAST -evalue 1e-20 -outfmt "6 qseqid sseqid pident score evalue" -num_threads 8 -out ${QUERY}_${DATABASE}_${BLAST}_hits.txt -db $DATABASE -query $QUERY

cut -f2 < ${QUERY}_${DATABASE}_${BLAST}_hits.txt | sort | uniq | while read ID
do
echo ${ID}
sed -n "/^>${ID}/,/>/p" < $DATABASE | sed '$d' >> ${QUERY}_${DATABASE}_${BLAST}_Sequences_of_best_hit.fas
done

#Blast the retrieved sequences against NCBI nr database
blastn -max_target_seqs 1 -outfmt "6 sacc qseqid sseqid length pident score evalue" -out ${QUERY}_${DATABASE}_${BLAST}_nr_BLAST_Matches.txt -db nt -num_threads 8  -query ${QUERY}_${DATABASE}_${BLAST}_Sequences_of_best_hit.fas

#Determine the taxonomic position of the hit
echo "Retrieving taxa from blast results"
cut -f1 < ${QUERY}_${DATABASE}_${BLAST}_nr_BLAST_Matches.txt | sort | uniq | while read GI
do
RetrieveTaxonPath ${GI} >> ${QUERY}_${DATABASE}_${BLAST}_Taxa_found.txt
done

### Cleaning up folder & Sorting results ###
rm -rf nt.*
rm -rf *.dmp nodesDB.txt nucl_gb.accession2taxid
mkdir ${QUERY}_${DATABASE}_${BLAST}_Results
mv *_hits.txt *_Sequences_of_best_hit.fas *_nr_BLAST_Matches.txt *_Taxa_found.txt ${QUERY}_${DATABASE}_${BLAST}_Results/
