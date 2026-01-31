#!/bin/sh

# (c) Torsten Hugo Struck @05.09.2024 

set -o errexit

### Usage of pipeline ###
# sh CountReadlength.sh [Name of fastq file]

#Please check if the patterns in the if-statements fit with your fastq format

### Variables ###
FASTQ=$1

printf "0\t" > temp.txt

while read -r LINE
do 
	if [[ $LINE =~ ^@ ]] && [[ $LINE =~ ccs$ ]]
		then printf "$LINE\t" >> temp.txt
	elif  [[ $LINE =~ ^\+$ ]]
		then printf "$LINE\n" >> temp.txt
	else printf "${#LINE}\t" >> temp.txt
	fi
done < $FASTQ

printf "\n" >> temp.txt

cut -f3 < temp.txt > ${FASTQ}_CountReadlength.txt

rm temp.txt