#!/bin/sh

# (c) Torsten Hugo Struck @23.04.2025, based on a script by Pia Merete Eriksen (i.e., repeat.sh) 

set -o errexit

### Usage of script ###
#This script does a de novo repeat modeling of a given genome
#Usage: sh RepeatModeling_v1.0.sh [Name of database] [Name of fasta file of assembled genome] [Number of threads] [Number of additional rounds]
#
# 1) Name of database: Name of the database that will be generated and then used; it can be any name (advisable is genus_species)
# 2) Name of fasta file of assembled genome: Name of the fasta file with the genome that shall be used for repeat modelling
# 3) Number of threads: Number of threads used for the analyses
# 4) Number of additional rounds: Number of additional rounds in the Repeatmodeling
#
#Example of usage: sh RepeatModeling_v1.0.sh Aegialites_55802 Aeg55_final_pri_asm.fa 16 2
#
# IMPORTANT NOTE: The parameters have to be provided in this order and a value for each parameter has to be provided!!!
# IMPORTANT NOTE: If you do changes to the script ALWAYS save the file under a different name. Do NOT modify this master script.

#load the RepeatModeler module
module --quiet purge
module load RepeatModeler/2.0.6-foss-2023a

#add an additional path for perl modules to the PERLLIB path
export PERL5LIB=$PERL5LIB:/storage/InvertOmics/00_Scripts/02_PerlModules

#setting the variables
database=$1 
assembly=$2
threads=$3
addrounds=$4

#Build the database
BuildDatabase -name $database $assembly

#Run repeat modeler
#The option LTRStruct runs the LTR structural discovery pipeline (LTR_Harvest and LTR_retreiver) and combines results with the RepeatScout/RECON pipeline
RepeatModeler -database $database -threads $threads -LTRStruct -numAddlRounds $addrounds
