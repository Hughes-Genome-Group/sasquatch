#!/bin/bash

##############################################################################################
##                                                                                          ##
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##                                                                                          ##
## Data preprocessing, background generation: submit genome-wide footprinting               ##
##                                                                                          ##
##############################################################################################

IDTAG="hg18_human_h_ery_1"

PIPE_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/pipeline/tissue_v2/background
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts

OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/${IDTAG}_minus_test

# Generated Strand-specific footprint files
FOOTPRINT_FILE_PLUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/${IDTAG}/footprints/${IDTAG}_ploidy_removed_Plus_merged.wig
FOOTPRINT_FILE_MINUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/${IDTAG}/footprints/${IDTAG}_ploidy_removed_Minus_merged.wig

# Wig border files indicating the first and last position with recorded DNase I cut sites (e.g. generate with ../scripts/get_footprint_wig_borders.pl )
BORDER_FILE_PLUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/${IDTAG}/footprints/${IDTAG}_ploidy_removed_Plus_merged_wig_borders
BORDER_FILE_MINUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/${IDTAG}/footprints/${IDTAG}_ploidy_removed_Minus_merged_wig_borders

#genome Build
BUILD='hg18'

#chose chr sizes and ploidy regions to filter
if [ "${BUILD}" == "hg19" ]
then
	REF_GENOME="/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
elif [ "${BUILD}" == "hg18" ]   
then
	REF_GENOME="/databank/raw/hg18_full/hg18_full.fa"
fi


###--------------------------
### Submit Footprint Counting 
### -------------------------
for k in 22
do
	chromosome=${k}

	prid=`qsub -N ncount_chr${chromosome}_${IDTAG} -v OUTPUT_DIR=${OUTPUT_DIR},SCRIPT_DIR=${SCRIPT_DIR},FOOTPRINT_FILE_PLUS=${FOOTPRINT_FILE_PLUS},FOOTPRINT_FILE_MINUS=${FOOTPRINT_FILE_MINUS},BORDER_FILE_PLUS=${BORDER_FILE_PLUS},BORDER_FILE_MINUS=${BORDER_FILE_MINUS},IDTAG=${IDTAG},chromosome=${chromosome},REF_GENOME=${REF_GENOME} ${PIPE_DIR}/runscript_background_processing.sh | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

	echo "Background cleavage counting Chromosome $chromosome as Job $prid submitted"

done
