#!/bin/bash


IDTAG="hg18_human_JH60"

PIPE_DIR=/hts/data4/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/pipeline/tissue_v2/background
SCRIPT_DIR=/hts/data4/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts

OUTPUT_DIR=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/${IDTAG}_minus_test

FOOTPRINT_FILE_PLUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/${IDTAG}/footprints/${IDTAG}_ploidy_removed_Plus_merged.wig
FOOTPRINT_FILE_MINUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/${IDTAG}/footprints/${IDTAG}_ploidy_removed_Minus_merged.wig

BORDER_FILE_PLUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/${IDTAG}/footprints/${IDTAG}_ploidy_removed_Plus_merged_wig_borders
BORDER_FILE_MINUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/${IDTAG}/footprints/${IDTAG}_ploidy_removed_Minus_merged_wig_borders

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
