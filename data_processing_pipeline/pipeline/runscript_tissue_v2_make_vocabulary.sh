#!usr/bin/bash

##############################################################################################
##                                                                                          ##
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##                                                                                          ##
## Data preprocessing: generate vocabulary with SFR per kmer                                ##
##                                                                                          ##
##############################################################################################

# Options: for Sun Grid Engine Cluster Batch Queue System
#$ -cwd
#$ -q batchq
#$ -m eas
#$ -j y
#$ -N saq_voc

#################
# Requirements: #
#################
#
# Does not require a sun grid engine cluster set-up
# Need Sasquatch scripts and the tissue of interest 
# preprocessed to have the kmer fiiles in place

###########
# Run as: #
# 
# sh   /path_to_sasquatch/Sasquatch/data_processing_pipeline/pipeline/runscript_tissue_v2_make_vocabulary.sh
# qsub /path_to_sasquatch/Sasquatch/data_processing_pipeline/pipeline/runscript_tissue_v2_make_vocabulary.sh

###############################
# Adjust Paths and Parameters #
###############################

# Set PIPELINE Paths ============================================= #
# Full path to your Sasquatch copy
SASQ_PATH=${SASQ_PATH}/path_to_sasquatch/Sasquatch

# Path to you Sasquatch R funtions copy 
COMMON_FUNCTIONS=${SASQ_PATH}/R_utility/functions_sasq_r_utility.R

# Path to all possible kmer files
KMER_IN=${SASQ_PATH}/data_processing_pipeline/kmers/Kmers_5_6_7_combined.txt

# Set basi parameters ============================================ #
# Organism [human, mouse]
ORGANISM="human"

# Defautlt "DNase"
DATA_TYPE="DNase"

#Tissue ID (Dsub directory nae)
TISSUE="my_new_tissue_rep1"

# Full path to where you stored your preprocessed data [default = "${SASQ_PATH}/${ORGANISM}/${DATA_TYPE}/"]
DATA_DIR=${SASQ_PATH}/${ORGANISM}/${DATA_TYPE}/${IDTAG}

#########################
# Auto Select pnorm tag #
#########################

# Select PNORM TAG
case "${ORGANISM}" in 

	human)
		PNORM_TAG="h_ery_1"
	;;

	mouse)
		PNORM_TAG="m_ery_1"
	;;

esac

# SET OUTPUT
OUTPUT_DIR=${DATA_DIR}/${TISSUE}

VOCAB=${OUTPUT_DIR}/vocabulary_${TISSUE}.txt

# override existing
#rm -f ${VOCAB}

# Submit Rscript ====================================

Rscript ${SCRIPT_DIR}/create_vocabulary.R ${COMMON_FUNCTIONS} ${TISSUE} ${ORGANISM} ${PNORM_TAG} ${DATA_DIR} ${DATA_TYPE} ${VOCAB} ${KMER_IN}

# =====================================
echo "Processing vocabulary file done"
