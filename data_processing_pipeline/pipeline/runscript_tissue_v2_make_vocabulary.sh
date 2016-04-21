#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M ron.schwessinger@ndcls.ox.ac.uk
#$ -m eas
#$ -j y
#$ -e /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -o /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -N voc_12878_uw_enc

#qsub /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/pipeline/runscript_tissue_v2_make_vocabulary.sh

SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts
COMMON_FUNCTIONS=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/R_utility/functions_sasq_r_utility.R

#input possible kmers
KMER_IN=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/kmers/Kmers_5_6_7_combined.txt

ORGANISM="human"

FRAG_TYPE="DNase"

#Tissue ID
TISSUE="ENCODE_UW_GM12878_merged"

#define DATA directory
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/${ORGANISM}/${FRAG_TYPE}/

# Select PNORM TAG
case "${ORGANISM}" in 

	human)
		case "${FRAG_TYPE}" in 
		
		DNase)	
			PNORM_TAG="h_ery_1"
		;;
		ATAC)
			PNORM_TAG="atac"
		;;
		esac
	;;

	mouse)
		case "${FRAG_TYPE}" in
		
		DNase)	
			PNORM_TAG="mm9"
		;;
		ATAC)
			PNORM_TAG="mm9_atac"
		;;
		esac
	;;

esac


# SET OUTPUT
OUTPUT_DIR=${DATA_DIR}/${TISSUE}

VOCAB=${OUTPUT_DIR}/vocabulary_${TISSUE}.txt

# override existing
#rm -f ${VOCAB}

# Submit Rscript ====================================

Rscript ${SCRIPT_DIR}/create_vocabulary.R ${COMMON_FUNCTIONS} ${TISSUE} ${ORGANISM} ${PNORM_TAG} ${DATA_DIR} ${FRAG_TYPE} ${VOCAB} ${KMER_IN}

# =====================================
echo "Processing _ vocabulary file done"