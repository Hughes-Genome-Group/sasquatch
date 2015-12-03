#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -e /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -o /t1-data1/WTSA_Dev/rschwess/clustereo

date

#clean up
rm -f ${REGIONS_FILE}
rm -f ${OUTPUT_DIR}/*.bam
rm -f ${OUTPUT_DIR}/*.bam.bai

#cd ${OUTPUT_DIR}
#rm -Rf ${OUTPUT_DIR}/footprints

