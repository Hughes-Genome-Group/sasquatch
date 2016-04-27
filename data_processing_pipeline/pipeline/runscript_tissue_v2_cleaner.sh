#!usr/bin/bash

##############################################################################################
##                                                                                          ##
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##                                                                                          ##
## Data preprocessing: deleate downloaded and intermediate files                            ##
##                                                                                          ##
##############################################################################################

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -j y
#$ -o /t1-data1/WTSA_Dev/rschwess/clustereo

date

#clean up
rm -f ${REGIONS_FILE}
rm -f ${OUTPUT_DIR}/*.bam
rm -f ${OUTPUT_DIR}/*.bam.bai

cd ${OUTPUT_DIR}

# optional delete footprint files
#rm -Rf ${OUTPUT_DIR}/footprints

