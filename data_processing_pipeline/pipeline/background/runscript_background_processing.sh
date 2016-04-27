#!usr/bin/bash

##############################################################################################
##                                                                                          ##
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##                                                                                          ##
## Data preprocessing, background: generate kmer based average profiles from naked DNA      ##
##                                                                                          ##
##############################################################################################

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -e /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -o /t1-data1/WTSA_Dev/rschwess/clustereo

#TEST ECHO
echo "OUTPUT_DIR	${OUTPUT_DIR}"
echo "FOOTPRINT_FILE_PLUS	${FOOTPRINT_FILE_PLUS}"
echo "FOOTPRINT_FILE_MINUS	${FOOTPRINT_FILE_MINUS}"
echo "BORDER_FILE_PLUS	 ${BORDER_FILE_PLUS}"
echo "BORDER_FILE_MINUS	 ${BORDER_FILE_MINUS}"
echo "IDTAG	${IDTAG}"
echo "chromosome	${chromosome}"
echo "REF_GENOME   ${REF_GENOME}"


#START

mkdir -p ${OUTPUT_DIR}/counts
cd ${OUTPUT_DIR}/counts

for i in 5 6 7
do

date 

grep -P "chr${chromosome}\s+" ${BORDER_FILE_PLUS} | perl ${SCRIPT_DIR}/count_kmers_naked_hashed.pl /t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/kmers/All_Possible_Sequences_${i}.txt ${FOOTPRINT_FILE_PLUS} ${i} - ${REF_GENOME} ${OUTPUT_DIR}/counts/kmer_${i}_${IDTAG}_plus_chromosome${chromosome}

echo "$i plus done"
date

grep -P "chr${chromosome}\s+" ${BORDER_FILE_MINUS} | perl ${SCRIPT_DIR}/count_kmers_naked_hashed_minus.pl /t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/kmers/All_Possible_Sequences_${i}.txt ${FOOTPRINT_FILE_MINUS} ${i} - ${REF_GENOME} ${OUTPUT_DIR}/counts/kmer_${i}_${IDTAG}_minus_chromosome${chromosome}

echo "$i minus done"
date

done





