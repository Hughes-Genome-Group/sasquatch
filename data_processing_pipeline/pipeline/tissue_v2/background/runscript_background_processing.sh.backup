#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -e /hts/data4/rschwess/clustereo
#$ -o /hts/data4/rschwess/clustereo

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

grep -P "chr${chromosome}\s+" ${BORDER_FILE_PLUS} | perl ${SCRIPT_DIR}/count_kmers_naked_hashed.pl /hts/data4/rschwess/dnase_motif_tissue/kmers/All_Possible_Sequences_${i}.txt ${FOOTPRINT_FILE_PLUS} ${i} - ${REF_GENOME} ${OUTPUT_DIR}/counts/kmer_${i}_${IDTAG}_plus_chromosome${chromosome}

echo "$i plus done"
date

grep -P "chr${chromosome}\s+" ${BORDER_FILE_MINUS} | perl ${SCRIPT_DIR}/count_kmers_naked_hashed_minus.pl /hts/data4/rschwess/dnase_motif_tissue/kmers/All_Possible_Sequences_${i}.txt ${FOOTPRINT_FILE_MINUS} ${i} - ${REF_GENOME} ${OUTPUT_DIR}/counts/kmer_${i}_${IDTAG}_minus_chromosome${chromosome}

echo "$i minus done"
date

done





