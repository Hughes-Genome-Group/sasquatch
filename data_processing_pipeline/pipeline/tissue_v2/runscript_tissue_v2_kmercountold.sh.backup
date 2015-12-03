#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -e /hts/data4/rschwess/clustereo
#$ -o /hts/data4/rschwess/clustereo

cd ${OUTPUT_DIR}/counts
date

### TEST ECHO ###

echo "OUTPUT_DIR = ${OUTPUT_DIR}"
echo "REGIONS_FILE_PLOIDY_FILTERED = ${REGIONS_FILE_PLOIDY_FILTERED}"
echo "REF_GENOME = ${REF_GENOME}"
echo "IDTAG = ${IDTAG}"
echo ""

#kmer count
for i in 5 6 7
do

	perl ${SCRIPT_DIR}/count_kmers_hashed.pl /hts/data4/rschwess/scripts/dnase_tissue_motif/workflow/kmers/All_Possible_Sequences_${i}.txt ${REGIONS_FILE_PLOIDY_FILTERED} ${OUTPUT_DIR}/footprints/${IDTAG}_Plus.wig ${i} ${REF_GENOME} >${OUTPUT_DIR}/counts/kmers_${i}_count_${IDTAG}_plus.txt

	echo "$i plus done"
	date

	perl ${SCRIPT_DIR}/count_kmers_hashed_minus.pl /hts/data4/rschwess/scripts/dnase_tissue_motif/workflow/kmers/All_Possible_Sequences_${i}.txt ${REGIONS_FILE_PLOIDY_FILTERED} ${OUTPUT_DIR}/footprints/${IDTAG}_Minus.wig ${i} ${REF_GENOME} >${OUTPUT_DIR}/counts/kmers_${i}_count_${IDTAG}_minus.txt

	echo "$i minus done"
	date

done
