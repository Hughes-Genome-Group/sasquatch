#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -e /hts/data4/rschwess/clustereo
#$ -o /hts/data4/rschwess/clustereo

cd ${COUNTS}
date

echo "OUTPUT_DIR = ${OUTPUT_DIR}"
echo "REGIONS_FILE_PLOIDY_FILTERED = ${REGIONS_FILE_PLOIDY_FILTERED}"
echo "REF_GENOME = ${REF_GENOME}"
echo "IDTAG = ${IDTAG}"
echo "PROPENSITY_PLUS = ${PROPENSITY_PLUS}"
echo "PROPENSITY_MINUS = ${PROPENSITY_MINUS}"
echo "Pnorm Propensities Source = ${pnormsource}"
echo "Count dir = ${COUNTS}"
echo ""

#kmer count
for i in 5 6 7
do

	perl ${SCRIPT_DIR}/count_kmers_pnorm.pl /hts/data4/rschwess/scripts/dnase_tissue_motif/workflow/kmers/All_Possible_Sequences_${i}.txt ${REGIONS_FILE_PLOIDY_FILTERED} ${OUTPUT_DIR}/footprints/${IDTAG}_Plus.wig ${i} ${REF_GENOME} ${PROPENSITY_PLUS} >${COUNTS}/kmers_${i}_count_${IDTAG}_pnorm_${pnormsource}_plus.txt


	echo "$i plus done"
	date

	perl ${SCRIPT_DIR}/count_kmers_pnorm_minus.pl /hts/data4/rschwess/scripts/dnase_tissue_motif/workflow/kmers/All_Possible_Sequences_${i}.txt ${REGIONS_FILE_PLOIDY_FILTERED} ${OUTPUT_DIR}/footprints/${IDTAG}_Minus.wig ${i} ${REF_GENOME} ${PROPENSITY_MINUS} >${COUNTS}/kmers_${i}_count_${IDTAG}_pnorm_${pnormsource}_minus.txt

	echo "$i minus done"
	date

done


