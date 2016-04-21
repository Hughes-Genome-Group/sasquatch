#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -e /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -o /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -N w4_bauer_kl7_updated_pipe

#qsub /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/workflow_scripts/Explanation_damage_multiple_seq_v3_pnorm_merged.sh

#run dissection and table creation on multiple sequences, assign a reference sequence and calulate the damage to 
#the footprint contributing kmers and rank the sequences according to that

#First set directories
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts		
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/validate_snps/bauer/query_macs_peaks

#select a background/normalization type (choose accordingly normalized kmer count files) there will pnorm dnase and pnorm atac files
#we have to check how much the different backgrounds differ to tell how much different norm types we need to have but at least one for atac and one for dnase will be needed (currently only dnase available)
# currently available "laza" (human fibroblast) "JH40" (human erythroid 40% mapped)
NORM_TYPE="h_ery_1"

ORGANISM="human"

FRAG_TYPE="DNase"

#select search mode exhaustive "e" or local "l"
SEARCH_MODE="e"

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
#for erythroid human example
TISSUE="ENCODE_CD34mobilized_UW"

#############
### INPUT ###
#############
#input file with sequences (first sequence is considered as the reference sequence,; relative to which all dmage is calculated later
sequence_list=${OUTPUT_DIR}/sas_in.txt

#TEMPORARY WRITTEN FILES'
#temporary sequence list per line
temp_sequence_list=${OUTPUT_DIR}/temp_seqs_input
#file to write temporary table with fsr calculations to
temp_fsr_file="${OUTPUT_DIR}/multiseq_fsr_table.txt"

#Specified output#
#output file listing the highest dmg and higes dmg causing variant per line
wrapper_dmg_outlist=${OUTPUT_DIR}/summary_dmg_table.txt

#length of sequences to split into
kl=7

#define DATA directory
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/${ORGANISM}/${FRAG_TYPE}/update_pipe

### BACKGROUNDs ###
case "${ORGANISM}" in 

human)

case "${FRAG_TYPE}" in 

DNase)	
#Specify the directory where the background files (naked dnaseI fibroblast cutting) are located,
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_h_ery_1/counts
#defining full paths to naked background count files
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_h_ery_1_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_h_ery_1_minus_merged
#Select accordingly normalized input files
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_${NORM_TYPE}_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_${NORM_TYPE}_minus.txt
;;
ATAC)
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH_atac/counts
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_atac_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_atac_minus_merged
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_minus.txt
;;
esac

;;

mouse)
case "${FRAG_TYPE}" in

DNase)	
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid/counts
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_mm9_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_mm9_minus_merged
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_mm9_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_mm9_minus.txt
;;
ATAC)
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid_atac/counts
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_mm9_atac_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_mm9_atac_minus_merged
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_minus.txt
;;
esac
;;

esac

#init. total running id
totalid=0

#init wrapper report file
echo -ne "tot.ID\tRef.seq\tVar.seq\tRef.kmer\tVar.kmer\tRef.SFR\tVar.SFR\ttotal.dmg\n" >${wrapper_dmg_outlist}

#read in every line of seqs_input and parse to fsr calulcation and perform dmg assessment
cat ${sequence_list} | while read line;
do

totalid=$(($totalid + 1))

echo $line | perl -ne 'my @a=split(/\s+/, $_); foreach (@a){ print $_."\n";}' >${temp_sequence_list}

#initialize internal id
id=0

#temp file
echo -n "" >$temp_fsr_file

#Go through list, calculate, stats for each and assemble in temp_table file
for sequence in `cat ${temp_sequence_list}`
do

#count up kmer id
id=$(($id + 1))

#echo "Processing ... $sequence"

echo -ne "$id\t$sequence" >>${temp_fsr_file}

# Split Sequence in kmers #
for kmer in `Rscript ${SCRIPT_DIR}/dissect_sequence.R ${sequence} ${kl} ${COMMON_FUNCTIONS}`
do

# 1 retrieve & store profiles #
#get the plus (sense) strand profile for the kmer from the tissue specific dnase footprint counts
profile_plus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_plus} ${COMMON_FUNCTIONS}`
profile_plus_count=`echo ${profile_plus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_plus=`echo ${profile_plus_out} | perl -ne '/:::(.+)/; print $1;'`
#same for antisense, minus strand, query from the tissue specific dnase footprint counts
profile_minus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_minus} ${COMMON_FUNCTIONS}`
#split into two variables
profile_minus_count=`echo ${profile_minus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_minus=`echo ${profile_minus_out} | perl -ne '/:::(.+)/; print $1;'`

#######
#merge#
#######
profile_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus} ${profile_minus} ${FRAG_TYPE} ${kl}`

# 2 normalize  and smooth profile with gaussian, estimate borders and calculate FSR score #
fsr_out=`Rscript ${SCRIPT_DIR}/fsr_calculation_pnorm_merge.R ${kmer} ${profile_merged} ${profile_plus} ${profile_minus} ${COMMON_FUNCTIONS}`
#output is again a string like: fsr_plus=x fsr_minus=x scale_plus=x scale_minus=x
fsr_merged=`echo $fsr_out} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`
fsr_plus=`echo $fsr_out} | perl -ne '/fsr_plus=(\d+\.?\d*)\s+/; print $1;'`
fsr_minus=`echo ${fsr_out} | perl -ne '/fsr_minus=(\d+\.?\d*)/; print $1;'`


echo -ne "\t${fsr_merged}" >>${temp_fsr_file}

done

echo "" >>${temp_fsr_file}

done

#run assessing dmg procedure
#either summing up 
if [ "$SEARCH_MODE" == "e" ]
then

  dmgout=`Rscript ${SCRIPT_DIR}/multi_seq_dmg_v3_summing_up.R ${temp_fsr_file} ${kl} ${totalid}`

#or single higehst scoring kmer dependend on run mode
elif [ "$SEARCH_MODE" == "l" ]
then

    dmgout=`Rscript ${SCRIPT_DIR}/multi_seq_dmg_v3.R ${temp_fsr_file} ${kl} ${totalid}`

else

  echo "Specify a proper Run mode!"
  exit 2;

fi

echo -ne "${totalid}\t$dmgout\n" >>${wrapper_dmg_outlist}

echo "$totalid done"

done

echo "$totalid done"

done

Rscript ${SCRIPT_DIR}/w4_order_summary.R ${wrapper_dmg_outlist} ${OUTPUT_DIR}/ordered


rm -f ${temp_sequence_list}
rm -f ${temp_fsr_file}


