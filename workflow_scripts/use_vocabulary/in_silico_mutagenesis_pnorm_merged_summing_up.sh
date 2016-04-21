#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -e /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -o /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -N inslico_mut_human_r2_hery

#qsub /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/workflow_scripts/use_vocabulary/in_silico_mutagenesis_pnorm_merged_summing_up.sh

#First set directories
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts		
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/tool_compare/sliding_mutation_compare_human_R2

#select a background/normalization type (choose accordingly normalized kmer count files) there will pnorm dnase and pnorm atac files
#we have to check how much the different backgrounds differ to tell how much different norm types we need to have but at least one for atac and one for dnase will be needed (currently only dnase available)
# currently available "laza" (human fibroblast) "JH40" (human erythroid 40% mapped)
NORM_TYPE="h_ery_1"

ORGANISM="human"

FRAG_TYPE="DNase"

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
#for erythroid human example
TISSUE="ENCODE_K562_UW_merged"

#############
### INPUT ###
#############
#input file with sequences (first sequence is considered as the reference sequence,; relative to which all dmage is calculated later
sequence_list=${OUTPUT_DIR}/r2_sasq_query

#TEMPORARY WRITTEN FILES'
#temporary sequence list per line
temp_sequence_list=${OUTPUT_DIR}/temp_seqs_input
#file to write temporary table with fsr calculations to
temp_fsr_file="${OUTPUT_DIR}/multiseq_fsr_table.${TISSUE}.txt"

#Final output#
#output file listing the highest dmg and higes dmg causing variant per line
wrapper_dmg_outlist=${OUTPUT_DIR}/summary_dmg_table_summing_up.${TISSUE}.txt

#length of sequences to split into
kl=7

#define DATA directory
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/${ORGANISM}/${FRAG_TYPE}

### Take Vocabulary File
VOCAB="${DATA_DIR}/${TISSUE}/vocabulary_${TISSUE}.txt"

### BACKGROUNDs ###
case "${ORGANISM}" in 
human)
case "${FRAG_TYPE}" in 
DNase)	
#Select accordingly normalized input files
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_${NORM_TYPE}_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_${NORM_TYPE}_minus.txt
;;
ATAC)
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_minus.txt
;;
esac
;;
mouse)
case "${FRAG_TYPE}" in
DNase)	
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_mm9_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_mm9_minus.txt
;;
ATAC)
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_minus.txt
;;
esac
;;

esac

#init. total running id
totalid=0
#init wrapper report file
echo -ne "tot.ID\ttotal.dmg\n" >${wrapper_dmg_outlist}

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

echo -ne "$id\t$sequence" >>${temp_fsr_file}

# Split Sequence in kmers #
for kmer in `Rscript ${SCRIPT_DIR}/dissect_sequence.R ${sequence} ${kl} ${COMMON_FUNCTIONS}`
do

##Cover Ns uncertanty#
case "$kmer" in 
	#if N found set total DMG to zero and continue with next kmer as there is nothing that can be calculated here
	*N*) echo -ne "${totalid}\tNA\n" >>${wrapper_dmg_outlist}
		break
	;;  
esac

##check if kmer already processed##
catch=`grep $kmer ${VOCAB}`

#string is not empty; kmer already found new line is old line 
sfr=`echo $catch | perl -ne '/^[A,C,G,T]+\s+(\d+\.?\d*)/; print $1;'`
#store merged i temp fsr file
echo -ne "\t${sfr}" >>${temp_fsr_file}

done

echo "" >>${temp_fsr_file}

done

dmgout=`Rscript ${SCRIPT_DIR}/multi_seq_dmg_insilico_mutagenesis_summing_up.R ${temp_fsr_file} ${kl} ${totalid}`
echo -ne "${totalid}\t$dmgout\n" >>${wrapper_dmg_outlist}

echo "$totalid bps processed"

done


rm -f ${temp_sequence_list}
rm -f ${temp_fsr_file}


