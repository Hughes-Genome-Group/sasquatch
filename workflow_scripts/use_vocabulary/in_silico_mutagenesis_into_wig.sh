#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -e /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -o /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -N mutagen_gataa_break_hery

#qsub /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/workflow_scripts/use_vocabulary/in_silico_mutagenesis_into_wig.sh

#run dissection and table creation on multiple sequences, assign a reference sequence and calulate the damage to 
#the footprint contributing kmers and rank the sequences according to REGIONS

#First set directories
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts	
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/tool_compare/bp_res_motif_disruption

mkdir -p ${OUTPUT_DIR}

#select a background/normalization type (choose accordingly normalized kmer count files) there will pnorm dnase and pnorm atac files
#we have to check how much the different backgrounds differ to tell how much different norm types we need to have but at least one for atac and one for dnase will be needed (currently only dnase available)
# currently available "laza" (human fibroblast) "JH40" (human erythroid 40% mapped)
NORM_TYPE="JH60"

ORGANISM="human"

FRAG_TYPE="DNase"

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
#for erythroid human example
TISSUE="human_erythroid_hg18"

#############
### INPUT ###
#############

#The QUERY#
#build
build="hg19"	

#chr
chromosome="chr1"

#input BED regions to process
IN_BED=${OUTPUT_DIR}/hg19_bed_in.bed

#Final WIG OUTPUT#
#output WIG file listing the highest dmg and higest dmg per chr pos
wig_abs_file=${OUTPUT_DIR}/sasq_insilico_mutagenesis_abs_${TISSUE}_${BUILD}_${chromosome}.wig
wig_max_file=${OUTPUT_DIR}/sasq_insilico_mutagenesis_max_${TISSUE}_${BUILD}_${chromosome}.wig

#init WIG file
echo -ne "track type=wiggle_0 name=SasQ_abs_DMG\n" >${wig_abs_file}
echo -ne "variableStep chrom=$chromosome\n" >>${wig_abs_file}

echo -ne "track type=wiggle_0 name=SasQ_max_DMG\n" >${wig_max_file}
echo -ne "variableStep chrom=$chromosome\n" >>${wig_max_file}

#length of sequences to split into
kl=7

#TEMPORARY WRITTEN FILES'
#file to write temporary table with fsr calculations to
#init temp kmer store file
echo -n "" >${OUTPUT_DIR}/temp_kmer_store
temp_sequence_list=${OUTPUT_DIR}/temp_seqs_input
temp_fsr_file=${OUTPUT_DIR}/multiseq_fsr_table.txt
wrapper_dmg_outlist=${OUTPUT_DIR}/summary_dmg_table.txt

#define DATA directory
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/${ORGANISM}/${FRAG_TYPE}/
### Take Vocabulary File
VOCAB="${DATA_DIR}/${TISSUE}/vocabulary_${TISSUE}.txt"

### BACKGROUNDs ###
case "${ORGANISM}" in 

human)
case "${FRAG_TYPE}" in 
DNase)	
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

echo "starting at:"
date

cat ${IN_BED} | while read line;
do

#chr pos to start
start=`echo ${line} | perl -ne '$_=~/^chr\d+\s+(\d+)/; print $1;'`
#chr pos to stop
stop=`echo ${line} | perl -ne '$_=~/\s+(\d+)$/; print $1;'`


#init. total running id
totalid=$start
stopid=$(( $start+12 ))

position=$(( $start+6 ))

until [ $stopid -ge $stop ]
do

#clear wrapper dmg outlist
echo -n "" >${wrapper_dmg_outlist}

###
perl ${SCRIPT_DIR}/arrange_insilico_mutagenesis_queries.pl `perl ${SCRIPT_DIR}/get_ref_sequence.pl ${build} ${chromosome} ${totalid} ${stopid} | tail -n 1` ${kl} | cut -f 3,4 >${OUTPUT_DIR}/seq_list
###

#read in every line of seqs_input and parse to fsr calulcation and perform dmg assessment
cat ${OUTPUT_DIR}/seq_list | while read line;
do

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

if [ ! "$catch" == "" ]; 
then 
	#string is not empty; kmer already found new line is old line 
	sfr=`echo $catch | perl -ne '/^[A,C,G,T]+\s+(\d+\.?\d*)/; print $1;'`
	#store merged i temp fsr file
	echo -ne "\t${sfr}" >>${temp_fsr_file}

else 

##ELSE## do the whole show
profile_plus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_plus} ${COMMON_FUNCTIONS}`
profile_plus=`echo ${profile_plus_out} | perl -ne '/:::(.+)/; print $1;'`
profile_minus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_minus} ${COMMON_FUNCTIONS}`
profile_minus=`echo ${profile_minus_out} | perl -ne '/:::(.+)/; print $1;'`
#merge#
profile_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus} ${profile_minus} ${FRAG_TYPE} ${kl}`
# 2 normalize  and smooth profile with gaussian, estimate borders and calculate FSR score #
fsr_out=`Rscript ${SCRIPT_DIR}/fsr_calculation_pnorm_merge.R ${kmer} ${profile_merged} ${profile_plus} ${profile_minus} ${COMMON_FUNCTIONS}`
#output is again a string like: fsr_plus=x fsr_minus=x scale_plus=x scale_minus=x
fsr_merged=`echo $fsr_out} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`

#store merged i temp fsr file
echo -ne "\t${fsr_merged}" >>${temp_fsr_file}

fi

done

echo "" >>${temp_fsr_file}

done

dmgout=`Rscript ${SCRIPT_DIR}/multi_seq_dmg_insilico_mutagenesis.R ${temp_fsr_file} ${kl} ${totalid}`
echo -ne "${totalid}\t$dmgout\n" >>${wrapper_dmg_outlist}

done

#get highest absolute of the variables and assemble wig entry
perl /t1-data1/WTSA_Dev/rschwess/scripts/dnase_tissue_motif/feed_highest_max_sasq_dmg_to_wig.pl $position ${wrapper_dmg_outlist} >>${wig_max_file}

perl /t1-data1/WTSA_Dev/rschwess/scripts/dnase_tissue_motif/feed_highest_abs_sasq_dmg_to_wig.pl $position ${wrapper_dmg_outlist} >>${wig_abs_file}

totalid=$(( $totalid + 1))
stopid=$(( $stopid + 1))
position=$(( $position + 1))
#until loop finished
done
done


echo "finished at:"
date

rm -f ${temp_sequence_list}
rm -f ${temp_fsr_file}


