#!usr/bin/bash

#run dissection and table creation on multiple sequences, assign a reference sequence and calulate the damage to 
#the footprint contributing kmers and rank the sequences according to that

#First set directories
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts			
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/idx_duke_testout/merge_strands

BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_h_ery_1/counts
BASE_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/background/hg19_ploidy_removed/basecomp

#Specify the location of the tissue specific data (for example human erythroid is located here):
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human_ploidy_correct

#I use the following TAG
TAG="human_erythroid_hg18"		

#Arguments for the single steps as in single run 
#Arguments for HSCor
#bandwidth to use for gaussian smoothing (default = 10) 
bandwidth_hsc=10
#1=plot 0= no plotting
hsc_plot_flag=0
#dir to write plot if flagged only required if flagged
hsc_plot_dir=${OUTPUT_DIR}

#Arguments for DNase Background Cutting
#fraction to calculate iqr (default 4 = 25% IQR)
dnase_iqr_fraction=4	
#MSQ threshold to decide if to normalize, if to extend the norm window or if to flag as not normalized, default=0.09 but its empirically adjusted to 1-2 datasets so should somehow be adjustable later on
msq_scale_thresh=0.09
#maximal extension allowed for extending the normalization window, default = 5 
max_extension=5	

#Arguments Base Plot
#1=plot 0= no plotting
base_plot_flag=0
#plot dir only required if flagged
base_plot_dir=${OUTPUT_DIR}
#bandwidth to smooth base plot
base_bandwidth=3

#############
### INPUT ###
#############
#input file with sequences (first sequence is considered as the reference sequence,; relative to which all dmage is calculated later
sequence_list=${OUTPUT_DIR}/seqs_input

#TEMPORARY WRITTEN FILES'
#temporary sequence list per line
temp_sequence_list=${OUTPUT_DIR}/temp_seqs_input
#file to write temporary table with fsr calculations to
temp_fsr_file=${OUTPUT_DIR}/multiseq_fsr_table.txt

#Specified output#
#output file listing the highest dmg and higes dmg causing variant per line
wrapper_dmg_outlist=${OUTPUT_DIR}/summary_dmg_table.txt

#length of sequences to split into
kl=6

#flat stored count files
infile_plus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_plus.txt
infile_minus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_minus.txt
#naked counts background
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_h_ery_1_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_h_ery_1_minus_merged
#base composition background
infile_base_comp=${BASE_DIR}/kmer_${kl}_hg19_ploidy_removed_basecomp_plus_merged


#init. total running id
totalid=0

#init wrapper report file
echo -ne "tot.ID\tRef.seq\tVar.seq\tRef.kmer\tVar.kmer\tRef.FSR\tVar.FSR\ttotal.dmg\n" >${wrapper_dmg_outlist}

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

echo "Processing ... $sequence"

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
#same but for the naked dnase digested fibroblast background cutting
profile_naked_plus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_naked_plus} ${COMMON_FUNCTIONS}`
#split into two variables
profile_naked_plus_count=`echo ${profile_naked_plus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_naked_plus=`echo ${profile_naked_plus_out} | perl -ne '/:::(.+)/; print $1;'`
#same but for the naked dnase digested fibroblast background cutting
profile_naked_minus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_naked_minus} ${COMMON_FUNCTIONS}`
#split into two variables
profile_naked_minus_count=`echo ${profile_naked_minus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_naked_minus=`echo ${profile_naked_minus_out} | perl -ne '/:::(.+)/; print $1;'`
# 3 dnase cutting msq iqr  #
dnase_out=`Rscript ${SCRIPT_DIR}/dnase_background_cutting.R ${kmer} ${profile_naked_plus} ${dnase_iqr_fraction} ${msq_scale_thresh} ${max_extension} ${COMMON_FUNCTIONS}`
norm_skip_flag=`echo ${dnase_out} | perl -ne '/skip_flag=(\d+\.?\d*)\s+/; print $1;'`
norm_extension=`echo ${dnase_out} | perl -ne '/extension=(\d+\.?\d*)/; print $1;'`

#uses stored profiles and some vlaues for the normalization from hsc step
fsr_out=`Rscript ${SCRIPT_DIR}/fsr_calculation_merge.R ${kmer} ${profile_plus} ${profile_minus} ${profile_naked_plus} ${profile_naked_minus} ${norm_skip_flag} ${norm_extension} ${COMMON_FUNCTIONS}`

#FSR scores
fsr_merged=`echo $fsr_out} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`
fsr_plus=`echo $fsr_out} | perl -ne '/fsr_plus=([^\s]+)\s+/; print $1;'`	
fsr_minus=`echo ${fsr_out} | perl -ne '/fsr_minus=([^\s]+)\s+/; print $1;'`

#Assemble fsr values
echo -ne "\t${fsr_merged}" >>${temp_fsr_file}

done
#newline
echo "" >>${temp_fsr_file}

done

#run assessing dmg procedure
dmgout=`Rscript ${SCRIPT_DIR}/multi_seq_dmg_v3.R ${temp_fsr_file} ${kl} ${totalid}`
echo -ne "${totalid}\t$dmgout\n" >>${wrapper_dmg_outlist}

done
#order summary file by totaldmg
Rscript ${SCRIPT_DIR}/w4_order_summary.R ${wrapper_dmg_outlist} ${wrapper_dmg_outlist}

#remove temp files
rm -f ${temp_sequence_list}
rm -f ${temp_fsr_file}


