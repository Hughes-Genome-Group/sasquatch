#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -e /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -o /t1-data1/WTSA_Dev/rschwess/clustereo

SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts			
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/idx_duke_testout/merge_strands/ENCODE_H1_hESC_duke_merged/smoothed

#Specify the directory where the background files (naked dnaseI fibroblast cutting) are located,
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_h_ery_1/counts

#Specify the location of the tissue specific data (for example human erythroid is located here):
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human_ploidy_correct

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
#for erythroid human example
TAG="ENCODE_H1_hESC_duke_merged"

#select a background/normalization type (choose accordingly normalized kmer count files) there will pnorm dnase and pnorm atac files
#we have to check how much the different backgrounds differ to tell how much different norm types we need to have but at least one for atac and one for dnase will be needed (currently only dnase available)
# currently available "laza" (human fibroblast) "JH40" (human erythroid 40% mapped)
NORM_TYPE="h_ery_1"	

###INPUT###
#enter a k-mer of interest
kmer="CACGTG"

for kmer in CACGTG TGACGTC WGATAA CANNTG CMGCCC
do

mkdir -p $OUTPUT_DIR

#get & store kmer length
kl=`expr length $kmer`	

#Select accordingly normalized input files
infile_plus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_pnorm_${NORM_TYPE}_plus.txt
infile_minus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_pnorm_${NORM_TYPE}_minus.txt

#defining full paths to naked background count files
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_h_ery_1_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_h_ery_1_minus_merged

#get the plus (sense) strand profile for the kmer from the tissue specific dnase footprint counts
profile_plus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_plus} ${COMMON_FUNCTIONS}`
#output is a string like "count=10000:::14:::12:::.... 

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

#######
#merge#
#######
profile_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus} ${profile_minus}`
profile_naked_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_naked_plus} ${profile_naked_minus}`


#####################################
# 5 single plots normalized profile #
#####################################

single_plot_dir=${OUTPUT_DIR}
smooth_flag=0

#writes plot to triple_plot_dir
Rscript ${SCRIPT_DIR}/triple_plots_merging.R ${kmer} ${profile_merged} ${profile_plus} ${profile_minus} ${profile_plus_count} ${smooth_flag} ${single_plot_dir} ${COMMON_FUNCTIONS}

echo "$kmer done"

done

