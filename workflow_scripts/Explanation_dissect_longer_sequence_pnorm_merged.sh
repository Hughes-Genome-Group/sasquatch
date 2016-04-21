#!usr/bin/bash

#second functionality: split a longer sequence into kmers of a specified length adn run every kmer and assemble the stats in a table
#ideally link from every kmer latter to produce all plots as in the single run 
#but for efficiency, all plots are flagged 0 in the run through all splittet kmers

#First set directory for scripts, output directory and location of the common_functions.R file that contains 
#all shared functions sourced from some of the latter scripts:

SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts	
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/sliding_base_change/longer_CTCF_sites

# currently available "laza" (human fibroblast) "JH40" (human erythroid 40% mapped)
NORM_TYPE="h_ery_1"

ORGANISM="human"

FRAG_TYPE="DNase"

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
#for erythroid human example
TISSUE="human_erythroid_hg18"

#Arguments for the single steps as in single run 
#I'm initially have with them but it would be nice to keep them somehow easy accessible for future optimizing#

#############
### INPUT ###
#############

#enter a sequence
sequence="CCCAGGTGAGGCAGGCCAGCAAGTCCTAATTTCTGTGCTCCCGGTTCTGGCATCCAGCCAGCTTCACAAGGACCAGGGCGCTCAAGGGAAGCATTCCAAGGGGCACGGGGGAACTCTTCAAAGCCAACTACAAGTTCTGCCACCAGGGGGAGTTGGGGCCCCAGGCTGCCCTGGCACCCATCCCCTACCCCTGACTGCCACTCCTGTCTACAGGGTGTGACCGCGGGCAGGCTGCCAACTCTCCCCTTCATTGCCTACAGGTAGGGGACTGTGCCCCAGCTGTGGAGGCTGTGATGAGAA"	

#length of sequences to split into
kl=7

#file to write table to
table_file="${OUTPUT_DIR}/splitted_kmer_table_ctcf3.txt"


#define DATA directory
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/${ORGANISM}/${FRAG_TYPE}/

### BACKGROUNDs ###
case "${ORGANISM}" in 

human)

case "${FRAG_TYPE}" in 

DNase)	
#Specify the directory where the background files (naked dnaseI fibroblast cutting) are located,
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_h_ery_1/counts
#Select accordingly normalized input files
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_${NORM_TYPE}_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_${NORM_TYPE}_minus.txt
;;
ATAC)
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH_atac/counts
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_minus.txt
;;
esac

;;

mouse)
case "${FRAG_TYPE}" in

DNase)	
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid/counts
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_mm9_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_mm9_minus.txt
;;
ATAC)
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid_atac/counts
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_minus.txt
;;
esac
;;

esac

###########################
# Split Sequence in kmers #
###########################

#initialize id
id=0

#print table column names
echo -e "id\tkmer\tfsr.merged\tfsr_plus\tfsr_minus" >${table_file}


######################################################
# Process each Kmer plot free to get stats (noplots) #
######################################################

for kmer in `Rscript ${SCRIPT_DIR}/dissect_sequence.R ${sequence} ${kl} ${COMMON_FUNCTIONS}`
do

#count up kmer id 
id=$(($id + 1))

###############################
# 1 retrieve & store profiles #
###############################

#get the plus (sense) strand profile for the kmer from the tissue specific dnase footprint counts
profile_plus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_plus} ${COMMON_FUNCTIONS}`
#output is a string like "count=10000:::14:::12:::.... 

#I then split the string into two variables the counts and the profile(everything after counts=10000, seperated by :::)
#to avoid temporary files, the profiles are splitted at each ::: in each following Rscript that queries them
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

#######################
# 2 Borders SFR score #
#######################

#uses stored profiles and some vlaues for the normalization from hsc step
fsr_out=`Rscript ${SCRIPT_DIR}/fsr_calculation_pnorm_merge.R ${kmer} ${profile_merged} ${profile_plus} ${profile_minus} ${COMMON_FUNCTIONS}`
#output is again a string like: fsr_plus=x fsr_minus=x scale_plus=x scale_minus=x
fsr_merged=`echo $fsr_out} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`	
fsr_plus=`echo $fsr_out} | perl -ne '/fsr_plus=(\d+\.?\d*)\s+/; print $1;'`	
fsr_minus=`echo ${fsr_out} | perl -ne '/fsr_minus=(\d+\.?\d*)/; print $1;'`

#####################################################
# write new line for one kmer to the assembly table #
#####################################################
# fsr.plus fsr.minus scale.plus scale.minus hscor dnase.msq dnase.iqr base.msq base.iqr

echo -e "${kmer}\t${id}\t${fsr_merged}\t${fsr_plus}\t${fsr_minus}" >>${table_file}

echo "kmer processed... ${id}"

done


