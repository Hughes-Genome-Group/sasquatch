#!usr/bin/bash

#third funcionality: Dissect two list (of identifcal length preferred) and perform functionality 2 for both lists. Report a table with stats for each sequence and produce a comparison table by calculating the differences of the FSRatios 

#First set directory for scripts, output directory and location of the common_functions.R file that contains 
#all shared functions sourced from some of the latter scripts:

SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts		
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/idx_duke_testout

ORGANISM="human"

FRAG_TYPE="DNase"

#select a background/normalization type (choose accordingly normalized kmer count files) there will pnorm dnase and pnorm atac files
#we have to check how much the different backgrounds differ to tell how much different norm types we need to have but at least one for atac and one for dnase will be needed (currently only dnase available)
# currently available "laza" (human fibroblast) "JH40" (human erythroid 40% mapped)
NORM_TYPE="h_ery_1"

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
#for erythroid human example
TISSUE="ENCODE_HepG2_duke_merged"

#Specify outputfile names and location
table_file1=${OUTPUT_DIR}/seq1_table.txt
table_file2=${OUTPUT_DIR}/seq2_table.txt
compare_file=${OUTPUT_DIR}/compare_table.txt

#########
# INPUT #
#########

#enter Sequences
sequence1="GCCTATTCAACAC"
sequence2="GCCTATCCAACAC"

#length of sequences to split into
kl=6

##############################################################
# DATA & BACKGROUNDs select according to organism and tissue #
##############################################################

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
infile_plus=${DATA_DIR}/${TISSUE}/kmers_${kl}_count_${TISSUE}_pnorm_mm9_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/kmers_${kl}_count_${TISSUE}_pnorm_mm9_minus.txt
;;
ATAC)
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid_atac/counts
infile_plus=${DATA_DIR}/${TISSUE}/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_minus.txt
;;
esac
;;

esac

###########################
# Split Sequence in kmers #
###########################

#### First Sequence ###

#initialize id
id=0

#print table column names
echo -e "kmer\tid\tfsr_merged\tfsr_plus\tfsr_minus" >${table_file1}

#split sequence1 and perform stats caalculation for each kmer, writing everything into table1
for kmer in `Rscript ${SCRIPT_DIR}/dissect_sequence.R ${sequence1} ${kl} ${COMMON_FUNCTIONS}`
do

######################################################
# Process each Kmer plot free to get stats (noplots) #
######################################################

id=$(($id + 1))

###############################
# 1 retrieve & store profiles #
###############################

#currently i have multiple outputs from one Rscript asembled in a easy to split string and do this right after each script. This is surely inelegant but i didnt come up with something better so far. if you need it differently implemented let me know
profile_plus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_plus} ${COMMON_FUNCTIONS}`
profile_plus_count=`echo ${profile_plus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_plus=`echo ${profile_plus_out} | perl -ne '/:::(.+)/; print $1;'`

profile_minus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_minus} ${COMMON_FUNCTIONS}`
profile_minus_count=`echo ${profile_minus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_minus=`echo ${profile_minus_out} | perl -ne '/:::(.+)/; print $1;'`

#######
#merge#
#######
profile_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus} ${profile_minus} ${FRAG_TYPE} ${kl}`

###########################################################################################
# 2 normalize  and smooth profile with gaussian, estimate borders and calculate FSR score #
###########################################################################################
fsr_out=`Rscript ${SCRIPT_DIR}/fsr_calculation_pnorm_merge.R ${kmer} ${profile_merged} ${profile_plus} ${profile_minus} ${COMMON_FUNCTIONS}`
#output is again a string like: fsr_plus=x fsr_minus=x scale_plus=x scale_minus=x
fsr_merged=`echo $fsr_out} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`	
fsr_plus=`echo $fsr_out} | perl -ne '/fsr_plus=(\d+\.?\d*)\s+/; print $1;'`	
fsr_minus=`echo ${fsr_out} | perl -ne '/fsr_minus=(\d+\.?\d*)/; print $1;'`

#################################################################################################
# At this point all necessary stats are available and the rest are some more plotting functions #
#################################################################################################

# fsr.plus fsr.minus
echo -e "${kmer}\t${id}\t${fsr_merged}\t${fsr_plus}\t${fsr_minus}" >>${table_file1}

echo "kmer processed... ${id}"

done

echo "Sequence 1 complete"



### Second Sequence ###

#initialize id
id=0

#print table column names
echo -e "kmer\tid\tfsr_merged\tfsr_plus\tfsr_minus" >${table_file2}

#split sequence 2 into kmers and process each wwriting results up in table2
for kmer in `Rscript ${SCRIPT_DIR}/dissect_sequence.R ${sequence2} ${kl} ${COMMON_FUNCTIONS}`
do

######################################################
# Process each Kmer plot free to get stats (noplots) #
######################################################

id=$(($id + 1))

###############################
# 1 retrieve & store profiles #
###############################

#currently i have multiple outputs from one Rscript asembled in a easy to split string and do this right after each script. This is surely inelegant but i didnt come up with something better so far. if you need it differently implemented let me know
profile_plus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_plus} ${COMMON_FUNCTIONS}`
profile_plus_count=`echo ${profile_plus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_plus=`echo ${profile_plus_out} | perl -ne '/:::(.+)/; print $1;'`

profile_minus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_minus} ${COMMON_FUNCTIONS}`
profile_minus_count=`echo ${profile_minus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_minus=`echo ${profile_minus_out} | perl -ne '/:::(.+)/; print $1;'`

#######
#merge#
#######
profile_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus} ${profile_minus} ${FRAG_TYPE} ${kl}`
###########################################################################################
# 2 normalize  and smooth profile with gaussian, estimate borders and calculate FSR score #
###########################################################################################
fsr_out=`Rscript ${SCRIPT_DIR}/fsr_calculation_pnorm_merge.R ${kmer} ${profile_merged} ${profile_plus} ${profile_minus} ${COMMON_FUNCTIONS}`
#output is again a string like: fsr_plus=x fsr_minus=x scale_plus=x scale_minus=x
fsr_merged=`echo $fsr_out} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`	
fsr_plus=`echo $fsr_out} | perl -ne '/fsr_plus=(\d+\.?\d*)\s+/; print $1;'`	
fsr_minus=`echo ${fsr_out} | perl -ne '/fsr_minus=(\d+\.?\d*)/; print $1;'`

#################################################################################################
# At this point all necessary stats are available and the rest are some more plotting functions #
#################################################################################################

# fsr.plus fsr.minus
echo -e "${kmer}\t${id}\t${fsr_merged}\t${fsr_plus}\t${fsr_minus}" >>${table_file2}

echo "kmer processed... ${id}"

done

echo "Sequence 2 complete"


#create a comparison table 
Rscript ${SCRIPT_DIR}/compare_table_merge.R ${table_file1} ${table_file2} >${compare_file}






