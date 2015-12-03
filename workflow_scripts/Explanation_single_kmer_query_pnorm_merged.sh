#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -e /hts/data4/rschwess/clustereo
#$ -o /hts/data4/rschwess/clustereo

SCRIPT_DIR=/hts/data4/rschwess/Sasquatch_offline/Sasquatch/scripts		
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/hts/data4/rschwess/dnase_motif_tissue/idx_duke_testout/

ORGANISM="human"

FRAG_TYPE="DNase"

#I use the following to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
TISSUE="ENCODE_CD34plus_mobilized"

# currently available "laza" (human fibroblast) "JH60" (human erythroid 60% mapped)  "atac" for atac;  mm9 --> atac= atac_mm9
NORM_TYPE="JH60"

###INPUT###
#enter a k-mer of interest
kmer="GTTTGA"

#get & store kmer length
kl=`expr length $kmer`	

##############################################################
# DATA & BACKGROUNDs select according to organism and tissue #
##############################################################

#define DATA directory
DATA_DIR=/hts/data4/rschwess/database_assembly/idx_correct_assembly/${ORGANISM}/${FRAG_TYPE}/

### BACKGROUNDs ###
case "${ORGANISM}" in 

human)

case "${FRAG_TYPE}" in 

DNase)	
#Specify the directory where the background files (naked dnaseI fibroblast cutting) are located,
BACKGROUND_DIR=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH60/counts
#defining full paths to naked background count files
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_JH60_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_JH60_minus_merged
#Select accordingly normalized input files
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_${NORM_TYPE}_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_${NORM_TYPE}_minus.txt
;;
ATAC)
BACKGROUND_DIR=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH_atac/counts
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
BACKGROUND_DIR=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid/counts
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_mm9_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_mm9_minus_merged
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_mm9_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_mm9_minus.txt
;;
ATAC)
BACKGROUND_DIR=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid_atac/counts
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_mm9_atac_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_mm9_atac_minus_merged
infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_plus.txt
infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_minus.txt
;;
esac
;;

esac

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

#################
#merge profiles #
#################
profile_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus} ${profile_minus} ${FRAG_TYPE} ${kl}`
#ATAC to select the ordinary merging
profile_naked_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_naked_plus} ${profile_naked_minus} ${FRAG_TYPE} ${kl}`	

#####################
# 3 calc SFR(FSR)   #
#####################
fsr_out=`Rscript ${SCRIPT_DIR}/fsr_calculation_pnorm_merge.R ${kmer} ${profile_merged} ${profile_plus} ${profile_minus} ${COMMON_FUNCTIONS}`
#output is again a string like: fsr_plus=x fsr_minus=x scale_plus=x scale_minus=x
fsr_merged=`echo ${fsr_out} | ${PERL} -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`

fsr_plus=`echo ${fsr_out} | ${PERL} -ne '/fsr_plus=(\d+\.?\d*)\s+/; print $1;'`	
fsr_minus=`echo ${fsr_out} | ${PERL} -ne '/fsr_minus=(\d+\.?\d*)/; print $1;'`

#####################################
# 3 single plots normalized profile #
#####################################

single_plot_dir=${OUTPUT_DIR}
smooth_flag=0

#single plots norm profile
Rscript ${SCRIPT_DIR}/pnorm_single_plots_merge.R ${kmer} ${profile_merged} ${profile_plus_count} ${smooth_flag} ${OUTPUT_DIR} ${COMMON_FUNCTIONS}


#single plots norm profile
#Rscript ${SCRIPT_DIR}/pnorm_single_plots.R ${kmer} ${profile_plus} ${profile_minus} ${profile_plus_count} ${profile_minus_count} ${smooth_flag} ${OUTPUT_DIR} ${COMMON_FUNCTIONS}


#naked background cut data plot
#Rscript ${SCRIPT_DIR}/single_plots_additional/naked_cut_plot_single.R ${kmer} ${profile_naked_plus} ${profile_naked_minus} ${profile_naked_plus_count} ${profile_naked_minus_count} ${smooth_flag} ${single_plot_dir} ${COMMON_FUNCTIONS}


#####################################################
# (4) Write single Kmer results into a single table #
#####################################################

#header
echo -e "sfr.merged\tsfr.plus\tsfr.minus"  >${OUTPUT_DIR}/single_stats_table
#data
echo -e "${fsr_merged}\t${fsr_plus}\t${fsr_minus}" >>${OUTPUT_DIR}/single_stats_table


###################################################################
# (5) calculate heavy smooth correlation and save a plot if flagged #
###################################################################

#Arguments / I'm initially have with them but it would be nice to keep them somehow easy accessible for future optimizing#

#bandwidth to use for gaussian smoothing (default = 10) 
#bandwidth_hsc=10

##To plot flag: 1=plot 0=no plotting
#hsc_plot_flag=1

##directory to write the plot if flagged, only required if flagged
#hsc_plot_dir=${OUTPUT_DIR}

##execute hsc correlation and plotting
#hscor=`Rscript ${SCRIPT_DIR}/heavy_smooth_correlation.R $kmer ${bandwidth_hsc} ${profile_merged} ${profile_naked_merged} ${hsc_plot_flag} ${hsc_plot_dir}`
##returns a value into hscor and writes a plot "heavy_smooth_correlation_plot_WGATAA.png" in the plot directory if flagged



