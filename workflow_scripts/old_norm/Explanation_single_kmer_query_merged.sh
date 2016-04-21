### Explanation Example File for Single kmer Query Workflow ###

#The first functionality is to query a single 5, 6, 7 - mer from a specified tissue database 

#First set directory for scripts, output directory and location of the common_functions.R file that contains 
#all shared functions sourced from some of the latter scripts:
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts			
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/idx_duke_testout/merge_strands

#Specify the directory where the background files (base frequency and naked dnaseI fibroblast cutting) are located,
#Currently this links to files for chr1 only but I'm close to finish the complete genome as background.
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_h_ery_1/counts
BASE_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/background/hg19_ploidy_removed/basecomp

#Specify the location of the tissue specific data (for example human erythroid is located here):
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human_ploidy_correct

#DNase or ATAC
FRAG_TYPE="DNase"

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
#for erythroid human example
TAG="human_erythroid_hg18"

###INPUT###
#enter a k-mer of interest
kmer="WGATAA"

#get & store kmer length

kl=`expr length $kmer`	

#based on the DATA_DIR and the TAG, precise paths to the tissue specific count files are defined and stored in variables

infile_plus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_plus.txt
infile_minus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_minus.txt

#defining full paths to naked background count files
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_h_ery_1_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_h_ery_1_minus_merged

#defining full paths to background base frequency count files
infile_base_comp=${BASE_DIR}/kmer_${kl}_hg19_ploidy_removed_basecomp_plus_merged


###############################
# 1 retrieve & store profiles #
###############################

#currently I have multiple outputs from one Rscript assembled in a easy to split string and do this right after each script. 
#This is surely inelegant but I didn't come up with something better so far. If you need it differently implemented let me know

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
#merged#
profile_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus} ${profile_minus} ${FRAG_TYPE} ${kl}`
profile_naked_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_naked_plus} ${profile_naked_minus}`

###################################################################
# 2 calculate heavy smooth correlation and save a plot if flagged #
###################################################################

#Arguments / I'm initially have with them but it would be nice to keep them somehow easy accessible for future optimizing#

#bandwidth to use for gaussian smoothing (default = 10) 
bandwidth_hsc=10

#To plot flag: 1=plot 0=no plotting
hsc_plot_flag=1

#directory to write the plot if flagged, only required if flagged
hsc_plot_dir=${OUTPUT_DIR}

#execute hsc correlation and plotting
hscor=`Rscript ${SCRIPT_DIR}/heavy_smooth_correlation.R $kmer ${bandwidth_hsc} ${profile_merged} ${profile_naked_merged} ${hsc_plot_flag} ${hsc_plot_dir}`
#returns a value into hscor and writes a plot "heavy_smooth_correlation_plot_WGATAA.png" in the plot directory if flagged


#########################################################################
# 3 dnase cutting msq iqr  #
#########################################################################

#Arguments
#fraction to calculate iqr (default 4 = 25% IQR)
dnase_iqr_fraction=4	
#MSQ threshold to decide if to normalize, if to extend the norm window or if to flag as not normalized, default=0.09 but its empirically adjusted to 1-2 datasets so should somehow be adjustable later on
msq_scale_thresh=0.09
#maximal extension allowed for extending the normalization window, default = 5 
max_extension=5	

#execute dnase cut background MSQ and IQR calculation and retrieve normalization to use for later as norm.flag and extension to use
dnase_out=`Rscript ${SCRIPT_DIR}/dnase_background_cutting.R ${kmer} ${profile_naked_plus} ${dnase_iqr_fraction} ${msq_scale_thresh} ${max_extension} ${COMMON_FUNCTIONS}`
#returns a string with 4 values like: dnase_msq=x dnase_iqr=x norm_skip=x norm_extension=x"

#split up the string and save values in separate values, (probably pretty inelegant sorry)
dnase_msq=`echo ${dnase_out} | perl -ne '/dnase_msq=(\d+\.\d+)\s+/; print $1;'`	
dnase_iqr=`echo ${dnase_out} | perl -ne '/dnase_iqr=(\d+\.?\d*)\s+/; print $1;'`
#skip flag: 0 = normalize with default extension as usual; 1 = extend the norm scale window with extension and normalize; 2 = moderate extension still doesnt make the backgorund properly normalizable > dont normalize and flag this way
norm_skip_flag=`echo ${dnase_out} | perl -ne '/skip_flag=(\d+\.?\d*)\s+/; print $1;'`
norm_extension=`echo ${dnase_out} | perl -ne '/extension=(\d+\.?\d*)/; print $1;'`

###############################################################################
# 4 get base composition around kmer, store base profile and plot if flagged  #
###############################################################################

#Arguments
#Plot flag: 1=plot 0= no plotting
base_plot_flag=1
#plot dir only required if flagged
base_plot_dir=${OUTPUT_DIR}
#bandwidth to smooth base plot
base_bandwidth=3

#execute msq and iqr calculation for base composition and plotting if flagged
base_out=`Rscript ${SCRIPT_DIR}/base_composition.R ${kmer} ${infile_base_comp} ${base_plot_flag} ${base_plot_dir} ${base_bandwidth} ${COMMON_FUNCTIONS}`
#output is string like: base_msq=x base_iqr=x, if flagged writes plot: "base_composition_plot_WGATAA.png"

#split and save values in separate values
base_msq=`echo ${base_out} | perl -ne '/base_msq=(\d+\.\d+)\s+/; print $1;'`	
base_iqr=`echo ${base_out} | perl -ne '/base_iqr=(\d+\.?\d*)/; print $1;'`

###########################################################################################
# 5 normalize  and smooth profile with gaussian, estimate borders and calculate FSR score #
###########################################################################################
#uses stored profiles and some vlaues for the normalization from hsc step
fsr_out=`Rscript ${SCRIPT_DIR}/fsr_calculation_merge.R ${kmer} ${profile_plus} ${profile_minus} ${profile_naked_plus} ${profile_naked_minus} ${norm_skip_flag} ${norm_extension} ${COMMON_FUNCTIONS}`

#FSR scores
fsr_merged=`echo $fsr_out} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`
fsr_plus=`echo $fsr_out} | perl -ne '/fsr_plus=([^\s]+)\s+/; print $1;'`	
fsr_minus=`echo ${fsr_out} | perl -ne '/fsr_minus=([^\s]+)\s+/; print $1;'`

profile_merge=`echo ${fsr_out} | perl -ne '/profile=:(.+)/; print $1;'`

#Scale factors
scale_plus=`echo ${fsr_out} | perl -ne '/scale_plus=(\d+\.?\d*)\s+/; print $1;'`
scale_minus=`echo ${fsr_out} | perl -ne '/scale_minus=(\d+\.?\d*)\s+/; print $1;'`

#################################################################################################
# At this point all necessary stats are available and the rest are some more plotting functions #
#################################################################################################
#writing a single stats table file

#header
echo -e "fsr.merged\tfsr.plus\tfsr.minus\tscale.plus\tscale.minus\thscor\tdnase.msq\tdnase.iqr\tbase.msq\tbase.iqr" >${OUTPUT_DIR}/single_stats_table
#data
echo -e "${fsr_merged}\t${fsr_plus}\t${fs_minus}\t${scale_plus}\t${scale_minus}\t${hscor}\t${dnase_msq}\t${dnase_iqr}\t${base_msq}\t${base_iqr}" >>${OUTPUT_DIR}/single_stats_table

########################################################
# 6 get triple plot: real data, background, normalized #
########################################################
#Arguments choose if smoothed raw profiles should be printed or not default = 1 > smooth raw data for plots
smooth_flag=0
#directory to write plot to
triple_plot_dir=${OUTPUT_DIR}

#writes plot to triple_plot_dir
Rscript ${SCRIPT_DIR}/triple_plots.R ${kmer} ${profile_plus} ${profile_minus} ${profile_naked_plus} ${profile_naked_minus} ${profile_plus_count} ${profile_minus_count} ${profile_naked_plus_count} ${profile_naked_minus_count} ${smooth_flag} ${norm_skip_flag} ${norm_extension} ${triple_plot_dir} ${COMMON_FUNCTIONS}


##########################
# (7) single plot merged #
##########################

smooth_flag=0

#single plots norm profile
Rscript ${SCRIPT_DIR}/pnorm_single_plots_merge.R ${kmer} ${profile_merge} ${profile_plus_count} ${smooth_flag} ${OUTPUT_DIR} ${COMMON_FUNCTIONS}

