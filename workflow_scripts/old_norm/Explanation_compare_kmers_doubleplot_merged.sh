### Explanation Example File for Comparing two kmers in a double plot fashion ###

#First set directory for scripts, output directory and location of the common_functions.R file that contains 
#all shared functions sourced from some of the latter scripts:
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts			
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/idx_duke_testout/merge_strands

#Specify the directory where the background files (base frequency and naked dnaseI fibroblast cutting) are located,
#Currently this links to files for chr1 only but I'm close to finish the complete genome as background.
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH60/counts
BASE_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/background/hg19_ploidy_removed/basecomp

#Specify the location of the tissue specific data (for example human erythroid is located here):
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human_ploidy_correct

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data
TAG="human_erythroid_hg18"

base_plot_flag=0
base_plot_dir=${OUTPUT_DIR}

###INPUT###
#enter a k-mer of interest
kmer1="GAATC"	
kmer2="GTTTC"

#get & store kmer length
kl1=`expr length $kmer1`
kl2=`expr length $kmer2`	

if [ $kl1 = $kl2 ] 
then
kl=$kl1
else
echo "Kmers dont have the same length"
fi

#based on the DATA_DIR and the TAG, precise paths to the tissue specific count files are defined and stored in variables
infile_plus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_plus.txt
infile_minus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_minus.txt

#defining full paths to naked background count files
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_JH60_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_JH60_minus_merged
#defining full paths to background base frequency count files
infile_base_comp=${BASE_DIR}/kmer_${kl}_hg19_ploidy_removed_basecomp_plus_merged

###############################
# 1 retrieve & store profiles #
###############################
#get the plus (sense) strand profile for the kmer from the tissue specific dnase footprint counts
profile_plus_out1=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer1} ${infile_plus} ${COMMON_FUNCTIONS}`
profile_plus_out2=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer2} ${infile_plus} ${COMMON_FUNCTIONS}`
#output is a string like "count=10000:::14:::12:::.... 

#I then split the string into two variables the counts and the profile(everything after counts=10000, seperated by :::)
#to avoid temporary files, the profiles are splitted at each ::: in each following Rscript that queries them
profile_plus_count1=`echo ${profile_plus_out1} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_plus1=`echo ${profile_plus_out1} | perl -ne '/:::(.+)/; print $1;'`
profile_plus_count2=`echo ${profile_plus_out2} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_plus2=`echo ${profile_plus_out2} | perl -ne '/:::(.+)/; print $1;'`

#same for antisense, minus strand, query from the tissue specific dnase footprint counts
profile_minus_out1=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer1} ${infile_minus} ${COMMON_FUNCTIONS}`
profile_minus_out2=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer2} ${infile_minus} ${COMMON_FUNCTIONS}`

#split into two variables
profile_minus_count1=`echo ${profile_minus_out1} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_minus1=`echo ${profile_minus_out1} | perl -ne '/:::(.+)/; print $1;'`
profile_minus_count2=`echo ${profile_minus_out2} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_minus2=`echo ${profile_minus_out2} | perl -ne '/:::(.+)/; print $1;'`

#same but for the naked dnase digested fibroblast background cutting
profile_naked_plus_out1=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer1} ${infile_naked_plus} ${COMMON_FUNCTIONS}`
profile_naked_plus_out2=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer2} ${infile_naked_plus} ${COMMON_FUNCTIONS}`
#split into two variables
profile_naked_plus_count1=`echo ${profile_naked_plus_out1} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_naked_plus1=`echo ${profile_naked_plus_out1} | perl -ne '/:::(.+)/; print $1;'`
profile_naked_plus_count2=`echo ${profile_naked_plus_out2} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_naked_plus2=`echo ${profile_naked_plus_out2} | perl -ne '/:::(.+)/; print $1;'`

#same but for the naked dnase digested fibroblast background cutting
profile_naked_minus_out1=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer1} ${infile_naked_minus} ${COMMON_FUNCTIONS}`
profile_naked_minus_out2=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer2} ${infile_naked_minus} ${COMMON_FUNCTIONS}`
#split into two variables
profile_naked_minus_count1=`echo ${profile_naked_minus_out1} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_naked_minus1=`echo ${profile_naked_minus_out1} | perl -ne '/:::(.+)/; print $1;'`
profile_naked_minus_count2=`echo ${profile_naked_minus_out2} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_naked_minus2=`echo ${profile_naked_minus_out2} | perl -ne '/:::(.+)/; print $1;'`
#merged#
profile_merged1=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus1} ${profile_minus1}`
profile_naked_merged1=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_naked_plus1} ${profile_naked_minus1}`
profile_merged2=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus2} ${profile_minus2}`
profile_naked_merged2=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_naked_plus2} ${profile_naked_minus2}`

###################################################################
# 2 calculate heavy smooth correlation and save a plot if flagged #
###################################################################

#bandwidth to use for gaussian smoothing (default = 10) 
bandwidth_hsc=10

#directory to write the plot if flagged, only required if flagged
hsc_plot_dir=${OUTPUT_DIR}

#STart with calculation hsc for both kmers seperately for later puposes
hsc_plot_flag=0
#execute hsc correlation and plotting
hscor1=`Rscript ${SCRIPT_DIR}/heavy_smooth_correlation.R $kmer1 ${bandwidth_hsc} ${profile_merged1} ${profile_naked_merged1} ${hsc_plot_flag} ${hsc_plot_dir}`
hscor2=`Rscript ${SCRIPT_DIR}/heavy_smooth_correlation.R $kmer2 ${bandwidth_hsc} ${profile_merged2} ${profile_naked_merged2} ${hsc_plot_flag} ${hsc_plot_dir}`

#THEN: Caluclate (and plot) the HSCoreelaition of the 2 kmers against each other
hsc_plot_flag=1
hscoragainst=`Rscript ${SCRIPT_DIR}/heavy_smooth_correlation_for_compare.R $kmer1 $kmer2 ${bandwidth_hsc} ${profile_merged1} ${profile_merged2} ${hsc_plot_flag} ${hsc_plot_dir}`


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
dnase_out=`Rscript ${SCRIPT_DIR}/dnase_background_cutting.R ${kmer1} ${profile_naked_plus1} ${dnase_iqr_fraction} ${msq_scale_thresh} ${max_extension} ${COMMON_FUNCTIONS}`
#returns a string with 4 values like: dnase_msq=x dnase_iqr=x norm_skip=x norm_extension=x"

#split up the string and save values in separate values, (probably pretty inelegant sorry)
dnase_msq1=`echo ${dnase_out} | perl -ne '/dnase_msq=(\d+\.\d+)\s+/; print $1;'`	
dnase_iqr1=`echo ${dnase_out} | perl -ne '/dnase_iqr=(\d+\.?\d*)\s+/; print $1;'`
#skip flag: 0 = normalize with default extension as usual; 1 = extend the norm scale window with extension and normalize; 2 = moderate extension still doesnt make the backgorund properly normalizable > dont normalize and flag this way
norm_skip_flag1=`echo ${dnase_out} | perl -ne '/skip_flag=(\d+\.?\d*)\s+/; print $1;'`
norm_extension1=`echo ${dnase_out} | perl -ne '/extension=(\d+\.?\d*)/; print $1;'`

#second ker
dnase_out=`Rscript ${SCRIPT_DIR}/dnase_background_cutting.R ${kmer2} ${profile_naked_plus2} ${dnase_iqr_fraction} ${msq_scale_thresh} ${max_extension} ${COMMON_FUNCTIONS}`
#returns a string with 4 values like: dnase_msq=x dnase_iqr=x norm_skip=x norm_extension=x"

#split up the string and save values in separate values, (probably pretty inelegant sorry)
dnase_msq2=`echo ${dnase_out} | perl -ne '/dnase_msq=(\d+\.\d+)\s+/; print $1;'`	
dnase_iqr2=`echo ${dnase_out} | perl -ne '/dnase_iqr=(\d+\.?\d*)\s+/; print $1;'`
#skip flag: 0 = normalize with default extension as usual; 1 = extend the norm scale window with extension and normalize; 2 = moderate extension still doesnt make the backgorund properly normalizable > dont normalize and flag this way
norm_skip_flag2=`echo ${dnase_out} | perl -ne '/skip_flag=(\d+\.?\d*)\s+/; print $1;'`
norm_extension2=`echo ${dnase_out} | perl -ne '/extension=(\d+\.?\d*)/; print $1;'`

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
base_out=`Rscript ${SCRIPT_DIR}/base_composition.R ${kmer1} ${infile_base_comp} ${base_plot_flag} ${base_plot_dir} ${base_bandwidth} ${COMMON_FUNCTIONS}`
#split and save values in separate values
base_msq1=`echo ${base_out} | perl -ne '/base_msq=(\d+\.\d+)\s+/; print $1;'`	
base_iqr1=`echo ${base_out} | perl -ne '/base_iqr=(\d+\.?\d*)/; print $1;'`
#second kmer
base_out=`Rscript ${SCRIPT_DIR}/base_composition.R ${kmer2} ${infile_base_comp} ${base_plot_flag} ${base_plot_dir} ${base_bandwidth} ${COMMON_FUNCTIONS}`
#split and save values in separate values
base_msq2=`echo ${base_out} | perl -ne '/base_msq=(\d+\.\d+)\s+/; print $1;'`	
base_iqr2=`echo ${base_out} | perl -ne '/base_iqr=(\d+\.?\d*)/; print $1;'`

###########################################################################################
# 5 normalize  and smooth profile with gaussian, estimate borders and calculate FSR score #
###########################################################################################
#uses stored profiles and some vlaues for the normalization from hsc step
fsr_out1=`Rscript ${SCRIPT_DIR}/fsr_calculation_merge.R ${kmer} ${profile_plus1} ${profile_minus1} ${profile_naked_plus1} ${profile_naked_minus1} ${norm_skip_flag} ${norm_extension} ${COMMON_FUNCTIONS}`
fsr_out2=`Rscript ${SCRIPT_DIR}/fsr_calculation_merge.R ${kmer} ${profile_plus2} ${profile_minus2} ${profile_naked_plus2} ${profile_naked_minus2} ${norm_skip_flag} ${norm_extension} ${COMMON_FUNCTIONS}`

#output is again a string like: fsr_plus=x fsr_minus=x scale_plus=x scale_minus=x
#FSRatios
fsr_plus1=`echo $fsr_out1} | perl -ne '/fsr_plus=(\d+\.?\d*)\s+/; print $1;'`	
fsr_minus1=`echo ${fsr_out1} | perl -ne '/fsr_minus=(\d+\.?\d*)/; print $1;'`
fsr_merged1=`echo $fsr_out1} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`
#Scale factors
scale_plus1=`echo ${fsr_out1} | perl -ne '/scale_plus=(\d+\.?\d*)/; print $1;'`
scale_minus1=`echo ${fsr_out1} | perl -ne '/scale_minus=(\d+\.?\d*)/; print $1;'`
#FSRatios
fsr_plus2=`echo $fsr_out2} | perl -ne '/fsr_plus=(\d+\.?\d*)\s+/; print $1;'`	
fsr_minus2=`echo ${fsr_out2} | perl -ne '/fsr_minus=(\d+\.?\d*)/; print $1;'`
fsr_merged2=`echo $fsr_out2} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`
#Scale factors
scale_plus2=`echo ${fsr_out2} | perl -ne '/scale_plus=(\d+\.?\d*)/; print $1;'`
scale_minus2=`echo ${fsr_out2} | perl -ne '/scale_minus=(\d+\.?\d*)/; print $1;'`

#################################################################################################
# At this point all necessary stats are available and the rest are some more plotting functions #
#################################################################################################
#writing a single stats table file
#header
echo -e "fsr.merged\tfsr.plus\tfsr.minus\tscale.plus\tscale.minus\thscor\tdnase.msq\tdnase.iqr\tbase.msq\tbase.iqr" >${OUTPUT_DIR}/single_stats_table_${kmer1}
echo -e "fsr.merged\tfsr.plus\tfsr.minus\tscale.plus\tscale.minus\thscor\tdnase.msq\tdnase.iqr\tbase.msq\tbase.iqr" >${OUTPUT_DIR}/single_stats_table_${kmer2}
#data
echo -e "${fsr_merged1}\t${fsr_plus1}\t${fs_minus1}\t${scale_plus1}\t${scale_minus1}\t${hscor1}\t${dnase_msq1}\t${dnase_iqr1}\t${base_msq1}\t${base_iqr1}" >>${OUTPUT_DIR}/single_stats_table_${kmer1}
echo -e "${fsr_merged2}\t${fsr_plus2}\t${fs_minus2}\t${scale_plus2}\t${scale_minus2}\t${hscor2}\t${dnase_msq2}\t${dnase_iqr2}\t${base_msq2}\t${base_iqr2}" >>${OUTPUT_DIR}/single_stats_table_${kmer2}

################################################
# 6 single overlap plot of normalized profiles #
################################################
smooth_flag=1
#directory to write plot to
plot_dir=${OUTPUT_DIR}

Rscript ${SCRIPT_DIR}/overlap_plot_normalized_merge.R ${kmer1} ${kmer2} ${profile_plus1} ${profile_minus1} ${profile_naked_plus1} ${profile_naked_minus1} ${profile_plus2} ${profile_minus2} ${profile_naked_plus2} ${profile_naked_minus2} ${profile_plus_count1} ${profile_minus_count1} ${profile_plus_count2} ${profile_minus_count2} ${smooth_flag} ${norm_skip_flag1} ${norm_skip_flag2} ${norm_extension1} ${norm_extension2} ${plot_dir} ${COMMON_FUNCTIONS}

### OPTIONAL ###
########################################################
# 7 get triple plot: real data, background, normalized #
########################################################
#standrad triple plot for each kmer /// OPTIONAL
#Arguments choose if smoothed raw profiles should be printed or not default = 1 > smooth raw data for plots
#1#
smooth_flag=1
#directory to write plot to
triple_plot_dir=${OUTPUT_DIR}

#writes plot to triple_plot_dir
Rscript ${SCRIPT_DIR}/triple_plots.R ${kmer1} ${profile_plus1} ${profile_minus1} ${profile_naked_plus1} ${profile_naked_minus1} ${profile_plus_count1} ${profile_minus_count1} ${profile_naked_plus_count1} ${profile_naked_minus_count1} ${smooth_flag} ${norm_skip_flag1} ${norm_extension1} ${triple_plot_dir} ${COMMON_FUNCTIONS}

Rscript ${SCRIPT_DIR}/triple_plots.R ${kmer2} ${profile_plus2} ${profile_minus2} ${profile_naked_plus2} ${profile_naked_minus2} ${profile_plus_count2} ${profile_minus_count2} ${profile_naked_plus_count2} ${profile_naked_minus_count2} ${smooth_flag} ${norm_skip_flag2} ${norm_extension2} ${triple_plot_dir} ${COMMON_FUNCTIONS}


