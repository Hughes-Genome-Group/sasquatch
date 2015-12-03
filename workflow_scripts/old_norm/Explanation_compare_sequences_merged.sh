#!usr/bin/bash

#third funcionality: Dissect two list (of identifcal length preferred) and perform functionality 2 for both lists. Report a table with stats for each sequence and produce a comparison table by calculating the differences of the FSRatios 

#First set directory for scripts, output directory and location of the common_functions.R file that contains 
#all shared functions sourced from some of the latter scripts:
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts			
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/idx_duke_testout/merge_strands

#Specify the directory where the background files (base frequency and naked dnaseI fibroblast cutting) are located,
#Currently this links to files for chr1 only but I'm close to finish the complete genome as background.

#Specify the directory where the background files (naked dnaseI fibroblast cutting) are located,
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH60/counts
BASE_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/background/hg19_ploidy_removed/basecomp

#Specify the location of the tissue specific data (for example human erythroid is located here):
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human_ploidy_correct

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
#for erythroid human example
TAG="human_erythroid_hg18"		

#Specify outputfile names and location
table_file1=${OUTPUT_DIR}/seq1_table.txt
table_file2=${OUTPUT_DIR}/seq2_table.txt
compare_file=${OUTPUT_DIR}/compare_table.txt

#Arguments for the single steps as in single run 
#I'm initially have with them but it would be nice to keeo them somehow easy accessible for future optimizing#

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

#########
# INPUT #
#########

#enter Sequences
sequence1="CACGTGAATC"	
sequence2="CACGTGTTTC"	

#length of sequences to split into
kl=6

#flat stored count files
infile_plus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_plus.txt
infile_minus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_minus.txt
#naked counts background
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_JH60_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_JH60_minus_merged
#base composition background
infile_base_comp=${BASE_DIR}/kmer_${kl}_hg19_ploidy_removed_basecomp_plus_merged


###########################
# Split Sequence in kmers #
###########################
#### First Sequence ###

#initialize id
id=0

#print table column names
echo -e "kmer\tid\tfsr_merged\tfsr_plus\tfsr_minus\tscale_plus\tscale_minus\thscor\tdnase_msq\tdnase_iqr\tbase_msq\tbase_iqr" >${table_file1}

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

profile_naked_plus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_naked_plus} ${COMMON_FUNCTIONS}`
profile_naked_plus_count=`echo ${profile_naked_plus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_naked_plus=`echo ${profile_naked_plus_out} | perl -ne '/:::(.+)/; print $1;'`

profile_naked_minus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_naked_minus} ${COMMON_FUNCTIONS}`
profile_naked_minus_count=`echo ${profile_naked_minus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_naked_minus=`echo ${profile_naked_minus_out} | perl -ne '/:::(.+)/; print $1;'`
#merged#
profile_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus} ${profile_minus}`
profile_naked_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_naked_plus} ${profile_naked_minus}`

###################################################################
# 2 calculate heavy smooth correlation and save a plot if flagged #
###################################################################
#execute hsc correlation and plotting
hscor=`Rscript ${SCRIPT_DIR}/heavy_smooth_correlation.R $kmer ${bandwidth_hsc} ${profile_merged} ${profile_naked_merged} ${hsc_plot_flag} ${hsc_plot_dir}`

#########################################################################
# 3 dnase cutting msq iqr  #
#########################################################################

#execute dnase cut background MSQ and IQR caclulation and retrieve normalization to use for later as norm.flag and extension to use
dnase_out=`Rscript ${SCRIPT_DIR}/dnase_background_cutting.R ${kmer} ${profile_naked_plus} ${dnase_iqr_fraction} ${msq_scale_thresh} ${max_extension} ${COMMON_FUNCTIONS}`

#save values in seperate values since output is in one string (probably pretty inelegant sorry)
dnase_msq=`echo ${dnase_out} | perl -ne '/dnase_msq=([^\s]+)\s+/; print $1;'`	
dnase_iqr=`echo ${dnase_out} | perl -ne '/dnase_iqr=([^\s]+)\s+/; print $1;'`
#skip flag: 0 = normalize with default extension as usual; 1 = extend the norm sclae window with extension and normalize; 2 = moderate extension still doesnt make the backgorund properly normalizable > dont normalize and flag this way
norm_skip_flag=`echo ${dnase_out} | perl -ne '/skip_flag=([^\s]+)\s+/; print $1;'`
norm_extension=`echo ${dnase_out} | perl -ne '/extension=([^\s]+)\s+/; print $1;'`

#########################################################################
# 4 get base composition around kmer, store profile and plot i flagged  #
#########################################################################

#execute msq and iqr calculation for base composition and plotting if flagged
base_out=`Rscript ${SCRIPT_DIR}/base_composition.R ${kmer} ${infile_base_comp} ${base_plot_flag} ${base_plot_dir} ${base_bandwidth} ${COMMON_FUNCTIONS}`

#save values in seperate values since output is in one string (probably pretty inelegant sorry)
base_msq=`echo ${base_out} | perl -ne '/base_msq=([^\s]+)\s+/; print $1;'`	
base_iqr=`echo ${base_out} | perl -ne '/base_iqr=([^\s]+)\s+/; print $1;'`


###########################################################################################
# 5 normalize  and smooth profile with gaussian, estimate borders and calculate FSR score #
###########################################################################################
fsr_out=`Rscript ${SCRIPT_DIR}/fsr_calculation_merge.R ${kmer} ${profile_plus} ${profile_minus} ${profile_naked_plus} ${profile_naked_minus} ${norm_skip_flag} ${norm_extension} ${COMMON_FUNCTIONS}`

#FSR scores
fsr_merged=`echo $fsr_out} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`
fsr_plus=`echo $fsr_out} | perl -ne '/fsr_plus=([^\s]+)\s+/; print $1;'`	
fsr_minus=`echo ${fsr_out} | perl -ne '/fsr_minus=([^\s]+)\s+/; print $1;'`

#Scale factors
scale_plus=`echo ${fsr_out} | perl -ne '/scale_plus=([^\s]+)\s+/; print $1;'`
scale_minus=`echo ${fsr_out} | perl -ne '/scale_minus=([^\s]+)\s+/; print $1;'`

#################################################################################################
# At this point all necessary stats are available and the rest are some more plotting functions #
#################################################################################################

# fsr.plus fsr.minus scale.plus scale.minus hscor dnase.msq dnase.iqr base.msq base.iqr
echo -e "${kmer}\t${id}\t${fsr_merged}\t${fsr_plus}\t${fsr_minus}\t${scale_plus}\t${scale_minus}\t${hscor}\t${dnase_msq}\t${dnase_iqr}\t${base_msq}\t${base_iqr}" >>${table_file1}

echo "kmer processed... ${id}"

done

echo "Sequence 1 complete"



### Second Sequence ###

#initialize id
id=0

#print table column names
echo -e "kmer\tid\tfsr_merged\tfsr_plus\tfsr_minus\tscale_plus\tscale_minus\thscor\tdnase_msq\tdnase_iqr\tbase_msq\tbase_iqr" >${table_file2}

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

profile_naked_plus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_naked_plus} ${COMMON_FUNCTIONS}`
profile_naked_plus_count=`echo ${profile_naked_plus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_naked_plus=`echo ${profile_naked_plus_out} | perl -ne '/:::(.+)/; print $1;'`

profile_naked_minus_out=`Rscript ${SCRIPT_DIR}/retrieve_profile.R ${kmer} ${infile_naked_minus} ${COMMON_FUNCTIONS}`
profile_naked_minus_count=`echo ${profile_naked_minus_out} | perl -ne '/count=(\d+\.?\d*):::/; print $1;'`
profile_naked_minus=`echo ${profile_naked_minus_out} | perl -ne '/:::(.+)/; print $1;'`
#merged#
profile_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus} ${profile_minus}`
profile_naked_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_naked_plus} ${profile_naked_minus}`

###################################################################
# 2 calculate heavy smooth correlation and save a plot if flagged #
###################################################################
#execute hsc correlation and plotting
hscor=`Rscript ${SCRIPT_DIR}/heavy_smooth_correlation.R $kmer ${bandwidth_hsc} ${profile_merged} ${profile_naked_merged} ${hsc_plot_flag} ${hsc_plot_dir}`

#########################################################################
# 3 dnase cutting msq iqr  #
#########################################################################
#execute dnase cut background MSQ and IQR caclulation and retrieve normalization to use for later as norm.flag and extension to use
dnase_out=`Rscript ${SCRIPT_DIR}/dnase_background_cutting.R ${kmer} ${profile_naked_plus} ${dnase_iqr_fraction} ${msq_scale_thresh} ${max_extension} ${COMMON_FUNCTIONS}`
#save values in seperate values since output is in one string (probably pretty inelegant sorry)
dnase_msq=`echo ${dnase_out} | perl -ne '/dnase_msq=([^\s]+)\s+/; print $1;'`	
dnase_iqr=`echo ${dnase_out} | perl -ne '/dnase_iqr=([^\s]+)\s+/; print $1;'`
#skip flag: 0 = normalize with default extension as usual; 1 = extend the norm sclae window with extension and normalize; 2 = moderate extension still doesnt make the backgorund properly normalizable > dont normalize and flag this way
norm_skip_flag=`echo ${dnase_out} | perl -ne '/skip_flag=([^\s]+)\s+/; print $1;'`
norm_extension=`echo ${dnase_out} | perl -ne '/extension=([^\s]+)\s+/; print $1;'`

#########################################################################
# 4 get base composition around kmer, store profile and plot i flagged  #
#########################################################################
#execute msq and iqr calculation for base composition and plotting if flagged
base_out=`Rscript ${SCRIPT_DIR}/base_composition.R ${kmer} ${infile_base_comp} ${base_plot_flag} ${base_plot_dir} ${base_bandwidth} ${COMMON_FUNCTIONS}`
#save values in seperate values since output is in one string (probably pretty inelegant sorry)
base_msq=`echo ${base_out} | perl -ne '/base_msq=([^\s]+)\s+/; print $1;'`	
base_iqr=`echo ${base_out} | perl -ne '/base_iqr=([^\s]+)\s+/; print $1;'`

###########################################################################################
# 5 normalize  and smooth profile with gaussian, estimate borders and calculate FSR score #
###########################################################################################
fsr_out=`Rscript ${SCRIPT_DIR}/fsr_calculation_merge.R ${kmer} ${profile_plus} ${profile_minus} ${profile_naked_plus} ${profile_naked_minus} ${norm_skip_flag} ${norm_extension} ${COMMON_FUNCTIONS}`

#FSR scores
fsr_merged=`echo $fsr_out} | perl -ne '/fsr_merged=(\d+\.?\d*)\s+/; print $1;'`
fsr_plus=`echo $fsr_out} | perl -ne '/fsr_plus=([^\s]+)\s+/; print $1;'`	
fsr_minus=`echo ${fsr_out} | perl -ne '/fsr_minus=([^\s]+)\s+/; print $1;'`

#Scale factors
scale_plus=`echo ${fsr_out} | perl -ne '/scale_plus=([^\s]+)\s+/; print $1;'`
scale_minus=`echo ${fsr_out} | perl -ne '/scale_minus=([^\s]+)\s+/; print $1;'`
#perl -ne '/scale_minus=(\d+\.?\d*)/; print $1;'`#keep save

#################################################################################################
# At this point all necessary stats are available and the rest are some more plotting functions #
#################################################################################################

# fsr.plus fsr.minus scale.plus scale.minus hscor dnase.msq dnase.iqr base.msq base.iqr
echo -e "${kmer}\t${id}\t${fsr_merged}\t${fsr_plus}\t${fsr_minus}\t${scale_plus}\t${scale_minus}\t${hscor}\t${dnase_msq}\t${dnase_iqr}\t${base_msq}\t${base_iqr}" >>${table_file2}

echo "kmer processed... ${id}"

done

echo "Sequence 2 complete"


#create a comparison table 
Rscript ${SCRIPT_DIR}/compare_table_merge.R ${table_file1} ${table_file2} >${compare_file}






