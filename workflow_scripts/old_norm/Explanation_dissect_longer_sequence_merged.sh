#!usr/bin/bash

#second functionality: split a longer sequence into kmers of a specified length adn run every kmer and assemble the stats in a table
#ideally link from every kmer latter to produce all plots as in the single run 
#but for efficiency, all plots are flagged 0 in the run through all splittet kmers

SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts			
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/idx_duke_testout/merge_strands

#Specify the directory where the background files (base frequency and naked dnaseI fibroblast cutting) are located,
#Currently this links to files for chr1 only but I'm close to finish the complete genome as background.
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH60/counts

#Specify the location of the tissue specific data (for example human erythroid is located here):
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human_ploidy_correct

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
#for erythroid human example
TAG="human_erythroid_hg18"		

#Arguments for the single steps as in single run 
#I'm initially have with them but it would be nice to keeo them somehow easy accessible for future optimizing#

#Arguments for HSCor
#bandwidth to use for gaussian smoothing (default = 10) 
bandwidth_hsc=10
#1=plot 0= no plotting
hsc_plot_flag=0
#dir to write plot if flagged only required if flagged
hsc_plot_dir=${OUTPUT_DIR}

base_plot_flag=0
base_plot_dir=${OUTPUT_DIR}

#############
### INPUT ###
#############

#enter a sequence
sequence="CTATGCACGTG"	

#length of sequences to split into
kl=6

#flat stored count files
infile_plus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_plus.txt
infile_minus=${DATA_DIR}/${TAG}/counts/kmers_${kl}_count_${TAG}_minus.txt
#naked counts background
infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_JH60_plus_merged
infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_JH60_minus_merged

###########################
# Split Sequence in kmers #
###########################

#file to write table to
table_file="${OUTPUT_DIR}/splitted_kmer_table.txt"

#initialize id
id=0

#print table column names
echo -e "kmer\tid\tfsr_merged\tfsr_plus\tfsr_minus\tscale_plus\tscale_minus\thscor\tdnase_msq\tdnase_iqr\tbase_msq\tbase_iqr\tnorm.flag\tnorm.extension" >${table_file}


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

#currently i have multiple outputs from one Rscript asembled in a easy to split string and do this right after each script. This is surely inelegant but i didnt come up with something better so far. if you need it differently implemented let me know

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
profile_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus} ${profile_minus}`
profile_naked_merged=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_naked_plus} ${profile_naked_minus}`

###################################################################
# 2 calculate heavy smooth correlation and save a plot if flagged #
###################################################################
#execute hsc correlation and plotting
hscor=`Rscript ${SCRIPT_DIR}/heavy_smooth_correlation.R $kmer ${bandwidth_hsc} ${profile_merged} ${profile_naked_merged} ${hsc_plot_flag} ${hsc_plot_dir}`
#returns a value into hscor and writes a plot "heavy_smooth_correlation_plot_WGATAA.png" in the plot directory if flagged

#########################################################################
# 3 dnase cutting msq iqr  #
#########################################################################

#execute dnase cut background MSQ and IQR caclulation and retrieve normalization to use for later as norm.flag and extension to use
dnase_out=`Rscript ${SCRIPT_DIR}/dnase_background_cutting.R ${kmer} ${profile_naked_plus} ${dnase_iqr_fraction} ${msq_scale_thresh} ${max_extension} ${COMMON_FUNCTIONS}`
#returns a string with 4 values like: dnase_msq=x dnase_iqr=x norm_skip=x norm_extension=x"

#split and save values in seperate values since output is in one string (probably pretty inelegant sorry)
dnase_msq=`echo ${dnase_out} | perl -ne '/dnase_msq=(\d+\.\d+)\s+/; print $1;'`	
dnase_iqr=`echo ${dnase_out} | perl -ne '/dnase_iqr=(\d+\.?\d*)\s+/; print $1;'`
#skip flag: 0 = normalize with default extension as usual; 1 = extend the norm sclae window with extension and normalize; 2 = moderate extension still doesnt make the backgorund properly normalizable > dont normalize and flag this way
norm_skip_flag=`echo ${dnase_out} | perl -ne '/skip_flag=(\d+\.?\d*)\s+/; print $1;'`
norm_extension=`echo ${dnase_out} | perl -ne '/extension=(\d+\.?\d*)/; print $1;'`

#########################################################################
# 4 get base composition around kmer, store profile and plot i flagged  #
#########################################################################
#execute msq and iqr calculation for base composition
base_out=`Rscript ${SCRIPT_DIR}/base_composition.R ${kmer} ${infile_base_comp} ${base_plot_flag} ${base_plot_dir} ${base_bandwidth} ${COMMON_FUNCTIONS}`
#output is string like: base_msq=x base_iqr=x
#split string and save  values in seperate variables 
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
scale_plus=`echo ${fsr_out} | perl -ne '/scale_plus=(\d+\.?\d*)/; print $1;'`
scale_minus=`echo ${fsr_out} | perl -ne '/scale_minus=(\d+\.?\d*)/; print $1;'`

#####################################################
# write new line for one kmer to the assembly table #
#####################################################
# fsr.plus fsr.minus scale.plus scale.minus hscor dnase.msq dnase.iqr base.msq base.iqr

echo -e "${kmer}\t${id}\t${fsr_merged}\t${fsr_plus}\t${fsr_minus}\t${scale_plus}\t${scale_minus}\t${hscor}\t${dnase_msq}\t${dnase_iqr}\t${base_msq}\t${base_iqr}\t${norm_skip_flag}\t${norm_extension}" >>${table_file}

echo "kmer processed... ${id}"

done


