### Explanation Example File for Comparing two kmers in a double plot fashion ###

# !!! I pruned the stats calculation as only the SFR rations are shown in sasquatch at the moment and the other stuff is rather unnecessary with the new normalization

#First set directory for scripts, output directory and location of the common_functions.R file that contains 
#all shared functions sourced from some of the latter scripts:
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/scripts			
COMMON_FUNCTIONS=${SCRIPT_DIR}/common_functions.R
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/idx_duke_testout/merge_strands

#Specify the directory where the background files (base frequency and naked dnaseI fibroblast cutting) are located,
#Currently this links to files for chr1 only but I'm close to finish the complete genome as background.
BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg19_ploidy_removed/counts

ORGANISM="human"

FRAG_TYPE="ATAC"

#select a background/normalization type (choose accordingly normalized kmer count files) there will pnorm dnase and pnorm atac files
#we have to check how much the different backgrounds differ to tell how much different norm types we need to have but at least one for atac and one for dnase will be needed (currently only dnase available)
# currently available "laza" (human fibroblast) "JH40" (human erythroid 40% mapped)
NORM_TYPE="JH60"

#I use the following TAG to create tissue specific sub directories and name the files accordingly when creating the tissue specific data,
#so that by setting a data directory and the correct tag you get access to the tissue specific count files later.
#for erythroid human example
TAG="human_erythroid_hg18"		
			
###INPUT###
#enter a k-mer of interest

kmer1="TGACGTC"	

kmer2="TCACGTC"

#get & store kmer length
kl1=`expr length $kmer1`
kl2=`expr length $kmer2`	

if [ $kl1 = $kl2 ] 
then
kl=$kl1
else
echo "Kmers dont have the same length"
fi


##############################################################
# DATA & BACKGROUNDs select according to organism and tissue #
##############################################################

#define DATA directory
DATA_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human${ORGANISM}/${FRAG_TYPE}/

### BACKGROUNDs ###
case "${ORGANISM}" in 

	human)

		case "${FRAG_TYPE}" in 
		
		DNase)	
			#Specify the directory where the background files (naked dnaseI fibroblast cutting) are located,
			BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH60/counts
			#defining full paths to naked background count files
			infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_JH60_plus_merged
			infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_hg18_human_JH60_minus_merged
			#Select accordingly normalized input files
			infile_plus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_${NORM_TYPE}_plus.txt
			infile_minus=${DATA_DIR}/${TISSUE}/counts/kmers_${kl}_count_${TISSUE}_pnorm_${NORM_TYPE}_minus.txt
		;;
		ATAC)
			BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH_atac/counts
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
			BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid/counts
			infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_mm9_plus_merged
			infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_mm9_minus_merged
			infile_plus=${DATA_DIR}/${TISSUE}/kmers_${kl}_count_${TISSUE}_pnorm_mm9_plus.txt
			infile_minus=${DATA_DIR}/${TISSUE}/kmers_${kl}_count_${TISSUE}_pnorm_mm9_minus.txt
		;;
		ATAC)
			BACKGROUND_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid_atac/counts
			infile_naked_plus=${BACKGROUND_DIR}/kmer_${kl}_mm9_atac_plus_merged
			infile_naked_minus=${BACKGROUND_DIR}/kmer_${kl}_mm9_atac_minus_merged
			infile_plus=${DATA_DIR}/${TISSUE}/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_plus.txt
			infile_minus=${DATA_DIR}/${TISSUE}/kmers_${kl}_count_${TISSUE}_pnorm_atac_mm9_minus.txt
		;;
		esac
	;;

esac

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

#######
#merge#
#######
profile_merged1=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus1} ${profile_minus1} ${FRAG_TYPE} ${kl}`
profile_merged2=`Rscript ${SCRIPT_DIR}/merge_profiles.R ${profile_plus2} ${profile_minus2} ${FRAG_TYPE} ${kl}`

################################################
# 2 single overlap plot of normalized profiles #
################################################

smooth_flag=1

#directory to write plot to
plot_dir=${OUTPUT_DIR}

Rscript ${SCRIPT_DIR}/overlap_plot_pnorm_normalized_merge.R ${kmer1} ${kmer2} ${profile_merged1} ${profile_merged2} ${profile_plus_count1} ${profile_plus_count2} ${smooth_flag} ${plot_dir} ${COMMON_FUNCTIONS}

