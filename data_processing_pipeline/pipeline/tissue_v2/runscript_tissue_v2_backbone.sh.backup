#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -e /hts/data4/rschwess/clustereo
#$ -o /hts/data4/rschwess/clustereo
#$ -N bb_human_ery_macs2_atac_oldnorm

#qsub /hts/data4/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/pipeline/tissue_v2/runscript_tissue_v2_backbone.sh
	
#verbose option for testing issues
#set -x

#path to pipeline directory and additional scripts
PIPE_DIR=/hts/data4/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/pipeline/tissue_v2
SCRIPT_DIR=/hts/data4/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/

##################
# Set PARAMETERS #
##################

#specify if "mouse" or "human"
ORGANISM="human"	

#genome Build ["hg18", "hg19", "mm9"] currently choosable
BUILD='hg18'

#IDtag to produce output directory and name the files
IDTAG="erythroid_hg18_macs2"

#specify if DNaseI or ATAC data ("DNase" or "ATAC")
DATA_TYPE="ATAC"

#type of sequencing ["singleend" / "pairedend"] 
SEQ_TYPE="pairedend"

#set output directory
OUTPUT_DIR=/hts/data4/rschwess/database_assembly/idx_correct_assembly/${ORGANISM}/${DATA_TYPE}/${IDTAG}/

#path to aligned reads bam file
BAM_FILE="/hts/data1/msuciu/ATAC-Seq/hProg/JH_d10/pipe/JH_d10_cd34pos_ATAC/hg18/filtered_Sorted.bam"
 
#identifier name of the peak file to produce the ploidy filtered peaks
PEAK_NAME="macs2_peakcall"
#full path to peak file
PEAK_FILE="${OUTPUT_DIR}/macs2_peakcall.bed"

#the perl script handling the region file is a bit lazy coded. it basically decides based on the end of of the file name (.gff or (.bed or .narroweak)) in which columns is has to look for the chromosome and start and stop coordinates
REGIONS_FILE=${PEAK_FILE}
REGIONS_FILE_PLOIDY_FILTERED=${OUTPUT_DIR}/${PEAK_NAME}_ploidy_filtered.bed

# ============================================================================================================== #
# Select reference genome, ploidy regions and cut propensities according to organism and genome build aligned to #
# ============================================================================================================== #
case "${ORGANISM}" in
	
	human)
		#select ref genome ploidy regions and chromosome sizes
		case "${BUILD}" in

			hg19)
				REF_GENOME="/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
				BIGWIG_CHRSIZES='/hts/data0/config/bigwig/hg19_sizes.txt'
				#ploidy regions to filter
				PLOIDY_REGIONS='/hts/data4/rschwess/database_assembly/region_exclude/hg19/wgEncodeDukeMapabilityRegionsExcludable.bed'
			;;

			hg18)
				REF_GENOME="/databank/raw/hg18_full/hg18_full.fa"
				BIGWIG_CHRSIZES='/hts/data0/config/bigwig/hg18_sizes.txt'
				#ploidy regions to filter
				PLOIDY_REGIONS='/hts/data4/rschwess/database_assembly/region_exclude/hg18/wgEncodeDukeRegionsExcluded.bed'
			;;

		esac

		#source of cut propensities (ATAC or DNase1)
		case "${DATA_TYPE}" in
	
			DNase)
				pnormsource="JH60"
				PROPENSITY_PLUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH60/pnorm/hg18_human_JH60_propensities_plus_merged
				PROPENSITY_MINUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH60/pnorm/hg18_human_JH60_propensities_minus_merged
			
			;;

			ATAC)
				pnormsource="atac"
				PROPENSITY_PLUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH_atac/pnorm/cut_kmer_6_hg18_plus_merged_propensities
				PROPENSITY_MINUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH_atac/pnorm/cut_kmer_6_hg18_minus_merged_propensities

			;;
		esac

	;;

	mouse)

		case "${BUILD}" in

			mm9)	
				REF_GENOME="/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
				BIGWIG_CHRSIZES='/hts/data0/config/bigwig/mm9_sizes.txt'
				#ploidy regions to filter
				PLOIDY_REGIONS='/hts/data4/rschwess/database_assembly/region_exclude/mm9/Ploidy_mm9_sorted.bed'							
			;;
		esac	

		#source of cut propensities (ATAC or DNase1)
		case "${DATA_TYPE}" in
	
			DNase)
				pnormsource="mm9"
				PROPENSITY_PLUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid/pnorm/cut_kmer_6_mm9_plus_merged_propensities
				PROPENSITY_MINUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid/pnorm/cut_kmer_6_mm9_minus_merged_propensities
			
			;;

			ATAC)
				pnormsource="atac_mm9"
				PROPENSITY_PLUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid_atac/pnorm/cut_kmer_6_atac_mm9_plus_merged_propensities
				PROPENSITY_MINUS=/hts/data4/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid_atac/pnorm/cut_kmer_6_atac_mm9_minus_merged_propensities

			;;
		esac

	;;
esac


####################
# START PROCESSING #
####################

echo "starting at ..."
date

### --------------------------------------
### filter regions file for ploidy regions
### --------------------------------------

module add bedtools

#bedtools intersect -v -a ${REGIONS_FILE} -b ${PLOIDY_REGIONS} >${REGIONS_FILE_PLOIDY_FILTERED}

### -----------------
### Submit Footprints
### -----------------

#fpid=`qsub -N fp_${IDTAG} -v OUTPUT_DIR=${OUTPUT_DIR},DATA_TYPE=${DATA_TYPE},BAM_FILE=${BAM_FILE},SCRIPT_DIR=${SCRIPT_DIR},BIGWIG_CHRSIZES=${BIGWIG_CHRSIZES},SEQ_TYPE=${SEQ_TYPE},IDTAG=${IDTAG} ${PIPE_DIR}/runscript_tissue_v2_footprinting.sh | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

#echo "Footprinting Job $fpid submitted"

### ------------------------
### submit count k-mers old
### ------------------------

##make counts directory
mkdir -p ${OUTPUT_DIR}/counts	
 -hold_jid $fpid
countoldid=`qsub -N kmerco_${IDTAG} -v OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},SCRIPT_DIR=${SCRIPT_DIR},REF_GENOME=${REF_GENOME},IDTAG=${IDTAG} ${PIPE_DIR}/runscript_tissue_v2_kmercountold.sh  | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

echo "K-mer counting Job $countoldid submitted"

### ------------------------
### submit count k-mers pnorm
### ------------------------
#make counts directory

#sleep 15

#COUNTS=${OUTPUT_DIR}/counts
#mkdir -p ${COUNTS}
#	
#countpnormid=`qsub -N kmerpn_${IDTAG} -hold_jid $fpid -v  COUNTS=${COUNTS},OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},SCRIPT_DIR=${SCRIPT_DIR},REF_GENOME=${REF_GENOME},IDTAG=${IDTAG},PROPENSITY_PLUS=${PROPENSITY_PLUS},PROPENSITY_MINUS=${PROPENSITY_MINUS},pnormsource=${pnormsource} ${PIPE_DIR}/runscript_tissue_v2_kmercount_pnorm.sh  | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

#countpnormid=`qsub -N kmerpn_${IDTAG} -v  COUNTS=${COUNTS},OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},SCRIPT_DIR=${SCRIPT_DIR},REF_GENOME=${REF_GENOME},IDTAG=${IDTAG},PROPENSITY_PLUS=${PROPENSITY_PLUS},PROPENSITY_MINUS=${PROPENSITY_MINUS},pnormsource=${pnormsource} ${PIPE_DIR}/runscript_tissue_v2_kmercount_pnorm.sh  | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

echo "K-mer counting P-norm Job $countpnormid submitted"

### --------------------------------------
### count reads, peaks and reads in peaks
### --------------------------------------

#echo "Total reads: `samtools view ${BAM_FILE} | wc -l`" >${OUTPUT_DIR}/read_stats.txt
#echo "Number of Peaks:: `wc -l ${REGIONS_FILE_PLOIDY_FILTERED}`" >>${OUTPUT_DIR}/read_stats.txt 
#echo "Reads in Peaks:: `bedtools intersect -abam -wa -a ${BAM_FILE} -b ${REGIONS_FILE_PLOIDY_FILTERED} | wc -l`" >>${OUTPUT_DIR}/read_stats.txt 

#### -------------------------------------------------------------------------------
#### Submit cleaner .sh for removing temporary stored bam files after jobs finished
#### -------------------------------------------------------------------------------
# qsub -N cleaner_${IDTAG} -hold_jid $fpid,$countpnormid -v REGIONS_FILE=${REGIONS_FILE},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},BAM_FILE=${BAM_FILE} ${PIPE_DIR}/runscript_tissue_v2_cleaner.sh

