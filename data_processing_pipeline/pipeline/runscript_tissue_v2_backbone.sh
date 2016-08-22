#!/bin/bash

##############################################################################################
##                                                                                          ##
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##                                                                                          ##
## Data preprocessing: etting paramters and job submission backbone                         ##
##                                                                                          ##
##############################################################################################

# Options: for Sun Grid Engine Cluster Batch Queue System 
#$ -cwd
#$ -q batchq
#$ -m eas
#$ -j y
#$ -N sasq_data_preprocess_backbone

#################
# Requirements: # 
#################
#
# Scripts were written for use on a SUN Grid ENgine batch queue archtitecture 
# To run on other ssystems please adjust the submissions via "qsub" 
#
# Requires samtools and optionally bedtools for blacklisted region filtering
#
# Notes:
# The kmercount step takes the most times and ressources depending on yoursystem and input files 
# (bam file size  and number of regions of interest)
# We haven't tested the scripts on wide range of systems yet so the following estiamtes 
# are based on our systems only
# 
# Rough time estimates (depending on input and system):
# footprinting ~ 30 min - 1h
# kmer count ~ 8 - 24 h

###########
# Run as: #
###########

# unix like 
# sh /path_to_sasquatch/Sasquatch/data_processing_pipeline/pipeline/runscript_tissue_v2_backbone.sh
#
# SGE queue system
# qsub /path_to_sasquatch/Sasquatch/data_processing_pipeline/pipeline/runscript_tissue_v2_backbone.sh

### ======================================= ###
### START OF PATHS AND PARAMETERS TO ADJUST ###
### ======================================= ###

# set verbose option for testing and debugging
# set -x

# Set PIPELINE AND OUTPUT Paths ============================================= #
# Full path to your Sasquatch copy
SASQ_PATH=/path_to_sasquatch/Sasquatch

# path to pipeline runscripts [default = ${SASQ_PATH}/data_processing_pipeline/pipeline]
PIPE_DIR=${SASQ_PATH}/data_processing_pipeline/pipeline

# path to used pipeline scripts [default = ${SASQ_PATH}/data_processing_pipeline/scripts]
SCRIPT_DIR=${SASQ_PATH}/data_processing_pipeline/scripts

# set output directory [default = ${SASQ_PATH}/data/${ORGANISM}/${DATA_TYPE}/${IDTAG}]
OUTPUT_DIR=${SASQ_PATH}/${ORGANISM}/${DATA_TYPE}/${IDTAG}


# Set INPUT ================================================================= #
# path to aligned reads bam file
BAM_FILE=/t1-data1/WTSA_Dev/rschwess/data_storage/DNase_for_CGAT/Keratinocyte-R1-DNase.bam

# full path to peak file [.bed .gff or .narrowPeak supported] e.g. for .bed only 3 columns required
REGIONS_FILE=/t1-data1/WTSA_Dev/rschwess/data_storage/DNase_for_CGAT/Keratinocyte-R1-DNase.bed
# the perl script handling the region file decides based on the end of the file name (.gff or (.bed or .narroPeak)) in which columns it has to look for the chromosome and start and stop coordinates


# Set PARAMETERS ========================================================== #
# specify organism ["mouse" or "human"]
ORGANISM="human"	

# genome Build ["hg18", "hg19", or "mm9"]
BUILD='hg19'

# IDtag to produce output directory and name the files
IDTAG="my_new_tissue_rep1"

# specify if DNaseI or ATAC data [default = "DNase"] 
# ("ATAC" only supported for test purposes)
# you can select ATAC to adjust the fooptint signal for the reportet transposase shifts
# However we do not publish specific ATAC background and do support the use of ATAC data only for testing purposes
DATA_TYPE="DNase"

# select type of sequencing ["singleend" / "pairedend"]
SEQ_TYPE="singleend"

# select if to filter your bedfile with provided ploidic regions
# We recommend to filter your peaks with a blacklist of unmappable or highly repetitive ("ploidic") regions suited to you genome build
# example files:
# hg19: 
# hg18: 
# mm9: 
# Requires "bedtools" installed
FILTER_PLOIDY_REGIONS=true

# identifier name of the peak file to produce a ploidy filtered peaks [default = `basename ${PEAK_FILE} .bed`]
PEAK_NAME=`basename ${PEAK_FILE} .bed`
# how to call the ploidy filtered bed file
REGIONS_FILE_PLOIDY_FILTERED=${OUTPUT_DIR}/${PEAK_NAME}_ploidy_filtered.bed

# ============================================================================ #
# Configure Paths to: 
#	genomes reference sequence .fa files 
#	chr sizes files; two tab separated columns (chrN, size)
#	ploidic regions to exclude (bed format)
# only requires those specified that are atully used                                      
# ============================================================================ #

# hg19
REF_GENOME_hg19=/full_path_to_mm9genome_fasta/genome.fa
BIGWIG_CHRSIZES_hg19=/full_path_to_mm9_chrom_sizes/mm9_sizes.txt
PLOIDY_REGIONS_hg19=/full_path_to_mm9_ploidic_regions_file/blacklist_ploidy_mm9.bed

# hg18
REF_GENOME_hg18=/full_path_to_mm9genome_fasta/genome.fa
BIGWIG_CHRSIZES_hg18=/full_path_to_mm9_chrom_sizes/mm9_sizes.txt
PLOIDY_REGIONS_hg18=/full_path_to_mm9_ploidic_regions_file/blacklist_ploidy_mm9.bed

# mm9
REF_GENOME_mm9=/full_path_to_mm9genome_fasta/genome.fa
BIGWIG_CHRSIZES_mm9=/full_path_to_mm9_chrom_sizes/mm9_sizes.txt
PLOIDY_REGIONS_mm9=/full_path_to_mm9_ploidic_regions_file/blacklist_ploidy_mm9.bed

# =========================================================================== #
# Configure Paths to genome wide background DNase cut propensities 
# Download preprocessed background files from our webtool repository:
# http://apps.molbiol.ox.ac.uk/sasquatch/cgi-bin/foot.cgi
# and extract them for example into 
# /full_path_to_my_sasquatch_copy/Sasquatch/data/human/
# only requires those specified that are actully used
# Background files are only organism dependent
# =========================================================================== #

# human DNase
PROPENSITY_PLUS_human_dnase=/${SASQ_PATH}/data/human/Background_human_ery_hg18_1_propensities_plus_merged
PROPENSITY_MINUS_human_dnase=/${SASQ_PATH}/data/human/Background_human_ery_hg18_1_propensities_minus_merged

# mouse DNase
PROPENSITY_PLUS_mouse_dnase=/${SASQ_PATH}/data/mouse/Background_mouse_ery_mm9_1_propensities_plus_merged
PROPENSITY_MINUS_mouse_dnase=/${SASQ_PATH}/data/mouse/Background_mouse_ery_mm9_1_propensities_minus_merged

### =========================== ###
### END OF PARAMETERS TO ADJUST	###
### =========================== ###

# ============================================================================================================== #
# Auto select reference genome, mapability (ploidy) regions to remove from the footprint signal and              
# background cut propensities according to organism and genome build 						 
# ============================================================================================================== #
case "${ORGANISM}" in
	
	human)
		#select ref genome ploidy regions and chromosome sizes
		case "${BUILD}" in

			hg19)
				REF_GENOME=${REF_GENOME_hg19}
				BIGWIG_CHRSIZES=${BIGWIG_CHRSIZES_hg19}
				PLOIDY_REGIONS=${PLOIDY_REGIONS_hg19}
			;;

			hg18)
				REF_GENOME=${REF_GENOME_hg18}
				BIGWIG_CHRSIZES=${BIGWIG_CHRSIZES_hg18}
				PLOIDY_REGIONS=${PLOIDY_REGIONS_hg18}
			;;

		esac

		#source of background cut propensities (DNase)
		pnormsource="h_ery_1"
		PROPENSITY_PLUS=${PROPENSITY_PLUS_human_dnase}
		PROPENSITY_MINUS=${PROPENSITY_MINUS_human_dnase}			
		
	;;

	mouse)

		case "${BUILD}" in

			mm9)	
				REF_GENOME=${REF_GENOME_mm9}
				BIGWIG_CHRSIZES=${BIGWIG_CHRSIZES_mm9}
				PLOIDY_REGIONS=${PLOIDY_REGIONS_mm9}						
			;;
		esac	

		#source of cut propensities (DNase1)
		pnormsource="m_ery_1"
		PROPENSITY_PLUS=${PROPENSITY_PLUS_mouse_dnase}
		PROPENSITY_MINUS=${PROPENSITY_MINUS_mouse_dnase}			
	;;
esac



####################
# START PROCESSING #
####################

echo "starting at ..."
date

# ------------------------------------------------------------------------------
# 0) OPTIONAL filter regions file for ploidy regions (requires bedtools)
# ------------------------------------------------------------------------------

#rather specific to our server system /not required when bedtools is accessible
if [ "${FILTER_PLOIDIC_REGIONS}" = true ];
then
	bedtools intersect -v -a ${REGIONS_FILE} -b ${PLOIDY_REGIONS} >${REGIONS_FILE_PLOIDY_FILTERED}
else
	${REGIONS_FILE_PLOIDY_FILTERED} = ${PLOIDY_REGIONS}
fi

# -----------------------------------------------------------------------------
# 1) Submit Footprinting
# -----------------------------------------------------------------------------
fpid=`qsub -N fp_${IDTAG} -v OUTPUT_DIR=${OUTPUT_DIR},DATA_TYPE=${DATA_TYPE},BAM_FILE=${BAM_FILE},SCRIPT_DIR=${SCRIPT_DIR},BIGWIG_CHRSIZES=${BIGWIG_CHRSIZES},SEQ_TYPE=${SEQ_TYPE},IDTAG=${IDTAG} ${PIPE_DIR}/runscript_tissue_v2_footprinting.sh | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

echo "Footprinting Job $fpid submitted"

# ---------------------------------------------------------------------------
# 2) submit count k-mers with propensity norm
# ---------------------------------------------------------------------------
#make counts directory

#sleep to ensure footprinting has finsihed properly
sleep 15

COUNTS=${OUTPUT_DIR}/counts
mkdir -p ${COUNTS}
	
countpnormid=`qsub -N kmerpn_${IDTAG} -hold_jid $fpid -v SASQ_PATH=${SASQ_PATH},COUNTS=${COUNTS},OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},SCRIPT_DIR=${SCRIPT_DIR},REF_GENOME=${REF_GENOME},IDTAG=${IDTAG},PROPENSITY_PLUS=${PROPENSITY_PLUS},PROPENSITY_MINUS=${PROPENSITY_MINUS},pnormsource=${pnormsource} ${PIPE_DIR}/runscript_tissue_v2_kmercount_pnorm.sh  | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

##without hold
#countpnormid=`qsub -N kmerpn_${IDTAG} -v SASQ_PATH=${SASQ_PATH},COUNTS=${COUNTS},OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},SCRIPT_DIR=${SCRIPT_DIR},REF_GENOME=${REF_GENOME},IDTAG=${IDTAG},PROPENSITY_PLUS=${PROPENSITY_PLUS},PROPENSITY_MINUS=${PROPENSITY_MINUS},pnormsource=${pnormsource} ${PIPE_DIR}/runscript_tissue_v2_kmercount_pnorm.sh  | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

echo "K-mer counting P-norm Job $countpnormid submitted"

# -----------------------------------------------------------------------------
# 3) count reads, peaks and reads in peaks into stats file
# -----------------------------------------------------------------------------

echo "Total reads: `samtools view ${BAM_FILE} | wc -l`" >${OUTPUT_DIR}/read_stats.txt
echo "Number of Peaks:: `wc -l ${REGIONS_FILE_PLOIDY_FILTERED}`" >>${OUTPUT_DIR}/read_stats.txt 
echo "Reads in Peaks:: `bedtools intersect -wa -a ${BAM_FILE} -b ${REGIONS_FILE_PLOIDY_FILTERED} | wc -l`" >>${OUTPUT_DIR}/read_stats.txt 

### ---------------------------------------------------------------------------
### Submit cleaner .sh for removing temporary stored bam and footprint files after jobs finished
### Include if you like the input bam files and bed files and footprint wig files to be be removed once finished
### --------------------------------------------------------------------------------------------
### qsub -N cleaner_${IDTAG} -hold_jid $fpid,$countpnormid -v REGIONS_FILE=${REGIONS_FILE},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},BAM_FILE=${BAM_FILE} ${PIPE_DIR}/runscript_tissue_v2_cleaner.sh

