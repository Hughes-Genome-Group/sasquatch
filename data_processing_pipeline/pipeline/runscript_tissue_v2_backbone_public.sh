#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M ron.schwessinger@ndcls.ox.ac.uk
#$ -m eas
#$ -j y
#$ -o /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -N bb_ag04449_rep2

#qsub /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/pipeline/runscript_tissue_v2_backbone_public.sh

#set -x

#path to pipeline directory and additional scripts
PIPE_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/pipeline/
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/


##################
# Set PARAMETERS #
##################

#specify if "mouse" or "human"
ORGANISM="human"

#genome Build ["hg18", "hg19", "mm9"] currently choosable
BUILD='hg19'

#IDtag to produce output directory and name the files
IDTAG="ENCODE_UW_AG04449_rep2"

#specify if DNaseI or ATAC data ("DNase" or "ATAC")
DATA_TYPE="DNase"

#output directory
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/${ORGANISM}/${DATA_TYPE}/${IDTAG}

#type of sequencing ["singleend" / "pairedend"] 
SEQ_TYPE="singleend"


# ================================================================================ #
# PRESET SOME USEFULE FTP LINKS 						   #
# Source: ENCODE  "duke" "UW" "UW_DGF" "UW_mouse"; chose download paths 	   #
# ================================================================================ #

DWN_SOURCE="UW"

BAM_NAME="wgEncodeUwDnaseAg04449AlnRep2.bam"
PEAK_NAME="wgEncodeUwDnaseAg04449PkRep2"
PEAK_FILE=${PEAK_NAME}.narrowPeak

#select DOWNLOAD PATH accroding to DWN_SOURCE
case "${DWN_SOURCE}" in
	
	duke)
		#DUKe golden path
		DOWNLOAD_PATH='hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/'
	;;
	UW)
		#UW golden path
		DOWNLOAD_PATH='hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/'
	;;
	UW_DGF)
		#UW DGF golden path
		DOWNLOAD_PATH='hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDgf/'

	;;
	UW_mouse)
		#UW DGF golden path
		DOWNLOAD_PATH='hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeUwDnase/'

	;;
esac

#set future bamfile path
BAM_FILE=${OUTPUT_DIR}/${BAM_NAME}
#the perl script handling the region file is a bit lazy coded. it basically decides based on the end of of the file name (.gff or (.bed or .narroweak)) in which columns is has to look for the chromosome and start and stop coordinates
REGIONS_FILE=${OUTPUT_DIR}/${PEAK_FILE}
REGIONS_FILE_PLOIDY_FILTERED=${OUTPUT_DIR}/${PEAK_NAME}_ploidy_filtered.bed

# ======================================================================================= #
# Configure Paths to genomes and chr sizes files for required organisms and genome builds #
# only required those specified that are atully used                                      #
# ======================================================================================= #

#hg19
REF_GENOME_hg19=/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
BIGWIG_CHRSIZES_hg19=/t1-data/user/config/bigwig/hg19_sizes.txt
PLOIDY_REGIONS_hg19=/t1-data1/WTSA_Dev/rschwess/database_assembly/region_exclude/hg19/wgEncodeDukeMapabilityRegionsExcludable.bed

#hg18
REF_GENOME_hg18=/databank/raw/hg18_full/hg18_full.fa
BIGWIG_CHRSIZES_hg18=/t1-data/user/config/bigwig/hg18_sizes.txt
PLOIDY_REGIONS_hg18=/t1-data1/WTSA_Dev/rschwess/database_assembly/region_exclude/hg18/wgEncodeDukeRegionsExcluded.bed

#mm9
REF_GENOME_mm9=/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa
BIGWIG_CHRSIZES_mm9=/t1-data/user/config/bigwig/mm9_sizes.txt
PLOIDY_REGIONS_mm9=/t1-data1/WTSA_Dev/rschwess/database_assembly/region_exclude/mm9/Ploidy_mm9_sorted.bed

# ===================================================================================================== #
# Configure Paths to genome wide background DNase/ATAC cut propensities according to organism and build #
# only required those specified that are atully used                                      		#
# ===================================================================================================== #

#human DNase
PROPENSITY_PLUS_human_dnase=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_h_ery_1/pnorm/hg18_h_ery_1_propensities_plus_merged
PROPENSITY_MINUS_human_dnase=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_h_ery_1/pnorm/hg18_h_ery_1_propensities_minus_merged

#human ATAC
PROPENSITY_PLUS_human_atac==/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_h_ery_1_atac/pnorm/cut_kmer_6_hg18_h_ery_1_atac_plus_merged_propensities
PROPENSITY_MINUS_human_atac=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_h_ery_1_atac/pnorm/cut_kmer_6_hg18_h_ery_1_atac_minus_merged_propensities

#mouse DNase
PROPENSITY_PLUS_mouse_dnase=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_m_ery_1/pnorm/cut_kmer_6_mm9_m_ery_1_plus_merged_propensities
PROPENSITY_MINUS_mouse_dnase=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_m_ery_1/pnorm/cut_kmer_6_mm9_m_ery_1_minus_merged_propensities

#mouse ATAC
PROPENSITY_PLUS_mouse_atac=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_m_ery_1_atac/pnorm/cut_kmer_6_atac_mm9_m_ery_1_atac_plus_merged_propensities
PROPENSITY_MINUS_mouse_atac=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_m_ery_1_atac/pnorm/cut_kmer_6_atac_m_ery_1_atac_minus_merged_propensities


# ============================================================================================================== #
# Select reference genome, mapability (ploidy) regions to remove from the footprint signal and                   #
# background cut propensities according to organism and genome build 						 #
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

		#source of background cut propensities (ATAC or DNase1)
		case "${DATA_TYPE}" in
	
			DNase)
				pnormsource="h_ery_1"
				PROPENSITY_PLUS=${PROPENSITY_PLUS_human_dnase}
				PROPENSITY_MINUS=${PROPENSITY_MINUS_human_dnase}			
			;;

			ATAC)
				pnormsource="atac"
				PROPENSITY_PLUS=${PROPENSITY_PLUS_human_atac}
				PROPENSITY_MINUS=${PROPENSITY_MINUS_human_atac}
			;;
		esac

	;;

	mouse)

		case "${BUILD}" in

			mm9)	
				REF_GENOME=${REF_GENOME_mm9}
				BIGWIG_CHRSIZES=${BIGWIG_CHRSIZES_mm9}
				PLOIDY_REGIONS=${PLOIDY_REGIONS_mm9}						
			;;
		esac	

		#source of cut propensities (ATAC or DNase1)
		case "${DATA_TYPE}" in
	
			DNase)
				pnormsource="mm9"
				PROPENSITY_PLUS=${PROPENSITY_PLUS_mouse_dnase}
				PROPENSITY_MINUS=${PROPENSITY_MINUS_mouse_dnase}			
			;;

			ATAC)
				pnormsource="atac_mm9"
				PROPENSITY_PLUS=${PROPENSITY_PLUS_mouse_atac}
				PROPENSITY_MINUS=${PROPENSITY_MINUS_mouse_atac}
			;;
		esac

	;;
esac


####################
# START PROCESSING #
####################

echo "starting at ..."
date

#### --------------------------------------
#### Download Files
#### --------------------------------------
mkdir -p ${OUTPUT_DIR}

cd ${OUTPUT_DIR}

rsync -a -P rsync://${DOWNLOAD_PATH}/${BAM_NAME} ${OUTPUT_DIR}

rsync -a -P rsync://${DOWNLOAD_PATH}/${BAM_NAME}.bai ${OUTPUT_DIR}

rsync -a -P rsync://${DOWNLOAD_PATH}/${PEAK_FILE}.gz ${OUTPUT_DIR}

#wget -nv "http://${DOWNLOAD_PATH}/${BAM_NAME}"

#wget -nv "http://${DOWNLOAD_PATH}/${BAM_NAME}.bai"

#wget -nv "http://${DOWNLOAD_PATH}/${PEAK_FILE}.gz"

gzip -d ${OUTPUT_DIR}/${PEAK_FILE}.gz

### --------------------------------------
### filter regions file for ploidy regions
### --------------------------------------

#rather specific to our server system /not required when bedtools is accessible
module add bedtools

bedtools intersect -v -a ${REGIONS_FILE} -b ${PLOIDY_REGIONS} >${REGIONS_FILE_PLOIDY_FILTERED}

##### -----------------
##### Submit Footprints
##### -----------------

fpid=`qsub -N fp_${IDTAG} -v OUTPUT_DIR=${OUTPUT_DIR},DATA_TYPE=${DATA_TYPE},BAM_FILE=${BAM_FILE},SCRIPT_DIR=${SCRIPT_DIR},BIGWIG_CHRSIZES=${BIGWIG_CHRSIZES},SEQ_TYPE=${SEQ_TYPE},IDTAG=${IDTAG} ${PIPE_DIR}/runscript_tissue_v2_footprinting.sh | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

echo "Footprinting Job $fpid submitted"

### ------------------------
### submit count k-mers pnorm
### ------------------------
# make counts directory
COUNTS=${OUTPUT_DIR}/counts
mkdir -p ${COUNTS}

countpnormid=`qsub -N kmerpn_${IDTAG} -hold_jid $fpid -v COUNTS=${COUNTS},OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},SCRIPT_DIR=${SCRIPT_DIR},REF_GENOME=${REF_GENOME},IDTAG=${IDTAG},PROPENSITY_PLUS=${PROPENSITY_PLUS},PROPENSITY_MINUS=${PROPENSITY_MINUS},pnormsource=${pnormsource} ${PIPE_DIR}/runscript_tissue_v2_kmercount_pnorm.sh  | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

##without hold
# countpnormid=`qsub -N kmerpn_${IDTAG} -v COUNTS=${COUNTS},OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},SCRIPT_DIR=${SCRIPT_DIR},REF_GENOME=${REF_GENOME},IDTAG=${IDTAG},PROPENSITY_PLUS=${PROPENSITY_PLUS},PROPENSITY_MINUS=${PROPENSITY_MINUS},pnormsource=${pnormsource} ${PIPE_DIR}/runscript_tissue_v2_kmercount_pnorm.sh  | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

echo "K-mer counting P-norm Job $countpnormid submitted"

### -----------------------------------------------------
### count reads, peaks and reads in peaks into stats file
### -----------------------------------------------------

echo "Total reads: `samtools view ${BAM_FILE} | wc -l`" >${OUTPUT_DIR}/read_stats.txt
echo "Number of Peaks:: `wc -l ${REGIONS_FILE_PLOIDY_FILTERED}`" >>${OUTPUT_DIR}/read_stats.txt 
echo "Reads in Peaks:: `bedtools intersect -wa -a ${BAM_FILE} -b ${REGIONS_FILE_PLOIDY_FILTERED} | wc -l`" >>${OUTPUT_DIR}/read_stats.txt 

### --------------------------------------------------------------------------------------------
### Submit cleaner .sh for removing temporary stored bam and footprint files after jobs finished
### --------------------------------------------------------------------------------------------
qsub -N cleaner_${IDTAG} -hold_jid $fpid,$countpnormid -v OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE=${REGIONS_FILE},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},BAM_FILE=${BAM_FILE} ${PIPE_DIR}/runscript_tissue_v2_cleaner.sh

