#!usr/bin/bash

#$ -cwd
#$ -q batchq
#$ -M rschwess
#$ -m eas
#$ -e /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -o /t1-data1/WTSA_Dev/rschwess/clustereo
#$ -N bb_k562_upd_pipe

#qsub /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch_v19082015/data_processing_pipeline/pipeline/tissue_v2/runscript_tissue_v2_backbone_public_updated_2reps.sh

#set -x

#path to pipeline directory and additional scripts
PIPE_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/pipeline/tissue_v2
SCRIPT_DIR=/t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts

#specify if "mouse" or "human"
ORGANISM="human"

#genome Build ["hg18", "hg19", "mm9"] currently choosable
BUILD='hg19'

#IDtag to produce output directory and name the files
IDTAG="ENCODE_K562_UW"

#specify if DNaseI or ATAC data ("DNase" or "ATAC")
DATA_TYPE="DNase"

#output directory
OUTPUT_DIR=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/${ORGANISM}/${DATA_TYPE}/update_pipe/${IDTAG}

#type of sequencing ["singleend" / "pairedend"] 
SEQ_TYPE="singleend"

#Source "Duke" "UW" "UW_DGF"; download path is chosen accordingly
DWN_SOURCE="UW"

BAM_NAME1="wgEncodeUwDnaseK562AlnRep1.bam"
BAM_NAME2="wgEncodeUwDnaseK562AlnRep2.bam"
#set names to assign
BAM_MERGED="wgEncodeUwDnaseK562Alnmerged.bam"

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
BAM_FILE1=${OUTPUT_DIR}/${BAM_NAME1}
BAM_FILE2=${OUTPUT_DIR}/${BAM_NAME2}
BAM_FILE=${OUTPUT_DIR}/${BAM_MERGED}

REGIONS_FILE=${OUTPUT_DIR}/merged_filtered_peaks_${IDTAG}.bed
REGIONS_FILE_PLOIDY_FILTERED=${OUTPUT_DIR}/merged_filtered_peaks_${IDTAG}_ploidy_filtered.bed

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
				PLOIDY_REGIONS='/t1-data1/WTSA_Dev/rschwess/database_assembly/region_exclude/hg19/wgEncodeDukeMapabilityRegionsExcludable.bed'
			;;

			hg18)
				REF_GENOME="/databank/raw/hg18_full/hg18_full.fa"
				BIGWIG_CHRSIZES='/hts/data0/config/bigwig/hg18_sizes.txt'
				#ploidy regions to filter
				PLOIDY_REGIONS='/t1-data1/WTSA_Dev/rschwess/database_assembly/region_exclude/hg18/wgEncodeDukeMapabilityRegionsExcludable.bed'
			;;

		esac

		#source of cut propensities (ATAC or DNase1)
		case "${DATA_TYPE}" in
	
			DNase)
				pnormsource="JH60"
				PROPENSITY_PLUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH60/pnorm/hg18_human_JH60_propensities_plus_merged
				PROPENSITY_MINUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH60/pnorm/hg18_human_JH60_propensities_minus_merged
			
			;;

			ATAC)
				pnormsource="atac"
				PROPENSITY_PLUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH_atac/pnorm/cut_kmer_6_hg18_plus_merged_propensities
				PROPENSITY_MINUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/hg18_human_JH_atac/pnorm/cut_kmer_6_hg18_minus_merged_propensities

			;;
		esac

	;;

	mouse)

		case "${BUILD}" in

			mm9)	
				REF_GENOME="/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
				BIGWIG_CHRSIZES='/hts/data0/config/bigwig/mm9_sizes.txt'
				#ploidy regions to filter
				PLOIDY_REGIONS='/t1-data1/WTSA_Dev/rschwess/database_assembly/region_exclude/mm9/Ploidy_mm9_sorted.bed'							
			;;
		esac	

		#source of cut propensities (ATAC or DNase1)
		case "${DATA_TYPE}" in
	
			DNase)
				pnormsource="mm9"
				PROPENSITY_PLUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid/pnorm/cut_kmer_6_mm9_plus_merged_propensities
				PROPENSITY_MINUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid/pnorm/cut_kmer_6_mm9_minus_merged_propensities
			
			;;

			ATAC)
				pnormsource="atac_mm9"
				PROPENSITY_PLUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid_atac/pnorm/cut_kmer_6_atac_mm9_plus_merged_propensities
				PROPENSITY_MINUS=/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/background/mm9_mouse_erythroid_atac/pnorm/cut_kmer_6_atac_mm9_minus_merged_propensities

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
### Add Modules
### --------------------------------------

module add macs
module add bedtools
module add samtools/0.1.19

### --------------------------------------
### Download Files
### --------------------------------------
mkdir -p ${OUTPUT_DIR}

cd ${OUTPUT_DIR}

# #download bam files
rsync -a -P rsync://${DOWNLOAD_PATH}/${BAM_NAME1} ${OUTPUT_DIR}
rsync -a -P rsync://${DOWNLOAD_PATH}/${BAM_NAME2} ${OUTPUT_DIR}
rsync -a -P rsync://${DOWNLOAD_PATH}/${BAM_NAME1}.bai ${OUTPUT_DIR}
rsync -a -P rsync://${DOWNLOAD_PATH}/${BAM_NAME2}.bai ${OUTPUT_DIR}

echo "Downloading finished ..."

### --------------------------------------
### Call Peaks (MACS default) merge and filter for present in both
### --------------------------------------

cd ${OUTPUT_DIR}

#rep1 call
macs -t ${BAM_FILE1} --name Rep1 --format="BAM"
TEMP_REP_PEAK1=${OUTPUT_DIR}/Rep1_peaks.bed
if [ $? -ne 0 ]
then
  echo "Error in MACS peak calling.. not enough reads?"
  exit -1
fi
#filter for chrM
grep -v chrM ${TEMP_REP_PEAK1} >temp
mv -f temp ${TEMP_REP_PEAK1}

#rep2 call
macs -t ${BAM_FILE2} --name Rep2 --format="BAM"
TEMP_REP_PEAK2=${OUTPUT_DIR}/Rep2_peaks.bed
if [ $? -ne 0 ]
then
  echo "Error in MACS peak calling.. not enough reads?"
  exit -1
fi
#filter for chrM
grep -v chrM ${TEMP_REP_PEAK2} >temp
mv -f temp ${TEMP_REP_PEAK2}

#With both: merge and sort
cat ${TEMP_REP_PEAK1} ${TEMP_REP_PEAK2} >merged
sort -k1,1 -k2,2n merged >merged.sorted
bedtools merge -d 10 -i merged.sorted >merged.sorted.collapsed

#filter for presence i both peak files
bedtools intersect -u -a merged.sorted.collapsed -b ${TEMP_REP_PEAK1} >temp
bedtools intersect -u -a temp -b ${TEMP_REP_PEAK2} >${REGIONS_FILE}
if [ $? -ne 0 ]
then
  echo "Error in merging and intersecting"
  exit -1
fi

#remove temp files
rm -f merged merged.sorted.collapsed temp merged.sorted

echo "Peak Calling finished ..."

### --------------------------------------
### Merge Bam files
### --------------------------------------

samtools merge ${BAM_FILE} ${BAM_FILE1} ${BAM_FILE2}

echo "Bam Files mergeing finished ..."

### --------------------------------------
### filter regions file for ploidy regions
### --------------------------------------

bedtools intersect -v -a ${REGIONS_FILE} -b ${PLOIDY_REGIONS} >${REGIONS_FILE_PLOIDY_FILTERED}

echo "Bed File is filtered ..."

### -----------------
### Submit Footprints
### -----------------

fpid=`qsub -N fp_${IDTAG} -v OUTPUT_DIR=${OUTPUT_DIR},DATA_TYPE=${DATA_TYPE},BAM_FILE=${BAM_FILE},SCRIPT_DIR=${SCRIPT_DIR},BIGWIG_CHRSIZES=${BIGWIG_CHRSIZES},SEQ_TYPE=${SEQ_TYPE},IDTAG=${IDTAG} ${PIPE_DIR}/runscript_tissue_v2_footprinting.sh | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

echo "Footprinting Job $fpid submitted"

### ------------------------
### submit count k-mers pnorm
### ------------------------
make counts directory
COUNTS=${OUTPUT_DIR}/counts
mkdir -p ${COUNTS}

countpnormid=`qsub -N kmerpn_${IDTAG} -hold_jid $fpid -v COUNTS=${COUNTS},OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},SCRIPT_DIR=${SCRIPT_DIR},REF_GENOME=${REF_GENOME},IDTAG=${IDTAG},PROPENSITY_PLUS=${PROPENSITY_PLUS},PROPENSITY_MINUS=${PROPENSITY_MINUS},pnormsource=${pnormsource} ${PIPE_DIR}/runscript_tissue_v2_kmercount_pnorm.sh  | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

# countpnormid=`qsub -N kmerpn_${IDTAG} -v COUNTS=${COUNTS},OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE_PLOIDY_FILTERED=${REGIONS_FILE_PLOIDY_FILTERED},SCRIPT_DIR=${SCRIPT_DIR},REF_GENOME=${REF_GENOME},IDTAG=${IDTAG},PROPENSITY_PLUS=${PROPENSITY_PLUS},PROPENSITY_MINUS=${PROPENSITY_MINUS},pnormsource=${pnormsource} ${PIPE_DIR}/runscript_tissue_v2_kmercount_pnorm.sh  | perl -ne '$_=~/\s+(\d+)\s+/; print $1;'`

echo "K-mer counting P-norm Job $countpnormid submitted"

### --------------------------------------
### count reads, peaks and reads in peaks
### --------------------------------------

echo "Total reads: `samtools view ${BAM_FILE} | wc -l`" >${OUTPUT_DIR}/read_stats.txt
echo "Number of Peaks:: `wc -l ${REGIONS_FILE_PLOIDY_FILTERED}`" >>${OUTPUT_DIR}/read_stats.txt 
echo "Reads in Peaks:: `bedtools intersect -abam -wa -a ${BAM_FILE} -b ${REGIONS_FILE_PLOIDY_FILTERED} | wc -l`" >>${OUTPUT_DIR}/read_stats.txt 

### -------------------------------------------------------------------------------
### Submit cleaner .sh for removing temporary stored bam files after jobs finished
### -------------------------------------------------------------------------------

qsub -N cleaner_${IDTAG} -hold_jid $fpid,$countpnormid -v OUTPUT_DIR=${OUTPUT_DIR},REGIONS_FILE=${REGIONS_FILE},BAM_FILE=${BAM_FILE},BAM_FILE1=${BAM_FILE1},BAM_FILE2=${BAM_FILE2} ${PIPE_DIR}/runscript_tissue_v2_cleaner.sh

