##############################################################################################
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##############################################################################################
# Function: Data preprocessing; Merge kmer profiles from multiple replicates
# Usage: adjust parameters and copy to prompt
# Author: Ron Schwessinger
# Date: 27/04/2016

#!usr/bin/bash

datadir="/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human/DNase"

rep1="ENCODE_UW_Monocytes_CD14plus_RO1746_rep1"
rep2="ENCODE_UW_Monocytes_CD14plus_RO1746_rep2"

merged="ENCODE_UW_Monocytes_CD14plus_RO1746_merged"

pnorm_tag="h_ery_1"

# ==================================

mkdir -p ${datadir}/${merged}/counts

cd ${datadir}/${merged}/counts

for i in 5 6 7
do

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_pnorm_${pnorm_tag}_plus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_pnorm_${pnorm_tag}_plus.txt >kmers_${i}_count_${merged}_pnorm_${pnorm_tag}_plus.txt

echo "... $i plus done"

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_pnorm_${pnorm_tag}_minus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_pnorm_${pnorm_tag}_minus.txt >kmers_${i}_count_${merged}_pnorm_{pnorm_tag}_minus.txt

echo "... $i minus done"

done

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/merge_read_stats.pl ${datadir}/${rep1}/read_stats.txt ${datadir}/${rep2}/read_stats.txt | tee ${datadir}/${merged}/read_stats.txt

cd $datadir


### ====== ###
### 3 reps ###
### ====== ###

datadir="/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human/DNase"

rep1="ENCODE_Duke_HeLa_S3_rep1"
rep2="ENCODE_Duke_HeLa_S3_rep2"
rep3="ENCODE_Duke_HeLa_S3_rep3"

merged="ENCODE_Duke_HeLa_S3_merged"

pnorm_tag="h_ery_1"

# ==================================

mkdir -p ${datadir}/${merged}/counts

cd ${datadir}/${merged}/counts

for i in 5 6 7
do

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_pnorm_${pnorm_tag}_plus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_pnorm_${pnorm_tag}_plus.txt  ${datadir}/${rep3}/counts/kmers_${i}_count_${rep3}_pnorm_${pnorm_tag}_plus.txt >kmers_${i}_count_${merged}_pnorm_${pnorm_tag}_plus.txt

echo "... $i plus done"

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_pnorm_${pnorm_tag}_minus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_pnorm_${pnorm_tag}_minus.txt ${datadir}/${rep3}/counts/kmers_${i}_count_${rep3}_pnorm_${pnorm_tag}_minus.txt >kmers_${i}_count_${merged}_pnorm_${pnorm_tag}_minus.txt

echo "... $i minus done"

done

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/merge_read_stats.pl ${datadir}/${rep1}/read_stats.txt ${datadir}/${rep2}/read_stats.txt ${datadir}/${rep3}/read_stats.txt | tee ${datadir}/${merged}/read_stats.txt

cd  $datadir



