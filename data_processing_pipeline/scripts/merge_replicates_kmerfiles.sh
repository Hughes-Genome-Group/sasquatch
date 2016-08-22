##############################################################################################
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##############################################################################################
# Function: Data preprocessing; Merge kmer profiles from multiple replicates
# Usage: adjust parameters and copy to prompt
# Author: Ron Schwessinger
# Date: 27/04/2016

#!usr/bin/bash

# Select data path where tissue directories are stores [defailt = ./Sasquatch/data"
datadir="./Sasquatch/data"

# Define Tissue directory tag, single replicate dirs and name for the merged target
tag="Tissue_dir"

# declare as many reps as needed
rep1="${tag}_rep1"
rep2="${tag}_rep2"
rep3="${tag}_rep3"

merged="${tag}_merged"

# Select pnorm tag used for DNase propensity normalisation
pnorm_tag="h_ery_1"

# ==================================

### ===================== ###
### For 2 replicates run: ###
### ===================== ###

mkdir -p ${datadir}/${merged}/counts

cd ${datadir}/${merged}/counts

for i in 5 6 7
do

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_pnorm_${pnorm_tag}_plus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_pnorm_${pnorm_tag}_plus.txt >kmers_${i}_count_${merged}_pnorm_${pnorm_tag}_plus.txt

echo "... $i plus done"

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_pnorm_${pnorm_tag}_minus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_pnorm_${pnorm_tag}_minus.txt >kmers_${i}_count_${merged}_pnorm_${pnorm_tag}_minus.txt

echo "... $i minus done"

done

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch/data_processing_pipeline/scripts/merge_read_stats.pl ${datadir}/${rep1}/read_stats.txt ${datadir}/${rep2}/read_stats.txt | tee ${datadir}/${merged}/read_stats.txt

cd $datadir



### ===================== ###
### For 3 replicates run: ###
### ===================== ###

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


# For more rplicates adjust the commands 
