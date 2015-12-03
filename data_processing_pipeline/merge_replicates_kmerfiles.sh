#!usr/bin/bash

datadir="/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human/DNase"

rep1="ENCODE_GM12878_UW_peak_refined_rep1"
rep2="ENCODE_GM12878_UW_peak_refined_rep2"

merged="ENCODE_GM12878_UW_peak_refined_merged"

mkdir -p ${datadir}/${merged}/counts

cd ${datadir}/${merged}/counts

for i in 5 6 7
do

#perl /t1-data1/WTSA_Dev/rschwess/scripts/dnase_tissue_motif/workflow/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_plus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_plus.txt >kmers_${i}_count_${merged}_plus.txt

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch_v19082015/data_processing_pipeline/scripts/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_pnorm_JH60_plus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_pnorm_JH60_plus.txt >kmers_${i}_count_${merged}_pnorm_JH60_plus.txt

echo "... $i plus done"

#perl /t1-data1/WTSA_Dev/rschwess/scripts/dnase_tissue_motif/workflow/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_minus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_minus.txt >kmers_${i}_count_${merged}_minus.txt

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch_v19082015/data_processing_pipeline/scripts/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_pnorm_JH60_minus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_pnorm_JH60_minus.txt >kmers_${i}_count_${merged}_pnorm_JH60_minus.txt

echo "... $i minus done"

done

cd $datadir


### ====== ###
### 3 reps ###
### ====== ###

datadir="/t1-data1/WTSA_Dev/rschwess/database_assembly/idx_correct_assembly/human/dnase"

rep1="ENCODE_HepG2_duke_rep1"
rep2="ENCODE_HepG2_duke_rep2"
rep3="ENCODE_HepG2_duke_rep3"

merged="ENCODE_HepG2_duke_merged"

mkdir -p ${datadir}/${merged}/counts

cd ${datadir}/${merged}/counts

for i in 5 6 7
do

#perl /t1-data1/WTSA_Dev/rschwess/scripts/dnase_tissue_motif/workflow/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_plus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_plus.txt ${datadir}/${rep3}/counts/kmers_${i}_count_${rep3}_plus.txt >kmers_${i}_count_${merged}_plus.txt

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch_v19082015/data_processing_pipeline/scripts/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_pnorm_JH60_plus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_pnorm_JH60_plus.txt  ${datadir}/${rep3}/counts/kmers_${i}_count_${rep3}_pnorm_JH60_plus.txt >kmers_${i}_count_${merged}_pnorm_JH60_plus.txt

echo "... $i plus done"

#perl /t1-data1/WTSA_Dev/rschwess/scripts/dnase_tissue_motif/workflow/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_minus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_minus.txt ${datadir}/${rep3}/counts/kmers_${i}_count_${rep3}_minus.txt >kmers_${i}_count_${merged}_minus.txt

perl /t1-data1/WTSA_Dev/rschwess/Sasquatch_offline/Sasquatch_v19082015/data_processing_pipeline/scripts/merge_replicates_kmerfiles.pl ${datadir}/${rep1}/counts/kmers_${i}_count_${rep1}_pnorm_JH60_minus.txt ${datadir}/${rep2}/counts/kmers_${i}_count_${rep2}_pnorm_JH60_minus.txt ${datadir}/${rep3}/counts/kmers_${i}_count_${rep3}_pnorm_JH60_minus.txt >kmers_${i}_count_${merged}_pnorm_JH60_minus.txt

echo "... $i minus done"

done

cd  $datadir



