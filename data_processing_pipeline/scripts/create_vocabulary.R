##############################################################################################
##                                                                                          ##
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##                                                                                          ##
## This program is free software: you can redistribute it and/or modify                     ##
## it under the terms of the GNU General Public License as published by                     ##
## the Free Software Foundation, either version 3 of the License, or                        ##
## (at your option) any later version.                                                      ##
##                                                                                          ##
## This program is distributed in the hope that it will be useful,                          ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of                           ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                             ##
## GNU General Public License for more details.	                                            ##
##                                                                                          ##
## You should have received a copy of the GNU General Public License                        ##
## along with this program. If not, see <http://www.gnu.org/licenses/>.                     ##
##                                                                                          ##
## Contact: Ron Schwessinger, ron.schwessinger@ndcls.ox.ac.uk                               ##
##          CBRG, genmail@molbiol.ox.ac.uk                                                  ##
##          Jim Hughes, jim.hughes@imm.ox.ac.uk                                             ##
##                                                                                          ##
## Address: The Weatherall Institute of Molecular Medicine                                  ##
##          University of Oxford                                                            ##
##          John Radcliffe Hospital                                                         ##
##          Headington                                                                      ##
##          Oxford OX3 9DS                                                                  ##
##                                                                                          ##
##############################################################################################

# R script creating a vocabulary file per tissue
# Date: 12/04/2016
# Author: Ron Schweßinger

# Use SasQ batch to calculate the SFR of every possible, non-ambiguous k-mer in the tissue of interest
# 
# Input: 
#   functions: R functions to source
#   tissue: tissue ID of interest
#   frag.type: DNase/ATAC fragmentation type
#   pnorm.tag: propensity norm identifier
#   kmer.list: infile listing all possible kmers
#   organism: organism of origin (query against what (human/mouse))
#   data.dir: data directory 
#   output.file: file to write output to
# 
# Returns:
#   vocabulary file listing the SFR for each possible k-mer


# 1) Get Arguments and source ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

source(args[1])  #common R functions

tissue <- args[2]  #get input table

organism <- args[3]

pnorm.tag <- args[4]

data.dir <- args[5]

frag.type<- args[6]

output.file <- args[7]

kmer.list <- args[8]

# Parameter for testing
# source("/home/ron/fusessh/Sasquatch_offline/Sasquatch/R_utility/functions_sasq_r_utility.R")
# kmer.list <- "/home/ron/fusessh/Sasquatch_offline/Sasquatch/data_processing_pipeline/kmers/Kmers_5_6_7_combined.txt"
# organism <- "human"
# tissue <- "ENCODE_Duke_Cerebellum_OC_merged"
# output.file <- "/home/ron/fusessh/sandbox/vocab_sandbox"
# pnorm.tag <- "h_ery_1"
# frag.type <- "DNase"
# data.dir <- "/home/ron/fusessh/database_assembly/idx_correct_assembly/human/DNase/"

# library(pbapply)

# Adapted Function Versions ====================================================================
GrepProfileVocab <- function(kmer, strand, k5p, k5m, k6p, k6m, k7p, k7m){
  # grep 250 bp surrounding kmer strand specific profile and normalize for 
  # total cuts in 250 bp window (ovocabulary version with preread in kmer tables)
  #
  # Args:
  #   kmer: input kmer
  #   strand: plus/minus
  #   k5p: kmer count 5 plus strand
  #   k5m: kmer count 5 plus strand
  #   k6p: kmer count 6 plus strand
  #   k6m: kmer count 6 plus strand
  #   k7p: kmer count 7 plus strand
  #   k7m: kmer count 7 plus strand
  #
  # Returns:
  #   normalized profile (relative cut frequency) and kmer occurence count
  
  #decode  ambivalent fasta code char into kmer list
  kmer.list <- DecodeKmer(kmer)  
  kl <- nchar(kmer) #get kmer length  
  count <- 0 #initialise count variable
  profile <- rep(0,(250+kl)) #initialise empty profile
  
  if(strand != "plus" && strand != "minus"){ stop("Select plus or minus strand") }
     
  for(km in kmer.list){  #for all entries in the splitted up kmer list

     if(nchar(kmer) == 7){
       match <- if(strand == "plus") subset(k7p, k7p[, 1] == km) else subset(k7m, k7m[, 1] == km)
     }else if(nchar(kmer) == 6){
       match <- if(strand == "plus") subset(k6p, k6p[, 1] == km) else subset(k6m, k6m[, 1] == km)
     }else if(nchar(kmer) == 5){
       match <- if(strand == "plus") subset(k5p, k5p[, 1] == km) else subset(k5m, k5m[, 1] == km)
     }else{
       stop("No valid kmer length")
     }
    
    count <- count + match[1,2]	#get second column as the count
    split <- match[1,-c(1,2)] #remove kmer and count from split
    split <- as.numeric(split) #convert to numeric
    profile <- profile + split[c(26:(length(split)-25))]	#get 250 bp surrounding kmer from 300 bp profile
  }
  
  #normalise to all recorded cuts in the 250 surrounding region 
  #--> relative cut frequency
  profile <- (profile / sum(profile))
  
  #make list for reporting
  newlist <- list( "profile"=profile, "count"=count ) #assemble return list
  return(newlist)
  
}

GetFootprintVocab <- function(kmer, 
                         frag.type, 
                         smooth=TRUE, 
                         smooth.bandwidth=5,
                         k5p, k5m, k6p, k6m, k7p, k7m){
  # Wrapper to retrieve merged, pruned, smoothed profile of kmer
  #
  # Args:
  #   kmer: input kmer 

  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #   smooth: flag if to smooth the profile (TRUE or FALSE)
  #   smooth.bandwidth: bandwidth to use for normal kernel smoothing 
  #   (default 5 which is fixed for most worklfows)
  #   k5p: kmer count 5 plus strand
  #   k5m: kmer count 5 plus strand
  #   k6p: kmer count 6 plus strand
  #   k6m: kmer count 6 plus strand
  #   k7p: kmer count 7 plus strand
  #   k7m: kmer count 7 plus strand
  #
  # Returns:
  #   merged, pruned (smoothed) profile and count of the kmer occurence
  
  #get kmer length
  kl=nchar(kmer)

  #grep strand specific profiles & counts
  l.plus <- GrepProfileVocab(kmer, strand="plus", k5p=k5p, k5m=k5p, k6p=k6p, k6m=k6m, k7p=k7p, k7m=k7m)
  l.minus <- GrepProfileVocab(kmer, strand="minus", k5p=k5p, k5m=k5m, k6p=k6p, k6m=k6m, k7p=k7p, k7m=k7m)
  
  #merge profiles accoring to selected fragmentation type
  if(frag.type == "ATAC"){
    #sum up profile as average of both profiles
    profile.merge <- (l.plus$profile + l.minus$profile)/2
  }else if(frag.type == "DNase"){
    #get length of the kmer as middle separately
    middle <- ( l.plus$profile[c(126:(125+kl))] + l.minus$profile[c(126:(125+kl))] ) / 2 #calc average of both strands for the kmer length middle
    #merge stradn specific with average of middle and only strand specific flanks
    profile.merge <- c( l.plus$profile[c(1:125)], middle, l.minus$profile[c((125+kl+1):length(l.minus$profile))] )	#combine left flank from plus right flank from minus and kmer middle from the average
  }else{
    cat("Select according framgentation type \"DNase\" or \"ATAC\"")
  }
  
  #smooth if specified so
  if(smooth){
    profile.merge <- ksmooth(c(1:length(profile.merge)), profile.merge, kernel="normal", bandwidth=smooth.bandwidth)$y
  }
  
  return(list(profile=profile.merge, count=l.plus$count)) #return
  
}

GetSFRVocab <- function(kmer, frag.type, k5p, k5m, k6p, k6m, k7p, k7m){
  # Wrapper function to get the SFR ratio
  # If indicated and available, use the present vocabulary file to directly grep the SFR
  # Else get the average profile, estimate the borders and calculate the SFR
  # Note that for using the vocabulary file only nonambivlent DNA chars are allowed
  # For the alternative ambivalent chars are decoded
  #
  # Args:
  #   kmer: k-mer to query
  #   k5p: kmer count 5 plus strand
  #   k5m: kmer count 5 plus strand
  #   k6p: kmer count 6 plus strand
  #   k6m: kmer count 6 plus strand
  #   k7p: kmer count 7 plus strand
  #   k7m: kmer count 7 plus strand
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #
  # Returns:
  #   SFR

    #no vocabfile present --> run retrieving, border estimation and SFR calculation
    
    fp <- GetFootprintVocab(kmer=kmer, frag.type=frag.type, 
                       smooth=TRUE, k5p=k5p, k5m=k5p, k6p=k6p, k6m=k6m, k7p=k7p, k7m=k7m) #get smoothed profile and count
    
    sh <- SobelBorders(profile=fp$profile, kl=nchar(kmer)) #estimate optimal shoulder midpoints and ranges
    
    #calculate SFR
    if(sh$flag == TRUE){
      sfr <- CalcSFR(
        profile=fp$profile, 
        us.mid=sh$us, ds.mid=sh$ds,
        range.us=sh$range.us, range.ds=sh$range.ds
      )
    }else{
      print("Shoulders could not be estimated properly will return SFR 1.0")
      sfr <- 1.0
    }
    return(sfr)
    
} 



# START ==================================================================================

# 1) get kmer.list and sort accordingly
tab <- read.table(kmer.list, colClasses = c("character"), header=FALSE)
colnames(tab) <- c("kmer")

k5p <- read.table(paste0(data.dir, "/", tissue, "/counts/", "kmers_5_count_", tissue, "_pnorm_", pnorm.tag, "_plus.txt"))
k5m <- read.table(paste0(data.dir, "/", tissue, "/counts/", "kmers_5_count_", tissue, "_pnorm_", pnorm.tag, "_minus.txt"))
k6p <- read.table(paste0(data.dir, "/", tissue, "/counts/", "kmers_6_count_", tissue, "_pnorm_", pnorm.tag, "_plus.txt"))
k6m <- read.table(paste0(data.dir, "/", tissue, "/counts/", "kmers_6_count_", tissue, "_pnorm_", pnorm.tag, "_minus.txt"))
k7p <- read.table(paste0(data.dir, "/", tissue, "/counts/", "kmers_7_count_", tissue, "_pnorm_", pnorm.tag, "_plus.txt"))
k7m <- read.table(paste0(data.dir, "/", tissue, "/counts/", "kmers_7_count_", tissue, "_pnorm_", pnorm.tag, "_minus.txt"))

# 2) calculate SFRs
tab$SFR <- sapply(tab$kmer, function(x){ 
    x <- round(GetSFRVocab(as.character(x), frag.type, k5p=k5p, k5m=k5m, k6p=k6p, k6m=k6m, k7p=k7p, k7m=k7m), digits=5)
    return(x)
  })

# 3) Assemble and write
write.table(tab, file = output.file, quote=F, sep="\t", col.names=F, row.names=F)

# ==============================================================================  
  
print("R-script for vocabulary file creation has run successfully")
