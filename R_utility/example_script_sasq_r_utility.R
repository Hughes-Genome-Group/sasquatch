library(ggplot2)

source("/home/ron/fusessh/Sasquatch_offline/Sasquatch/R_utility/functions_sasq_r_utility.R")

data.dir <- "/home/ron/fusessh/database_assembly/idx_correct_assembly/human/DNase/"
background.dir <- "/home/ron/fusessh/database_assembly/idx_correct_assembly/background/"
background.tissue <- "hg18_human_JH60"
  
out.dir <- "/home/ron/Daten/WIMM/Sasquatch_working/temp_out"


#TODO(rschwess): Include example of how to retrieve single strand profiles

### === Set Some Parameters === ###
kmer <- "CACGTG"

tissue <- "human_erythroid_hg18"

frag.type <- "DNase"

smooth <- TRUE

vocab.file <- paste0(data.dir,"/",tissue,"/vocabulary_",tissue,".txt")

### START ###

### ===== TEST R BASIC FUNCTIONS ===== ###

#get the footprint
fp <- GetFootprint(kmer=kmer, tissue=tissue, data.dir=data.dir, frag.type=frag.type, smooth=F)

#get shoulder list (use smoothed profile or smooth within call)
sh <- SobelBorders(fp$profile, kl=nchar(kmer))
# sh <- SobelBorders(SmoothProfile(fp$profile, 5), kl=nchar(kmer))

#make single profile plot
p <- PlotSingle(profile=fp$profile, kl=nchar(kmer), plot.shoulders=TRUE, shoulders=sh, ylim=c(0,0.01))
p

#make pruned profile plot
p <- PlotSingle(profile=fp$profile, kl=nchar(kmer), plot.shoulders=F, shoulders=sh, ylim=c(0,0.01), xlim=c(-50,50))
p

#plot overlap
kmer1 <- "AGATAA"
kmer2 <- "GGATAA"
fp1 <- GetFootprint(kmer="WGATAA", tissue=tissue, data.dir=data.dir, frag.type=frag.type, smooth=T)
fp2 <- GetFootprint(kmer=kmer2, tissue=tissue, data.dir=data.dir, frag.type=frag.type, smooth=T)

p <- PlotOverlap(fp1$profile, fp2$profile, kmer1, kmer2)
p


#disect a longer sequence
seq <- "GGATATGATAGATACCT"

dl <- DissectSequence(seq, 7, list=TRUE)


### ===== TEST R WRAPPER FUNCTIONS ===== ###

#wrapper to get SFR
sfr <- GetSFR(kmer, tissue, data.dir, vocab.flag=TRUE, vocab.file=vocab.file, frag.type="DNase")

#single plot wrapper
p <- PlotSingleKmer(kmer="TGATAA", tissue=tissue, data.dir=data.dir, frag.type=frag.type, 
                    smooth=TRUE, plot.shoulders=F, ylim=c(0,0.01), xlim=c(-75,75))
p

#wrapper for overlap from kmers only
p <- PlotOverlapKmers(
  kmer1="CACGTG", kmer2="CACGTT", tissue1=tissue, tissue2=tissue, data.dir, frag.type="DNase", 
  smooth=F, ylim=c(0,0.01), xlim=c(-75,75)
  )
p

#wrapper to query longer sequence
dl <- QueryLongSequence(sequence=seq, kl=7, tissue=tissue, data.dir=data.dir, 
                        vocab.flag=TRUE, vocab.file=vocab.file, frag.type=frag.type, 
                        plots=TRUE, smooth=smooth, plot.shoulders=TRUE, ylim=c(0,0.01), xlim=c(-125,125)
                        )

#wrapper for Ref-Var Batch query

#make example dataframe
tdf <- data.frame(
        id=c(1,2,3), 
        ref=c("ATAGATAATCGCT", "ATAGATAATCGCT", "ATAGATAATCGCT"),
        var=c("ATAGATCATCGCT", "ATAGATTATCGCT", "ATAGATGATCGCT")
        )

bcomp <- RefVarBatch(ref.var.df=tdf, kl=7, damage.mode="exhaustive", 
                        tissue, data.dir, vocab.flag=TRUE, vocab.file=vocab.file, frag.type=frag.type)


#jaspar batch


#wrapper to compare two sequences
sequence1 <- "ATAGATAATCGCT"
sequence2 <- "ATAGATCATCGCT"

comp <- CompareSequences(sequence1=sequence1, sequence2=sequence2, kl=6, damage.mode="exhaustive", 
                             tissue=tissue, data.dir=data.dir, vocab.flag=TRUE, 
                             vocab.file=vocab.file, frag.type=frag.type, plots=FALSE)


#wrapper to get strand specific footprint profiles for tisue or background
sfp <- GetFootprintStrand(kmer="WGATAA", tissue=tissue, data.dir=data.dir, frag.type=frag.type, smooth=TRUE, smooth.bandwidth=5, background.flag=FALSE)

bfp <- GetFootprintStrand(kmer="WGATAA", tissue=background.tissue, data.dir=background.dir, frag.type=frag.type, smooth=TRUE, smooth.bandwidth=5, background.flag=TRUE)


splots <- PlotSingleStrands(kmer="WGATAA", tissue = background.tissue, data.dir = background.dir, frag.type = frag.type,
                            smooth=TRUE, background.flag = TRUE)



#

