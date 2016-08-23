##############################################################################################
##                                                                                          ##
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##                                                                                          ##
##############################################################################################

###############################
## Example R-Script          ##
## Author: Ron Schwessinger  ##
###############################

source("/home/ron/fusessh/Sasquatch_offline/Sasquatch/R_utility/functions_sasq_r_utility.R")

### === Set Some Initial Parameters ===============================================================

data.dir <- "/home/ron/fusessh/database_assembly/idx_correct_assembly/human/DNase/"
pnorm.tag <- "h_ery_1" #identifier for the propensity source used

# only for background plots
background.dir <- "/home/ron/fusessh/database_assembly/idx_correct_assembly/background/"
background.tissue <- "hg18_human_h_ery_1"

# output directory for tables and plots  
out.dir <- "/home/ron/Daten/WIMM/Sasquatch_working/"

# select tissue (e.g. "list.files(data.dir)")
tissue <- "WIMM_primary_erythroid_Fibach_Fade8"

# select fragmentation type
frag.type <- "DNase"


### ===== TEST R BASIC FUNCTIONS ==================================================================


# single k-mers analysis  ---------------------------------------------------------------------

# select a k-mer of interest
kmer <- "CGCATGC"

# get the footprint 
fp <- GetFootprint(kmer=kmer, tissue=tissue, data.dir=data.dir, pnorm.tag=pnorm.tag, frag.type=frag.type, smooth=TRUE)
#returns list object with $profile and $count


# estimate the shoulders from the profile (use smoothed profile or smooth within call!)
sh <- SobelBorders(fp$profile, kl=nchar(kmer))
# returns list object: 
#   $us and $ds shoulder positions upstream and downstream respectively
#   $range.us $range.ds range(size) of the respective shoulder; 
#   $flag TRUE/FALSE indicating if shoulders could be estimated 


# make single, merged profile plot 
p <- PlotSingle(profile=fp$profile, 
                kl=nchar(kmer), 
                plot.shoulders=TRUE, 
                shoulders = sh, 
                ylim=c(0,0.0123))
plot(p)

# further example: make a pruned profile plot with no shoulders
p <- PlotSingle(profile=fp$profile, 
                kl=nchar(kmer), 
                plot.shoulders=FALSE, 
                ylim=c(0,0.01), 
                xlim=c(-50,50))
plot(p)

# Get overlap of profiles ----------------------------------------------------------------------

kmer1 <- "WGATAA" #note FASTA ambiguous code is supported
kmer2 <- "WGATTA"

fp1 <- GetFootprint(kmer=kmer1, tissue=tissue, data.dir=data.dir, pnorm.tag=pnorm.tag, frag.type=frag.type, smooth=T)
fp2 <- GetFootprint(kmer=kmer2, tissue=tissue, data.dir=data.dir, pnorm.tag=pnorm.tag, frag.type=frag.type, smooth=T)

# make an overlap plot
p <- PlotOverlap(
  fp1$profile, 
  fp2$profile, 
  kmer1, 
  kmer2,
  fp1$count,
  fp2$count,
  ymode="separate", 
  xlim=c(-80,80), 
  ylim=c(0.0020, 0.0100),
  plot.shoulders=FALSE
  )

p

# dissect a longer sequence  ----------------------------------------------------------------------
seq <- "GGATATGATAGATACCT"

# helper function to dissect sequence into list of 1 bp sliding k-mers
dl <- DissectSequence(seq, kl=7, list=FALSE)
# returns: list (list=TRUE) or vector (list = FALSE)
# use with lapply or sapply how you like

### ===== TEST R WRAPPER FUNCTIONS =================================================================

# library(RColorBrewer)
# color.store <- brewer.pal(3,"Set1")

# Wrapper to get SFR from k-mer and tissues --> returns single SFR value --------------------------
sfr <- GetSFR(kmer="CACGTG", 
              tissue=tissue, 
              data.dir=data.dir,
              pnorm.tag=pnorm.tag,
              vocab.flag=F, 
              frag.type="DNase")

# Wrapper for single plot
p <- PlotSingleKmer(kmer="GGCGGG", tissue=tissue, data.dir=data.dir, pnorm.tag=pnorm.tag, frag.type=frag.type, 
                    smooth=FALSE, plot.shoulders=FALSE, ylim=c(0,0.01), xlim=c(-70,70),
                    color="black")
p

# Example to save a plot --------------------------------------------------------------------------
ggsave(p, filename=file.path(out.dir,
       paste0("single_profile_", kmer, "_humane_erythroid.svg")), 
       width=10, height=10/2.5)

# Wrapper for overlap plots from k-mers -----------------------------------------------------------
p <- PlotOverlapKmers(
  kmer1="CACGTG", kmer2="CACGTT", 
  tissue1=tissue, tissue2=tissue, 
  data.dir=data.dir, pnorm.tag=pnorm.tag, frag.type="DNase", 
  smooth=TRUE, plot.shoulders = TRUE,
  ylim=c(0,0.01), xlim=c(-75,75)
  )
p

# Wrapper to query longer sequence ----------------------------------------------------------------
dl <- QueryLongSequence(sequence=seq, 
                        kl=7, 
                        tissue=tissue, 
                        data.dir=data.dir,
                        pnorm.tag=pnorm.tag,
                        vocab.flag=TRUE, 
                        frag.type=frag.type, 
                        plots=FALSE, 
                        smooth=TRUE, 
                        plot.shoulders=TRUE, 
                        ylim=c(0,0.01), 
                        xlim=c(-125,125)
                        )

# Wrapper for Ref-Var Batch query -----------------------------------------------------------------
# make example dataframe
tdf <- data.frame(
        id=c("1", "2", "3"), 
        ref=c("ATAGATAATCGCT", "ATAGATAATCGCT", "ATATATTCTCGCT"),
        var=c("ATAGATCATCGCT", "ATAGATTATCGCT", "ATAGATGATCGCT")
        )

bcomp <- RefVarBatch(ref.var.df=tdf, 
                     kl=7, 
                     damage.mode="exhaustive", 
                     tissue=tissue, 
                     data.dir=data.dir,
                     pnorm.tag=pnorm.tag,
                     vocab.flag=TRUE, 
                     frag.type=frag.type)

# Meet old JASPAR -----------------------------------------------------------------------------------
#load Rdata object storing the jaspar 2016 pwms (all versions)
library(Biostrings)
library(TFBSTools)

load("/home/ron/fusessh/database_assembly/jaspar/jaspar2016.human.9606.all.versions")

# Single JASPAR query
QueryJaspar(sequence="AGATAATAG", threshold=0.8, pwm.data=human.pwm)

# Wrapper for batch quary a batch Ref Var Dataframe
jbcomp <- QueryJasparBatch(df=bcomp, damage.threshold=0.3, match.threshold=0.8, pwm.data=human.pwm)
  

# Wrapper to compare two sequences -------------------------------------------------------------------
comp <- CompareSequences(
  sequence1="CAGTTTTATGAGG", 
  sequence2="CAGTTTCATGAGG", 
  kl=7,
  data.dir = data.dir,
  pnorm.tag = pnorm.tag,
  damage.mode = "exhaustive",
  tissue = tissue, 
  vocab.flag = TRUE,
  frag.type = "DNase", 
  plots = "highest"
  )

# Wrapper to get strand specific footprint profiles for tissue or background -------------------------
sfp <- GetFootprintStrand(kmer="WGATAA", tissue=tissue, data.dir=data.dir, pnorm.tag = pnorm.tag, frag.type=frag.type, smooth=TRUE, smooth.bandwidth=5, background.flag=FALSE)

bfp <- GetFootprintStrand(kmer="WGATAA", tissue=background.tissue, data.dir=background.dir, pnorm.tag = pnorm.tag, frag.type=frag.type, smooth=TRUE, smooth.bandwidth=5, background.flag=TRUE)

# plot single strands
splots <- PlotSingleStrands(kmer="WGATAA", tissue = background.tissue, data.dir = background.dir, frag.type = frag.type,
                            smooth=TRUE, background.flag = TRUE)

# Insilico mutations ---------------------------------------------------------------------------------

# will split the sequence into windows matching to the selected k-mer length kl
# e.g. for kl=7 it wil split the sequence into 1 bp sliding windows of 13 bp length

# !!!Note!!! The first base position where the damage is predcited is the ("kl"th) position
# in the sequence. E.g. for kl=7 the position value should refer to the 7th base in the sequence. 
# Vice Versa the sequence input should start kl-1 bp before your base position of interest and the 
# last kl-1 bp positions will not be analysed explicitly

df.insilico <- InSilicoMutation(sequence="GTGCCCGCATGTGCTTATTTCTGCAAAAATAAACCATGGCAGG", 
                                kl=7, 
                                chr="chr1",
                                position=13330000,
                                report="all",
                                damage.mode="exhaustive",
                                tissue=tissue,
                                data.dir=data.dir,
                                pnorm.tag = pnorm.tag,
                                vocab.flag=TRUE,
                                frag.type=frag.type,
                                progress.bar=TRUE
                                )
head(df.insilico)

# Note: progress.bar = TRUE will require the package "pbapply" which visualized the progress
# install packagew with install.packages("pbapply") or set to FALSE

# Make a InSilicoMutationplot from the processed in silico mutation data.frame
rp <- InSilicoMutationPlot(df.insilico)
rp

# Full example for in silico mutation data frame and InSilicoMutationplot --------------
library(ggplot2)
library(RColorBrewer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)
genome <- BSgenome.Hsapiens.UCSC.hg18

# To get the next 30 bases starting from start.pos
chr <- "chr16"
start.pos <- 145852
end.pos <- start.pos + 30

# Get sequence
seq <- as.character(getSeq(genome, "chr16", start=start.pos-6, end=end.pos+6))

df.insilico <- InSilicoMutation(sequence=seq, 
                         kl=7, 
                         chr="chr16",
                         position=start.pos,
                         report="all",
                         damage.mode="exhaustive",
                         tissue="human_erythroid_hg18",
                         data.dir=data.dir,
                         pnorm.tag = pnorm.tag,
                         vocab.flag=TRUE,
                         frag.type=frag.type,
                         progress.bar = TRUE
                         )

rp <- InSilicoMutationPlot(df.insilico, ylim=c(-4,4))
rp

# Manual alternative --------------------------------------------------------------------------------

# Get mutation data frame
d <- GetPossibleMutations(sequence=c("AGGGATACGTAGACGGTGTAA"), kl=7, chr="chrX", position=1345990)

# calculate damage using apply and the more basic functions
d$damage <- apply(d, 1, function(x) CompareSequences(sequence1=x[5], 
                                                     sequence2=x[6], 
                                                     kl=7, 
                                                     damage.mode="exhaustive", 
                                                     tissue=tissue, 
                                                     data.dir=data.dir,
                                                     pnorm.tag = pnorm.tag,
                                                     vocab.flag=TRUE, 
                                                     frag.type=frag.type, 
                                                     plots=FALSE
                                                     )$summary$total.damage
                  )

# Preload data for faster processing --------------------------------------------------------------------------------

# Preload vocab file
data.dir <- "/home/ron/fusessh/database_assembly/idx_correct_assembly/human/DNase/"
pnorm.tag <- "h_ery_1" #identifier for the propensity source used

# output directory for tables and plots  
out.dir <- "/home/ron/Daten/WIMM/Sasquatch_working/"

# select tissue (e.g. "list.files(data.dir)")
tissue <- "WIMM_primary_erythroid_Fibach_Fade8"

# select fragmentation type
frag.type <- "DNase"

vocab <- PreLoadVocab(data.dir, tissue)
profiles <- PreLoadKmerProfiles(7, data.dir, tissue, pnorm.tag)

# Apply functions with preload vocabulary file
# 1) Get Footprint
fp <- GetFootprint("CACGTGG", tissue, data.dir, pnorm.tag, frag.type, smooth=T, preload=T, preload.profiles = profiles)

# 2) GetSFR
sfr.v <- GetSFR(kmer="CGCATGC", tissue=tissue, data.dir=data.dir, pnorm.tag=pnorm.tag, vocab.flag=T, frag.type="DNase", preload=T, preload.vocab=vocab, preload.profiles=profiles)
sfr.p <- GetSFR(kmer="CGCATGC", tissue=tissue, data.dir=data.dir, pnorm.tag=pnorm.tag, vocab.flag=F, frag.type="DNase", preload=T, preload.vocab=vocab, preload.profiles=profiles)
sfr.v
sfr.p

# 3) Dissect lOnger Sequence
ds1 <- QueryLongSequence("CCGCGCTTATGTACC", 7, tissue, data.dir, pnorm.tag, vocab.flag = F, frag.type = "DNase", preload=F)
ds2 <- QueryLongSequence("CCGCGCTTATGTACC", 7, tissue, data.dir, pnorm.tag, vocab.flag = F, frag.type = "DNase", preload=T, preload.vocab = vocab, preload.profiles = profiles)
ds1
ds2

#compare time
system.time(QueryLongSequence("CCGCGCTTATGTACC", 7, tissue, data.dir, pnorm.tag, vocab.flag = F, frag.type = "DNase", preload=F))
system.time(QueryLongSequence("CCGCGCTTATGTACC", 7, tissue, data.dir, pnorm.tag, vocab.flag = F, frag.type = "DNase", preload=T, preload.vocab = vocab, preload.profiles = profiles))

system.time(QueryLongSequence("CCGCGCTTATGTACC", 7, tissue, data.dir, pnorm.tag, vocab.flag = T, frag.type = "DNase", preload=F))
system.time(QueryLongSequence("CCGCGCTTATGTACC", 7, tissue, data.dir, pnorm.tag, vocab.flag = T, frag.type = "DNase", preload=T, preload.vocab = vocab, preload.profiles = profiles))

# 4) Compare Sequences
comp <- CompareSequences(
  sequence1="CAGTTTCATGAGG", 
  sequence2="CAGTTTTATGAGG", 
  kl=7,
  data.dir=data.dir,
  pnorm.tag = pnorm.tag,
  damage.mode="exhaustive",
  tissue=tissue, 
  vocab.flag=T,
  frag.type="DNase", 
  plots="highest",
  preload=T,
  preload.vocab=vocab,
  preload.profiles=profiles
)

# 5) RefVarBatch
tdf <- data.frame(
  id=c("1", "2", "3"), 
  ref=c("ATAGATAATCGCT", "ATAGATAATCGCT", "ATATATTCTCGCT"),
  var=c("ATAGATCATCGCT", "ATAGATTATCGCT", "ATAGATGATCGCT")
)

bcomp <- RefVarBatch(ref.var.df=tdf, 
                     kl=7, 
                     damage.mode="exhaustive", 
                     tissue=tissue, 
                     data.dir=data.dir,
                     pnorm.tag=pnorm.tag,
                     vocab.flag=TRUE, 
                     frag.type=frag.type,
                     preload=TRUE,
                     preload.vocab = vocab,
                     preload.profiles = profiles)

# 6) In Silico Mutation
df.insilico <- InSilicoMutation(sequence="GTGCCCGCATGTGCTTATTTCTGCAAAAATAAACCATGGCAGG", 
                                kl=7, 
                                chr="chr1",
                                position=13330000,
                                report="all",
                                damage.mode="exhaustive",
                                tissue=tissue,
                                data.dir=data.dir,
                                pnorm.tag = pnorm.tag,
                                vocab.flag=TRUE,
                                frag.type=frag.type,
                                progress.bar=TRUE,
                                preload=TRUE,
                                preload.vocab = vocab,
                                preload.profiles = profiles
)
head(df.insilico)



p <- InSilicoMutationPlot(df.insilico)
p
