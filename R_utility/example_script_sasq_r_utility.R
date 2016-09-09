##############################################################################################
##                                                                                          ##
## Sasquatch, Sequence based prediction of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##                                                                                          ##
##############################################################################################

###############################
## Example R-Script          ##
## Author: Ron Schwessinger  ##
###############################

source("./R_utility/functions_sasq_r_utility.R")

# Note: the example script runs with the pre-processed tissue dummy data provided with the code distribution. 
# To run different sequences on different tissues, please download your data of interest and the appropriate 
# background data from the webtool address: http://apps.molbiol.ox.ac.uk/sasquatch/cgi-bin/foot.cgi
# and extract them into the data repository 

### === Set Some Initial Parameters ===========================================

data.dir <- "./data/human/DNase/"  # data storage where downloaded / pre-processed data were extracted 

pnorm.tag <- "h_ery_1" # identifier for the propensity source used ["h_ery_1" = human, "m_ery_1" = mouse]

# select tissue (e.g. "list.files(data.dir)")
tissue <- "Dummy_tissue_example"

# select fragmentation type ["DNase" or "ATAC" (currently only for testing purposes)]
frag.type <- "DNase"

# Paths to background data only needs specifying when tryin to plot background plot. Every other normalisation is already done. 
background.dir <- "./data/human/background/"
background.tissue <- "Background_dummy_h_ery_1"

# output directory for table and plots  
out.dir <- "../sasq_sandbox/"


### ===== TEST SasQ R BASIC FUNCTIONS =========================================

# We will first run through the Basic SasQ functions, 
# bare in mind that there are wrappers for all common tasks, discussed afterwards. 

# 1) single k-mers analysis  --------------------------------------------------

# select a k-mer of interest
kmer <- "CGCATGC"

# get the footprint 
fp <- GetFootprint(kmer=kmer, tissue=tissue, data.dir=data.dir, pnorm.tag=pnorm.tag, frag.type=frag.type, smooth=TRUE)
# returns list object with $profile and $count

# estimate the shoulders from the profile (use a smoothed profile or smooth within call!)
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

# further example: make a pruned profile plot with no shoulders plotted
p <- PlotSingle(profile=fp$profile, 
                kl=nchar(kmer), 
                plot.shoulders=FALSE, 
                ylim=c(0,0.01), 
                xlim=c(-50,50))
plot(p)


# 2) Get overlap of profiles --------------------------------------------------

kmer1 <- "TGACTCA"
kmer2 <- "TGAGTCA"

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
  ymode="separate"
  )

plot(p)


# 3) dissect a longer sequence  ----------------------------------------------------------------------
seq <- "GGATATGATAGATACCT"
# helper function to dissect sequence into list of 1 bp sliding k-mers
dl <- DissectSequence(seq, kl=7, list=FALSE)
# returns: list (list=TRUE) or vector (list = FALSE)
# use with lapply or sapply or ever you how you like


### ===== TEST R WRAPPER FUNCTIONS =================================================================

# color.store <- brewer.pal(3,"Set1")

# vocab.flag = TRUE
# indicates that we have a precalculated vocabulary file (SFR for every unique k-mer) present in the repository subdirectory
# this is true for most all the tissues we distribute but can be switched of in case you haven't processed that file yet.
# using the vocabulary file speeds the calculation up significantly. 

# Wrapper to get SFR from k-mer and tissues --> returns single SFR value --------------------------
sfr <- GetSFR(kmer="CACGTG", 
              tissue=tissue, 
              data.dir=data.dir,
              pnorm.tag=pnorm.tag,
              vocab.flag=TRUE,  
              frag.type="DNase")

# Wrapper for single plot
s <- PlotSingleKmer(kmer="GGCGGG", tissue=tissue, data.dir=data.dir, pnorm.tag=pnorm.tag, frag.type=frag.type, 
                    smooth=FALSE, plot.shoulders=FALSE, ylim=c(0,0.01), xlim=c(-70,70),
                    color="black")
plot(s)

# Example to save a plot --------------------------------------------------------------------------
ggsave(p, filename=file.path(out.dir, "single_profile_GGCGGG_humane_erythroid.png"), width=10, height=10/2.5)


# Wrapper for overlap plots from k-mers -----------------------------------------------------------
o <- PlotOverlapKmers(
  kmer1="TGACTCA", kmer2="TGAGTCA",
  tissue1=tissue, tissue2=tissue, 
  data.dir=data.dir, pnorm.tag=pnorm.tag, frag.type="DNase"
  )
plot(o)

# Wrapper to query longer sequence ----------------------------------------------------------------
dl <- QueryLongSequence(sequence="CACGTGG", 
                        kl=6, 
                        tissue=tissue, 
                        data.dir=data.dir,
                        pnorm.tag=pnorm.tag,
                        vocab.flag=TRUE, 
                        frag.type=frag.type
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
# load Rdata object storing the jaspar 2016 pwms (all versions)
library(Biostrings)
library(TFBSTools)

load("data/jaspar/jaspar2016.human.9606.all.versions")

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
  plots = FALSE
  )

# Wrapper to get strand specific footprint profiles for tissue or background -------------------------
sfp <- GetFootprintStrand(kmer="CACGTG", tissue=tissue, data.dir=data.dir, pnorm.tag = pnorm.tag, frag.type=frag.type, smooth=TRUE, smooth.bandwidth=5, background.flag=FALSE)

bfp <- GetFootprintStrand(kmer="CACGTG", tissue=background.tissue, data.dir=background.dir, pnorm.tag = pnorm.tag, frag.type=frag.type, smooth=TRUE, smooth.bandwidth=5, background.flag=TRUE)

# plot single strands (of background)
splots <- PlotSingleStrands(kmer="CACGTG", tissue = tissue, data.dir = data.dir, frag.type = frag.type, pnorm.tag = pnorm.tag,
                            smooth=TRUE, background.flag = FALSE)
bplots <- PlotSingleStrands(kmer="CACGTG", tissue = background.tissue, data.dir = background.dir, frag.type = frag.type,
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
# installthe package with install.packages("pbapply") or set to FALSE

# Make a InSilicoMutationplot from the processed in silico mutation data.frame
rp <- InSilicoMutationPlot(df.insilico)
plot(rp)

# Full example for in silico mutation data frame and InSilicoMutationplot --------------
# to get a sequence of interest load a reference genome
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)
genome <- BSgenome.Hsapiens.UCSC.hg18

# Fet 30 bases starting from start.pos of interest
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
                         tissue=tissue,
                         data.dir=data.dir,
                         pnorm.tag = pnorm.tag,
                         vocab.flag=TRUE,
                         frag.type=frag.type,
                         progress.bar = TRUE
                         )

rp <- InSilicoMutationPlot(df.insilico, ylim=c(-4,4))
plot(rp)


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
# We can preload the whole vocabulary file or kmer based profiles files into memory.
# This speeds up analysis over multiple k-mers, longer sequences, batches of variants and in silico mutations significantly. 
# We then provide the preloaded data to the respective functions and indicate what we have preloaded and provided. 
# Note: Different functions require different data to be preloaded. While running over many k-mers or longer sequences can be 
# done nicely with only the vocabulary file preloaded, everythin that tries to plto profiles will need the profiles loaded to 
# profit from the speed-up.

# To preload the vocabulary file
vocab <- PreLoadVocab(data.dir, tissue)

# To preload the cut profiles (specify which kmer size you want to preload)
profiles.6mers <- PreLoadKmerProfiles(6, data.dir, tissue, pnorm.tag)
profiles.7mers <- PreLoadKmerProfiles(7, data.dir, tissue, pnorm.tag)

# Apply functions with preload vocabulary file
# 1) Get Footprint
fp <- GetFootprint("CACGTG", tissue, data.dir, pnorm.tag, frag.type, smooth=T, preload=T, preload.profiles = profiles.6mers)

# 2) GetSFR from vocabulary directly or get the profiles and calculate it on the fly
sfr.v <- GetSFR(kmer="CGCATGC", tissue=tissue, data.dir=data.dir, pnorm.tag=pnorm.tag, vocab.flag=T, frag.type="DNase", preload=T, preload.vocab=vocab)
sfr.p <- GetSFR(kmer="CGCATGC", tissue=tissue, data.dir=data.dir, pnorm.tag=pnorm.tag, vocab.flag=F, frag.type="DNase", preload=T, preload.profiles=profiles)
sfr.v
sfr.p

# 3) Dissect longer Sequence
ds <- QueryLongSequence("CCGCGCTTATGTACC", 7, tissue, data.dir, pnorm.tag, vocab.flag = TRUE, frag.type = "DNase", preload=T, preload.vocab = vocab)
ds


#compare time
system.time(QueryLongSequence("CCGCGCTTATGTACC", 7, tissue, data.dir, pnorm.tag, vocab.flag = T, frag.type = "DNase", preload=F))
system.time(QueryLongSequence("CCGCGCTTATGTACC", 7, tissue, data.dir, pnorm.tag, vocab.flag = T, frag.type = "DNase", preload=T, preload.vocab = vocab))

# 4) Compare Sequences
comp <- CompareSequences(
  sequence1="CAGTTTCATGAGG", 
  sequence2="CAGTTTTATGAGG", 
  kl=7,
  data.dir=data.dir,
  pnorm.tag = pnorm.tag,
  damage.mode="exhaustive",
  tissue=tissue, 
  vocab.flag=TRUE,
  frag.type="DNase", 
  plots=FALSE,
  preload=TRUE,
  preload.vocab=vocab
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
                     preload.vocab = vocab)

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
                                preload.vocab = vocab
)
head(df.insilico)

p <- InSilicoMutationPlot(df.insilico)
plot(p)
