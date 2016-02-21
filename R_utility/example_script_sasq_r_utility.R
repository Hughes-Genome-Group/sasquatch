library(ggplot2)

source("/home/ron/fusessh/Sasquatch_offline/Sasquatch/R_utility/functions_sasq_r_utility.R")

data.dir <- "/home/ron/fusessh/database_assembly/idx_correct_assembly/human/DNase/"
background.dir <- "/home/ron/fusessh/database_assembly/idx_correct_assembly/background/"
background.tissue <- "hg18_human_JH60"
  
out.dir <- "/home/ron/Daten/WIMM/Sasquatch_working/"


#TODO(rschwess): Include example of how to retrieve single strand profiles

### === Set Some Parameters === ###
kmer <- "AGATAA"

tissue <- "human_erythroid_hg18"

frag.type <- "DNase"

smooth <- TRUE

# vocab.file <- paste0(data.dir,"/",tissue,"/vocabulary_",tissue,".txt")

### START ###

### ===== TEST R BASIC FUNCTIONS ===== ###
kmer <- "CACGTG"
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
kmer1 <- "WGATAA"
kmer2 <- "WGATTA"
fp1 <- GetFootprint(kmer=kmer1, tissue=tissue, data.dir=data.dir, frag.type=frag.type, smooth=T)
fp2 <- GetFootprint(kmer=kmer2, tissue=tissue, data.dir=data.dir, frag.type=frag.type, smooth=T)

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



#disect a longer sequence
seq <- "GGATATGATAGATACCT"

dl <- DissectSequence(seq, 7, list=TRUE)


### ===== TEST R WRAPPER FUNCTIONS ===== ###

library(RColorBrewer)
color.store <- brewer.pal(3,"Set1")

#wrapper to get SFR
sfr <- GetSFR(kmer="CACGTG", tissue="human_erythroid_hg18", data.dir=data.dir, vocab.flag=T, frag.type="DNase")

kmer = "AGATCA"

#single plot wrapper
p <- PlotSingleKmer(kmer=kmer, tissue=tissue, data.dir=data.dir, frag.type=frag.type, 
                    smooth=TRUE, plot.shoulders=F, ylim=c(0,0.01), xlim=c(-70,70),
                    color="red")
p

ggsave(p, filename=paste0("/home/ron/Dokumente/phd_application/interviews/presentation/figures/single_profile_", kmer, "_humane_erythroid.svg"), width=10, height=10/2.5)

#wrapper for overlap from kmers only
p <- PlotOverlapKmers(
  kmer1="CACGTG", kmer2="CACGTT", tissue1=tissue, tissue2=tissue, data.dir, frag.type="DNase", 
  smooth=T, ylim=c(0,0.01), xlim=c(-75,75), plot.shoulders = TRUE
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
        id=c("1", "2", "3"), 
        ref=c("ATAGATAATCGCT", "ATAGATAATCGCT", "ATATATTCTCGCT"),
        var=c("ATAGATCATCGCT", "ATAGATTATCGCT", "ATAGATGATCGCT")
        )

bcomp <- RefVarBatch(ref.var.df=tdf, kl=7, damage.mode="exhaustive", 
                        tissue, data.dir, vocab.flag=TRUE, frag.type=frag.type)


#jaspar batch
#load Rdata object storing the jaspar 2014 pwms (all versions)
library(Biostrings)
library(TFBSTools)

load("/home/ron/fusessh/database_assembly/jaspar/jaspar2014.human.9606.all.versions")

#single Jaspar query
QueryJaspar(sequence="AGATAATAG", threshold=0.8, pwm.data=pwm.in)

#warpper for batch quary a batch Ref Var Dataframe
#test with bcomp
jbcomp <- QueryJasparBatch(df=bcomp, damage.threshold=0.3, match.threshold=0.8, pwm.data=human.pwm)
  
#wrapper to compare two sequences
sequence1 <- "CAGTTTCATGAGG"
sequence2 <- "CAGTTTTATGAGG"

comp <- CompareSequences(
  sequence1=sequence1, 
  sequence2=sequence2, 
  kl=7,
  data.dir=data.dir,
  damage.mode="exhaustive",
  tissue="DNase_He_refined_LNCaP_50U_50_100bp_L_D", 
  vocab.flag=FALSE,
  frag.type="DNase", 
  plots="highest"
  )


#wrapper to get strand specific footprint profiles for tisue or background
sfp <- GetFootprintStrand(kmer="WGATAA", tissue=tissue, data.dir=data.dir, frag.type=frag.type, smooth=TRUE, smooth.bandwidth=5, background.flag=FALSE)

bfp <- GetFootprintStrand(kmer="WGATAA", tissue=background.tissue, data.dir=background.dir, frag.type=frag.type, smooth=TRUE, smooth.bandwidth=5, background.flag=TRUE)


splots <- PlotSingleStrands(kmer="WGATAA", tissue = background.tissue, data.dir = background.dir, frag.type = frag.type,
                            smooth=TRUE, background.flag = TRUE)



##############################################################

#insilico mutations

#get mutation data frame
d <- GetPossibleMutations(sequence=c("AGGGATACGTAGACGGTGTAAACCCGTGCATAGTAGA"), kl=7, chr="chrX", position=1345990)

d$damage <- apply(d, 1, function(x) CompareSequences(sequence1=x[5], sequence2=x[6], kl=7, damage.mode="exhaustive", 
                 tissue=tissue, data.dir=data.dir, vocab.flag=TRUE, 
                 vocab.file=vocab.file, frag.type=frag.type, plots=FALSE)$summary$total.damage )

#wrapper for mutating everything
d1 <- InSilicoMutation(  sequence="AGGGATACGTAGACGGGTGT", 
                                kl=7, 
                                chr="chr1",
                                position=13330000,
                                report="all",
                                damage.mode="exhaustive",
                                tissue=tissue,
                                data.dir=data.dir,
                                vocab.flag=TRUE,
                                frag.type=frag.type
                                )

rp <- RainbowPlot(d1)
rp


library(ggplot2)
library(RColorBrewer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)
genome <- BSgenome.Hsapiens.UCSC.hg18

seq <- as.character(getSeq(genome, "chr16", start=145846, end=145949))

r2.df <- InSilicoMutation(sequence=seq, 
                         kl=7, 
                         chr="chr16",
                         position=145852,
                         report="all",
                         damage.mode="exhaustive",
                         tissue="human_erythroid_hg18",
                         data.dir=data.dir,
                         vocab.flag=TRUE,
                         frag.type=frag.type
)

rp <- RainbowPlot(r2.df, ylim=c(-4,4)) + geom_vline(xintercept=145912, linetype="dashed")
rp

