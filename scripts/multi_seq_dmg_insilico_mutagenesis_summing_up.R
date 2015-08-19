#Compare SFR values from multiple sequences given an expected number of footprint contributing kmers 

#get Arguments
args <- commandArgs(trailingOnly = TRUE)

intable <- read.table(args[1], header=FALSE)	#fsr per kmer position along the sequence one sequence per line

exp.fp.kmers <- 1 #expected number of footprint contributing kmers

kl <- as.numeric(args[2])

total.id <- as.numeric(args[3]) #total id running through complete seqs input

#temp input for testing
#intable <- read.table('/hts/data4/rschwess/dnase_motif_tissue/idx_duke_testout/merge_strands/multiseq_fsr_table.txt', header=FALSE)
#exp.fp.kmers <- 1
#kl=6
#total.id=1

#store seq-id identifiers 
seq.id <- as.data.frame(intable[,c(1:2)])
colnames(seq.id) <- c("id", "sequence")

#trim fsr only table
fsr.table <- intable[,-c(1:2)]
fsr.table <- as.data.frame(fsr.table)

#define reference sequence (firstline)
ref.fsr <- as.data.frame(fsr.table[1,])
var.fsr <- as.data.frame(fsr.table[2,])

#calculate dmg
dmg <- ref.fsr - var.fsr
sumof <- sum(dmg)


#print sum of dmg dmg 
cat(as.numeric(sumof))


