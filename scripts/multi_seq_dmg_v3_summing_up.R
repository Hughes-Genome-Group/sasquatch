#Compare SFR values from multiple sequences given an expected number of footprint contributing kmers 

#get Arguments
args <- commandArgs(trailingOnly = TRUE)

intable <- read.table(args[1], header=FALSE)	#fsr per kmer position along the sequence one sequence per line

exp.fp.kmers <- 1 #expected number of footprint contributing kmers

kl <- as.numeric(args[2])

total.id <- as.numeric(args[3]) #total id running through complete seqs input

#temp input for testing
#intable <- read.table('/t1-data1/WTSA_Dev/rschwess/dnase_motif_tissue/idx_duke_testout/merge_strands/multiseq_fsr_table.txt', header=FALSE)
#exp.fp.kmers <- 1
#kl=6
#total.id=1

#store seq-id identifiers 
seq.id <- as.data.frame(intable[,c(1:2)])
colnames(seq.id) <- c("id", "sequence")

#trim fsr only table
fsr.table <- intable[,-c(1:2)]
fsr.table <- as.data.frame(fsr.table)
colnames(fsr.table) <-c(1:dim(fsr.table)[2])	

#define reference sequence (firstline)
ref.seq <- as.character(seq.id[1,2])
ref.fsr <- as.data.frame(fsr.table[1,])
colnames(ref.fsr) <-c(1:dim(ref.fsr)[2])
#split up ref seq in kmers of length kl, store reference kmers to focus
temp.array <- strsplit(ref.seq, split="")[[1]]
ref.array <- c()
for(i in c(1:(length(temp.array)-(kl-1))) ){
		ref.array <- c( ref.array, substring(ref.seq, i, (i+kl-1)) )
}

#define reference sequence (firstline)
var.seq <- as.character(seq.id[2,2])
var.fsr <- as.data.frame(fsr.table[2,])
colnames(var.fsr) <-c(1:dim(var.fsr)[2])
#split up ref seq in kmers of length kl, store reference kmers to focus
temp.array <- strsplit(var.seq, split="")[[1]]
var.array <- c()
for(i in c(1:(length(temp.array)-(kl-1))) ){
		var.array <- c( var.array, substring(var.seq, i, (i+kl-1)) )
}

###use expected number of fp kmers to rank and fix the kmers to consider for calculation

#calculate dmg
dmg <- ref.fsr - var.fsr
sumof <- sum(dmg)

#rank kmers according to fsr
ref.tag <- "int"
var.tag <- "int"
temp.i <- 0
while(ref.tag == var.tag){	#go deacreasingly down until kmers carry variance / if no variance report last with 0 0  values
	temp.i <- temp.i + 1
	chosen.ref.kmers <- sort(abs(dmg), decreasing=TRUE)[temp.i]
	chosen.ids <- as.numeric(colnames(chosen.ref.kmers)) #ids to fix

	ref.tag <- unlist( lapply(chosen.ids, function(x) paste0(ref.array[x])) )	#count up kmers as kmer id
	var.tag <- unlist( lapply(chosen.ids, function(x) paste0(var.array[x])) )	#count up kmers as kmer id

	if(temp.i == length(ref.fsr)) { break;}

}

#for each sequence, for each chosen kmer
#calculate the difference in fsr
diff.tab <- data.frame( fsr.table[,chosen.ids] )
colnames(diff.tab) <- chosen.ids

for(i in c(1:dim(diff.tab)[2])){
	temp <- diff.tab[1,i]
	diff.tab[,i] <- temp - diff.tab[,i]
}


out.tab <- data.frame(diff.tab)		#for out table
# out.tab[,(dim(out.tab)[2] + 1)] <- apply(diff.tab,1,function(x) x <- sum(x)/exp.fp.kmers )	#calulcate dmg per focused kmer
# out.tab[,(dim(out.tab)[2] + 1)] <- apply(diff.tab,1,function(x) x <- sum(x))	#calculate total dmg per sequence
out.tab <- cbind(seq.id,out.tab)	#add id and sequences

#get focused fsrs
ref.fsr <- fsr.table[1,chosen.ids]
var.fsr <- fsr.table[2,chosen.ids]

#print Ref Var 
cat( paste( as.character(out.tab$sequence[1]), as.character(out.tab$sequence[2]), as.character(ref.tag), as.character(var.tag), as.numeric(ref.fsr), as.numeric(var.fsr), as.numeric(sumof), sep="\t"))


