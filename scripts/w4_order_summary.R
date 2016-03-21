#read in summary tab and sort by total dmg and print

#get Arguments
args <- commandArgs(trailingOnly = TRUE)

intable <- read.table(args[1], header=TRUE)	#fsr per kmer position along the sequence one sequence per line

out <- args[2]

intable <- intable[order(intable[,dim(intable)[2]], decreasing=TRUE),] 

write.table(intable, file=out, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
