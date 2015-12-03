#get Arguments
args <- commandArgs(trailingOnly = TRUE)

table1 <- read.table(args[1], header=TRUE)	#kmer of interest

table2 <- read.table(args[2], header=TRUE)	#kmer of interes

#cut only kmer id fsr
table1 <- table1[,c(1:5)]
table2 <- table2[,c(1:5)]

#make new compare table 
tab <- cbind(table1[,c(2,1,3,4,5)], table2[,c(1,3,4,5)])

#add difference in fsr columns
tab[,(dim(tab)[2]+1)] <- tab[,3] - tab[,7]
tab[,(dim(tab)[2]+1)] <- tab[,4] - tab[,8]
tab[,(dim(tab)[2]+1)] <- tab[,5] - tab[,9]

#add flag that indicates if kmers are the same or suffer a variation (only snps supported more or less)
tab[,(dim(tab)[2]+1)] <- as.character(tab[,2]) != as.character(tab[,5])
#set new column names
names(tab) <- c("id", "kmer1", "SFR_merged1", "SFR_plus1", "SFR_minus1", "kmer2", "SFR_merged2", "SFR_plus2", "SFR_minus2", "DMG_merged", "DMG_plus", "DMG_minus", "var_flag")

#order table by difference 
tab <- tab[order(tab[,10],decreasing=TRUE),]

#return print
write.table(tab, file="", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


