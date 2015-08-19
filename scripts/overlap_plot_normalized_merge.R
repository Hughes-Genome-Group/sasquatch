suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[21])	#source common functions

kmer1 <- args[1]	#kmer of interest
kmer2 <- args[2]

profile.plus1 <- as.numeric( strsplit(args[3], ":")[[1]] )	#data plus strand profile #split the combined string at :s
profile.minus1 <- as.numeric( strsplit(args[4], ":")[[1]] )	#data minus strand profile #split the combined string at :s
profile.naked.plus1 <- as.numeric( strsplit(args[5], ":")[[1]] )	#background plus strand profile #split the combined string at :s
profile.naked.minus1 <- as.numeric( strsplit(args[6], ":")[[1]] )	#background minus strand profile #split the combined string at :s
profile.plus2 <- as.numeric( strsplit(args[7], ":")[[1]] )	#data plus strand profile #split the combined string at :s
profile.minus2 <- as.numeric( strsplit(args[8], ":")[[1]] )	#data minus strand profile #split the combined string at :s
profile.naked.plus2 <- as.numeric( strsplit(args[9], ":")[[1]] )	#background plus strand profile #split the combined string at :s
profile.naked.minus2 <- as.numeric( strsplit(args[10], ":")[[1]] )	#background minus strand profile #split the combined string at :s

#occurence counts
count.plus1 <- args[11]
count.minus1 <- args[12]
count.plus2 <- args[13]
count.minus2 <- args[14]

smooth.flag <- as.numeric(args[15])	#to smooth or not to smooth

skip.flag1 <- args[16]
skip.flag2 <- args[17]

extension1 <- as.numeric(args[18])
extension2 <- as.numeric(args[19])
plot.dir <- args[20]


# START #
kl=nchar(kmer1)

#make plots in lists
profile.over.list <- plot.profile.single.norm.overlap.merge(kmer1, kmer2, profile.plus1, profile.plus2, profile.minus1, profile.minus2, profile.naked.plus1, profile.naked.plus2, profile.naked.minus1, profile.naked.minus2, skip.flag1, skip.flag2, extension1, extension2, count.plus1, count.plus2, count.minus1, count.minus2)

#save as single overlap plot

#plot is stored in theplot object #my current solution is to plot it to the desired dir if Flag TRUE /not sure hwo you want to handle plots
#plus
png.file=paste0(plot.dir,'/compare_overlap_plot_',kmer1,'_',kmer2,'_plus.png')
png(png.file, width=900, height=300)
		
		print( profile.over.list$pp ) #window normalized

suppress.output <- dev.off()


