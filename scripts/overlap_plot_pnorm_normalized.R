suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[13])	#source common functions

kmer1 <- args[1]	#kmer of interest
kmer2 <- args[2]

profile.plus1 <- as.numeric( strsplit(args[3], ":")[[1]] )	#data plus strand profile #split the combined string at :s
profile.minus1 <- as.numeric( strsplit(args[4], ":")[[1]] )	#data minus strand profile #split the combined string at :s

profile.plus2 <- as.numeric( strsplit(args[5], ":")[[1]] )	#data plus strand profile #split the combined string at :s
profile.minus2 <- as.numeric( strsplit(args[6], ":")[[1]] )	#data minus strand profile #split the combined string at :s

#occurence counts
count.plus1 <- args[7]
count.minus1 <- args[8]
count.plus2 <- args[9]
count.minus2 <- args[10]

smooth.flag <- as.numeric(args[11])	#to smooth or not to smooth
plot.dir <- args[12]

# START #
kl=nchar(kmer1)

#smooth if flagged so
if(smooth.flag == 1){
	profile.plus1 <- ksmooth(c(1:length(profile.plus1)), profile.plus1, kernel="normal", bandwidth=5)$y
	profile.plus2 <- ksmooth(c(1:length(profile.plus2)), profile.plus2, kernel="normal", bandwidth=5)$y
	profile.minus1 <- ksmooth(c(1:length(profile.minus1)), profile.minus1, kernel="normal", bandwidth=5)$y
	profile.minus2 <- ksmooth(c(1:length(profile.minus2)), profile.minus2, kernel="normal", bandwidth=5)$y
}

#plotnorm
profile.plus1 <- profile.plus1/sum(profile.plus1)
profile.plus2 <- profile.plus2/sum(profile.plus2)
profile.minus1 <- profile.minus1/sum(profile.minus1)
profile.minus2 <- profile.minus2/sum(profile.minus2)

#make plots in lists
profile.over.list <- plot.profile.single.pnorm.overlap(kmer1, kmer2, profile.plus1, profile.plus2, profile.minus1, profile.minus2, count.plus1, count.plus2, count.minus1, count.minus2)

#save as single overlap plot

#plot is stored in theplot object #my current solution is to plot it to the desired dir if Flag TRUE /not sure hwo you want to handle plots
#plus
png.file=paste0(plot.dir,'/compare_overlap_plot_',kmer1,'_',kmer2,'_plus.png')
png(png.file, width=900, height=300)
		
		print( profile.over.list$pp ) #window normalized

suppress.output <- dev.off()

#minus
png.file=paste0(plot.dir,'/compare_overlap_plot_',kmer1,'_',kmer2,'_minus.png')
png(png.file, width=900, height=300)

		print( profile.over.list$pm ) #window normalized

suppress.output <- dev.off()



