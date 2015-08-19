suppressMessages(library(ggplot2))

#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[8])	#get common functions from file

kmer <- args[1]	#kmer of interest

profile.plus <- as.numeric( strsplit(args[2], ":")[[1]] )	#data plus strand profile #split the combined string at :s

profile.minus <- as.numeric( strsplit(args[3], ":")[[1]] )	#data minus strand profile #split the combined string at :s

#occurence counts
count.plus <- args[4]
count.minus <- args[5]

smooth.flag <- args[6]	#to smooth or not to smooth

plot.dir <- args[7]

# START #
kl=nchar(kmer)

#make plots in lists
profile.norm.list <- plot.profile.pnorm.data(kmer, profile.plus=profile.plus, profile.minus=profile.minus, count.plus=count.plus, count.minus=count.minus, smooth.flag=smooth.flag)

#plot is stored in theplot object #my current solution is to plot it to the desired dir if Flag TRUE /not sure hwo you want to handle plots
#plus
png.file=paste0(plot.dir,'/norm_data_plot_',kmer,'_plus.png')
png(png.file, width=900, height=300)

		print( profile.norm.list$pp )	#real

suppress.output <- dev.off()

#minus
png.file=paste0(plot.dir,'/norm_data_plot_',kmer,'_minus.png')
png(png.file, width=900, height=300)


		print( profile.norm.list$pm )	#real

suppress.output <- dev.off()



