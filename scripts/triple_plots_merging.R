suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[8])	#source common functions

kmer <- args[1]	#kmer of interest

profile.merged <- as.numeric( strsplit(args[2], ":")[[1]] )

profile.plus <- as.numeric( strsplit(args[3], ":")[[1]] )

profile.minus <- as.numeric( strsplit(args[4], ":")[[1]] )

#occurence counts
count.plus <- args[5]

smooth.flag <- as.numeric(args[6])	#to smooth or not to smooth
plot.dir <- args[7]

# START #
kl=nchar(kmer)

profile.merged.list <- plot.profile.pnorm.data.single(kmer, profile.plus=profile.merged, count.plus=count.plus, smooth.flag=smooth.flag)
profile.plus.list <- plot.profile.pnorm.data.single(kmer, profile.plus=profile.plus, count.plus=count.plus, smooth.flag=smooth.flag)
profile.minus.list <- plot.profile.pnorm.data.single(kmer, profile.plus=profile.minus, count.plus=count.plus, smooth.flag=smooth.flag)

#save as triple plots
png.file=paste0(plot.dir,'/triple_plot_',kmer,'_merging.png')
png(png.file, width=900, height=900)
		
		grid.newpage() # Open a new page on grid device
		pushViewport(viewport(layout = grid.layout(3, 1)))
		print( profile.merged.list$pp, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))	#real
		print( profile.plus.list$pp, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) #naked
		print( profile.minus.list$pp, vp = viewport(layout.pos.row = 3, layout.pos.col = 1)) #window normalized

suppress.output <- dev.off()

