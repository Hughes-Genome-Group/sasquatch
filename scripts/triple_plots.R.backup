suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[14])	#source common functions

kmer <- args[1]	#kmer of interest

profile.plus <- as.numeric( strsplit(args[2], ":")[[1]] )	#data plus strand profile #split the combined string at :s

profile.minus <- as.numeric( strsplit(args[3], ":")[[1]] )	#data minus strand profile #split the combined string at :s

profile.naked.plus <- as.numeric( strsplit(args[4], ":")[[1]] )	#background plus strand profile #split the combined string at :s

profile.naked.minus <- as.numeric( strsplit(args[5], ":")[[1]] )	#background minus strand profile #split the combined string at :s

#occurence counts
count.plus <- args[6]
count.minus <- args[7]
count.naked.plus <- args[8]
count.naked.minus <- args[9]

smooth.flag <- as.numeric(args[10])	#to smooth or not to smooth

skip.flag <- args[11]

extension <- as.numeric(args[12])

plot.dir <- args[13]

# START #
kl=nchar(kmer)

#make plots in lists
profile.raw.list <- plot.profile.raw.data(kmer, profile.plus=profile.plus, profile.minus=profile.minus, count.plus=count.plus, count.minus=count.minus, smooth.flag=smooth.flag)

profile.naked.list <- plot.profile.background.cutting(kmer, profile.plus=profile.naked.plus, profile.minus=profile.naked.minus, count.plus=count.naked.plus, count.minus=count.naked.minus, smooth.flag=smooth.flag)

profile.norm.list <- plot.profile.normalized(kmer, profile.plus, profile.minus, profile.naked.plus, profile.naked.minus, skip.flag, extension)

#save as triple plots

#plot is stored in theplot object #my current solution is to plot it to the desired dir if Flag TRUE /not sure hwo you want to handle plots
#plus
png.file=paste0(plot.dir,'/triple_plot_',kmer,'_plus.png')
png(png.file, width=900, height=900)
		
		grid.newpage() # Open a new page on grid device
		pushViewport(viewport(layout = grid.layout(3, 1)))
		print( profile.raw.list$pp, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))	#real
		print( profile.naked.list$pp, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) #naked
		print( profile.norm.list$pp, vp = viewport(layout.pos.row = 3, layout.pos.col = 1)) #window normalized

suppress.output <- dev.off()

#minus
png.file=paste0(plot.dir,'/triple_plot_',kmer,'_minus.png')
png(png.file, width=900, height=900)

		grid.newpage() # Open a new page on grid device
		pushViewport(viewport(layout = grid.layout(3, 1)))
		print( profile.raw.list$pm, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))	#real
		print( profile.naked.list$pm, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) #naked
		print( profile.norm.list$pm, vp = viewport(layout.pos.row = 3, layout.pos.col = 1)) #window normalized

suppress.output <- dev.off()



