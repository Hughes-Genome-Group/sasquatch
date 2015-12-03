suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))

#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[25])	#source common functions

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
count.naked.plus1 <- args[13]
count.naked.minus1 <- args[14]
count.plus2 <- args[15]
count.minus2 <- args[16]
count.naked.plus2 <- args[17]
count.naked.minus2 <- args[18]

smooth.flag <- as.numeric(args[19])	#to smooth or not to smooth

skip.flag1 <- args[20]
skip.flag2 <- args[21]

extension1 <- as.numeric(args[22])
extension2 <- as.numeric(args[23])
plot.dir <- args[24]


# START #
kl=nchar(kmer1)

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



