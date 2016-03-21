#calculate heavy smooth correlation and create plot if required
#get Arguments
args <- commandArgs(trailingOnly = TRUE)

kmer=args[1]	#kmer of interest
bandwidth=args[2]	#bandwidth for heavy gaussian smoothing # default = 10

profile=args[3]		#retrieved and stored profile
profile.naked=args[4]	#naked background profile

toplot.flag=args[5]		#TRUE/FALSE if to process and retrieve the plots or not 

if(toplot.flag == 1){	#require only if flag is TRUE
	plot.dir=args[6]	#directory to write the plot to 
}

#split up into profiles from unix stored string
profile <- strsplit(profile, ":")[[1]]
profile.naked <- strsplit(profile.naked, ":")[[1]]

pp <- ksmooth(c(1:length(profile)), profile, kernel="normal", bandwidth=bandwidth)$y	#smooth
pn <- ksmooth(c(1:length(profile.naked)), profile.naked, kernel="normal", bandwidth=bandwidth)$y	#smooth

pp <- pp / sum(pp)
pn <- pn / sum(pn)

hscorrelation <- cor(pp, pn)	#calculate correlation

cat(hscorrelation)	# return heavy smooth correlation value

if(toplot.flag  == 1){	#only if toplotflag = TRUE

	suppressMessages(library(ggplot2))
	#assemble in dataframe for plotting 
	df <- data.frame(x=c(1:length(pp)), value=c(pp,pn), species=c(rep("data profile", length(pp)),rep("background profile", length(pn))))

	theplot <- ggplot(df, aes(x=x, y=value , group=species, colour=species),  xlab="bp idx", ylab="heavy smoothed cut probability", main=paste0("Heavy Smoothed Profiles ",kmer)) + geom_line( ) + ylim(0,0.01)

	#plot is stored in theplot object #my current solution is to plot it to the desired dir if Flag TRUE /not sure hwo you want to handle plots
	png.file=paste0(plot.dir,'/heavy_smooth_correlation_plot_',kmer,'.png')
	png(png.file, width=900, height=300)
		print(theplot)
	suppress.output <- dev.off()
}	



