#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[6]) #source common functions from file

kmer <- args[1]	#kmer of interest

infile.base <- args[2]	#base frequency input file

toplot.flag <- args[3]	#0 if no potting 1 if to plot

if(toplot.flag == 1){
	plot.dir <- args[4]	#plot dir onl required if lot flag true
	bandwidth <- as.numeric(args[5])	#bandwidth to smooth base composition plot default=3
}


#START#

#get base composition ratio for each base
get.base.list <- get.comp(kmer, infile.base)

#calculate median sqqare error  and iqr for bases
msq.base.list <- msq.base(kmer, get.base.list$a, get.base.list$c, get.base.list$g, get.base.list$t)	#get msq msqprofile and iqr

#return msq and iqr
cat(paste0("base_msq=", msq.base.list$err.sum, " base_iqr=", msq.base.list$err.iqr))

#if plot flag ture produce a plot out of that
if(toplot.flag == 1){

	suppressMessages(library(ggplot2))
			
	base.comp.plot <- print.base.comp(kmer=kmer, bandwidth=bandwidth, a=get.base.list$a, c=get.base.list$c, g=get.base.list$g, t=get.base.list$t) #maake the plot

	#plot is stored in plot object #my current solution is  again to plot it to the desired dir if Flag TRUE /not sure hwo you want to handle plots
	png.file=paste0(plot.dir,'/base_composition_plot_',kmer,'.png')
	png(png.file, width=900, height=300)
		print(base.comp.plot)
	suppress.output <- dev.off()

}


