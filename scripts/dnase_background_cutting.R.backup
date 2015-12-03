#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[6]) #source common functions from file

kmer <- args[1]	#kmer of interest

profile <- args[2]	#dnase background cut frequency input file
profile <- strsplit(profile, ":")[[1]]	#split the combined string at :s
profile <- as.numeric(profile)

fraction <- as.numeric(args[3])	#fraction to calculate iqr (4 = 25% IQR)
scale.thresh <- as.numeric(args[4])		#MSQ threshold to decide if to normalize, if to extend teh norm window or if to flag as not normalized, default=0.09 but its empirically adjusted to 1-2 datasets so should somehow be adjustable later on
 
max.extension <- as.numeric(args[5])	#maximal extension allowed for extending the normalization window, default = 5 


#START#
msq.dnase.list <- msq.dnase(kmer, profile, fraction=fraction, scale.thresh=scale.thresh, max.ext=max.extension)	

#return msq and iqr
cat(paste0("dnase_msq=", msq.dnase.list$err.sum, " dnase_iqr=", msq.dnase.list$err.iqr, " skip_flag=", msq.dnase.list$skip.flag, " extension=", msq.dnase.list$extension))

