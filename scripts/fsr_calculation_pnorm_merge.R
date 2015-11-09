#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[5]) #source common functions from file

kmer <- args[1]	#kmer of interest

profile.merged <- as.numeric( strsplit(args[2], ":")[[1]] )	#data strand merged profile #split the combined string at :s

profile.plus <- as.numeric( strsplit(args[3], ":")[[1]] )	#data plus strand profile #split the combined string at :s

profile.minus <- as.numeric( strsplit(args[4], ":")[[1]] )	#data minus strand profile #split the combined string at :s

# START #
kl=nchar(kmer)

#gaussian smooth normalied profiles
profile.merged <- ksmooth(c(1:length(profile.merged)), profile.merged, kernel="normal", bandwidth=5)$y
profile.merged <- profile.merged/sum(profile.merged)

profile.plus <- ksmooth(c(1:length(profile.plus)), profile.plus, kernel="normal", bandwidth=5)$y
profile.plus <- profile.plus/sum(profile.plus)

profile.minus <- ksmooth(c(1:length(profile.minus)), profile.minus, kernel="normal", bandwidth=5)$y
profile.minus <- profile.minus/sum(profile.minus)

#estimate footprint, borders and best shoulder ranges
border.list.plus <- sobel.borders(profile=profile.plus, kl=kl)
border.list.minus <- sobel.borders(profile=profile.minus, kl=kl)
border.list.merged <- sobel.borders(profile=profile.merged, kl=kl)

#calculate adjusted SFR
if(border.list.merged$flag){
	fsr.merged <- adj.fsr(profile=profile.merged, kl=kl, usb=border.list.merged$us, dsb=border.list.merged$ds, usr=border.list.merged$range.us, dsr=border.list.merged$range.ds)
}else{	
	fsr.merged="1.0" 
}

if(border.list.plus$flag){
	fsr.plus <- adj.fsr(profile=profile.plus, kl=kl, usb=border.list.plus$us, dsb=border.list.plus$ds, usr=border.list.plus$range.us, dsr=border.list.plus$range.ds)
}else{	
	fsr.plus="1.0" 
}

if(border.list.minus$flag){
	fsr.minus <- adj.fsr(profile=profile.minus, kl=kl, usb=border.list.minus$us, dsb=border.list.minus$ds, usr=border.list.minus$range.us, dsr=border.list.plus$range.ds)
}else{	
	fsr.minus="1.0" 
}

#return / print
cat(paste0("fsr_merged=", fsr.merged," fsr_plus=", fsr.plus, " fsr_minus=", fsr.minus))





