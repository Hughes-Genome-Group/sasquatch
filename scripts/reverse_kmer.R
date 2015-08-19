args <- commandArgs(trailingOnly = TRUE)

kmer <- args[1]

#tiny function to create the reverse complement (inverted) of a kmer
mersplit <- unlist(strsplit(kmer,''))
#create reverse complement
rc <- paste0(rev(unlist(lapply(mersplit,function(x){
	if(x == 'A'){ x='T'; }
	else if(x == 'G'){ x='C'; }
	else if(x == 'C'){ x='G'; }
	else if(x == 'T'){ x='A'; }
	else if(x == 'W'){ x='W'; }
	else if(x == 'S'){ x='S'; }
	else if(x == 'N'){ x='N'; }
	else if(x == 'K'){ x='M'; }
	else if(x == 'M'){ x='K'; }
	else if(x == 'R'){ x='Y'; }
	else if(x == 'Y'){ x='R'; }
	else if(x == 'T'){ x='A'; }
	else if(x == 'D'){ x='H'; }
	else if(x == 'H'){ x='D'; }
	else if(x == 'B'){ x='V'; }
	else if(x == 'V'){ x='B'; }

}))),collapse="")

cat(rc)



