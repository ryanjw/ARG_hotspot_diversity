# variable and thing have to be in quotes
column_splitter<-function(dataset, variable, thing, place, record_split){
dataset$splitter<-NA
counter<-1
for(i in 1:dim(dataset)[1]){
	dataset$splitter[i]<-unlist(strsplit(as.vector(dataset$variable[i]),fixed=T,split=thing))[place]
	counter<-counter+1
	if(counter > record_split){
		print(i/dim(dataset)[1])
		counter<-0
	}	
}
return(dataset)
}
