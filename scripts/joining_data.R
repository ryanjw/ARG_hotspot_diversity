# This script is to join all the different datasets

#Setting working directory
setwd("~/Google Drive/postdoc_writing/ARG/ARG_hotspots/")

#Reading in all necessary libraries
source("./scripts/libraries.R")

#First need to call in function that corrects some issues with column names being offset
source("./scripts/row_column_correction.R")


# pulling in new arg datasets
args<-read.table("./data/raw_data/arg-summary.tsv",header=T,sep="\t",check.names=F)
args2<-read.table("./data/raw_data/arg-summary2.tsv",header=T,sep="\t",check.names=F)

names(args)[-287]<-colnames(args)[-1]
args<-args[,-287]

names(args)[-184]<-colnames(args)[-1]
args<-args[,-184]

args<-cbind(rownames(args),args)
args2<-cbind(rownames(args2),args2)
names(args)[1]<-"ARG"
names(args2)[1]<-"ARG"

args_melt<-melt(args,id="ARG")
args2_melt<-melt(args2,id="ARG")

new_args_bind<-rbind(args_melt, args2_melt)

# pulling in old arg datasets
data<-read.table("./data/raw_data/summary-count.tsv",header=T,sep="\t",check.names=F)
names(data)[-873]<-names(data)[-1]
data<-data[,-873]
data<-cbind(rownames(data),data)
names(data)[1]<-"ARG"
data_melt<-melt(data, id="ARG")

args_super_bind<-rbind(new_args_bind,data_melt)
args_super_bind$dup<-args_super_bind$variable

# some of the unique id's have a .1, this is going to be dealt with

args_super_bind_rm<-subset(args_super_bind, value > 0)
dim(args_super_bind_rm)
counter<-0
for(i in 1:dim(args_super_bind_rm)[1]){
	
	thing<-rev(unlist(strsplit(as.vector(args_super_bind_rm$dup[i]),fixed=T,split=".")))[1]	
	args_super_bind_rm$splitter[i]<-thing
	counter<-counter+1
	if(counter > 6000){
		print(i/dim(args_super_bind_rm)[1])
	}
}
