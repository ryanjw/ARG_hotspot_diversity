# This script is to join all the different datasets

#Setting working directory
setwd("~/Google Drive/postdoc_writing/ARG/ARG_hotspot_diversity")

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
data<-read.table("./data/raw_data/arg-summary-count.tsv",header=T,sep="\t",check.names=F)
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
fulldim<-dim(args_super_bind_rm)[1]
args_super_bind_rm$splitter<-NA
for(i in 1:fulldim){
	
	args_super_bind_rm$splitter[i]<-rev(unlist(strsplit(as.vector(args_super_bind_rm$dup[i]),fixed=T,split=".")))[1]	
	counter<-counter+1
	if(counter > 50000){
		print(i/fulldim)
		counter<-0
	}
}

# splitting out groups that have the 1 now after the function above and turning them into a factor
args_subset_1<-subset(args_super_bind_rm, splitter==1)
args_subset_1$variable<-factor(args_subset_1$variable)

# now taking all of the ones that do not have that and changing the id's to have the proper .3 at the end
args_subset_0<-subset(args_super_bind_rm, splitter != 1)
args_subset_0$dup<-factor(paste(args_subset_0$dup,".3",sep=""))

#now putting everything together with the proper decimals
args_proper_dec<-rbind(args_subset_1,args_subset_0)

# can write this point to avoid doing this again
# write.table(args_proper_dec, "./args_longform_possible_duplicates.txt",row.names=F,sep="\t")
# write.table(args_proper_dec, "~/Desktop/temp/args_longform_possible_duplicates.txt",row.names=F,sep="\t")
# now we are going to make sure we have a count per ARG per sample under "dup", then we will make a data.frame with fun.aggregate=mean.  If there are duplicate samples, we will therefore get the average between the two (they would be the same number)
args_summary<-ddply(args_proper_dec, .(dup,ARG),summarise,.progress="text",
count=mean(value)
)
args_cast<-dcast(args_summary, dup~ARG, value.var="count",fun.aggregate=mean,fill=0)



# Now doing this with intI1

int1<-read.table("./data/raw_data/int-summary.tsv",header=T,sep="\t",check.names=F)
int2<-read.table("./data/raw_data/int-summary2.tsv",header=T,sep="\t",check.names=F)

int1_corr<-row_column_correction(int1,"IntI1")
int2_corr<-row_column_correction(int2,"IntI1")

int1_melt<-melt(int1_corr, id="IntI1")
int2_melt<-melt(int2_corr, id="IntI1")

int1_orig<-read.table("./data/raw_data/int1-summary-count.tsv",header=T,sep="\t",check.names=F,strip.white=TRUE,row.names=1)
int1_orig_corr<-row_column_correction(int1_orig, "IntI1")
int1_orig_melt<-melt(int1_orig_corr,id="IntI1" )

intI1_super_bind<-rbind(int1_orig_melt,int1_melt,int2_melt)
intI1_super_bind$splitter<-NA
intI1_super_bind_rm<-subset(intI1_super_bind)

counter<-0
for(i in 1:dim(intI1_super_bind_rm)[1]){
	
	thing<-rev(unlist(strsplit(as.vector(intI1_super_bind_rm$variable[i]),fixed=T,split=".")))[1]	
	intI1_super_bind_rm$splitter[i]<-thing
	counter<-counter+1
	if(counter > 6000){
		print(i/dim(intI1_super_bind_rm)[1])
	}
}

int_subset_1<-subset(intI1_super_bind_rm, splitter=="1")
int_subset_0<-subset(intI1_super_bind_rm, splitter != "1")
int_subset_0$variable<-factor(paste(int_subset_0$variable,".3",sep=""))
int_proper_dec<-rbind(int_subset_1,int_subset_0)
int_summary<-ddply(int_proper_dec, .(variable,IntI1),summarise,.progress="text",
count=mean(value)
)

funky_ints<-read.csv("./data/funky-ints.csv",header=F)
head(int_summary)
int_summary_slim<-int_summary[! int_summary$IntI1 %in% funky_ints$V1,]
dim(int_summary_slim)
int_cast<-dcast(int_summary, variable~IntI1, value.var="count",fun.aggregate=mean,fill=0)



# now doing with recas
rec1<-read.table("./data/raw_data/rec-summary.tsv",header=T,sep="\t",check.names=F)
rec2<-read.table("./data/raw_data/rec-summary2.tsv",header=T,sep="\t",check.names=F)

rec1_corr<-row_column_correction(rec1,"reca")
rec2_corr<-row_column_correction(rec2,"reca")

rec1_corr<-row_column_correction(rec1,"reca")
rec2_corr<-row_column_correction(rec2,"reca")

rec1_melt<-melt(rec1_corr, id="reca")
rec2_melt<-melt(rec2_corr, id="reca")

head(rec1_melt)
rec1_melt<-subset(rec1_melt, value > 0)
rec1_melt$V1<-NA
for(i in 1:dim(rec1_melt)[1]){
	p1<-unlist(strsplit(as.vector(rec1_melt$reca[i]),fixed=TRUE,split="|"))
	rec1_melt$V1[i]<-paste("sp",p1[2],p1[3],sep="|")
}

rec1_slim<-rec1_melt[! rec1_melt$V1 %in% nono_list$V1,]
dim(rec1_melt)
head(rec1_melt)
length(unique(rec1_melt$reca))
length(unique(nono_list$V1))
length(unique(rec1_slim$reca))


rec2_melt<-subset(rec2_melt, value > 0)
rec2_melt$V1<-NA
for(i in 1:dim(rec2_melt)[1]){
	p1<-unlist(strsplit(as.vector(rec2_melt$reca[i]),fixed=TRUE,split="|"))
	rec2_melt$V1[i]<-paste("sp",p1[2],p1[3],sep="|")
}

rec2_slim<-rec2_melt[! rec2_melt$V1 %in% nono_list$V1,]

# need to get orig orig set


reca<-read.table("./data/raw_data/summary-count.tsv",header=T,sep="\t",check.names=F)
reca_orig_corr<-row_column_correction(reca,"reca")
reca_orig_melt<-melt(reca_orig_corr,id="reca" )
reca_orig_melt<-subset(reca_orig_melt, value > 0)
reca_orig_melt$V1<-NA
for(i in 1:dim(reca_orig_melt)[1]){
	p1<-unlist(strsplit(as.vector(reca_orig_melt$reca[i]),fixed=TRUE,split="|"))
	reca_orig_melt$V1[i]<-paste("sp",p1[2],p1[3],sep="|")
}

reca_orig_slim<-reca_orig_melt[! reca_orig_melt$V1 %in% nono_list$V1,]
reca_super_bind<-rbind(reca_orig_slim,rec1_slim,rec2_slim)

reca_super_bind$splitter<-NA
reca_super_bind_rm<-subset(reca_super_bind, value > 0)

nono_list<-read.table("./data/funky-recas.txt",sep="\t",strip.white=TRUE)
as.matrix(head(nono_list))

# counter<-0
# fulldim<-dim(reca_super_bind_rm)[1]
# for(i in 3420:fulldim){
	
	# reca_super_bind_rm$splitter[i]<-rev(unlist(strsplit(as.character(reca_super_bind_rm$variable[i]),fixed=T,split=".")))[1]	
	# counter<-counter+1
	# if(counter > 1000){
		# print(i/fulldim)
		# counter<-0
	# }
# }
#this is a faster way
df<-data.frame(as.vector(unique(reca_super_bind_rm$variable)))
names(df)[1]<-"variable"
df$splitter<-(NA)
str(df)
counter<-0
fulldim<-dim(df)[1]
for(i in 1:fulldim){
	
	df$splitter[i]<-rev(unlist(strsplit(as.character(df$variable[i]),fixed=T,split=".")))[1]	
	counter<-counter+1
	if(counter > 10){
		print(i/fulldim)
		counter<-0
	}
}

reca_super_bind_rm_m<-merge(reca_super_bind_rm[,-5],df,by="variable")


reca_subset_1<-subset(reca_super_bind_rm_m, splitter=="1")
reca_subset_0<-subset(reca_super_bind_rm_m, splitter != "1")
length(unique(reca_subset_1$variable))
length(unique(reca_subset_0$variable))
head(reca_subset_0)
reca_subset_0$variable<-factor(paste(reca_subset_0$variable,".3",sep=""))

reca_proper_dec<-rbind(reca_subset_1,reca_subset_0)
head(reca_proper_dec)
length(unique(reca_proper_dec$variable))
head(int_proper_dec)

reca_summary<-ddply(reca_proper_dec, .(variable,reca),summarise,.progress="text",
count=mean(value)
)
head(reca_summary)
dim(reca_summary)
length(unique(int_summary$variable))
reca_cast<-dcast(reca_summary, variable~reca, value.var="count",fun.aggregate=mean,fill=0)


