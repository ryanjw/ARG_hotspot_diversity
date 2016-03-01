# this script contains other iterations of co-occurrence analysis.  While not contained in the paper, for transparency of the analytical process, they were shared

meta<-read.csv("~/Desktop/big_ARG_data/mgrast-allwgs.csv")
head(meta)
length(unique(meta$id))
meta$id<-as.factor(meta$id)
rec1<-read.table("~/Desktop/big_ARG_data/rec-summary.tsv",header=T,sep="\t",check.names=F)
rec2<-read.table("~/Desktop/big_ARG_data/rec-summary2.tsv",header=T,sep="\t",check.names=F)

rec1_corr<-row_column_correction(rec1,"reca")
rec2_corr<-row_column_correction(rec2,"reca")

rec1_corr<-row_column_correction(rec1,"reca")
rec2_corr<-row_column_correction(rec2,"reca")

rec1_melt<-melt(rec1_corr, id="reca")
rec2_melt<-melt(rec2_corr, id="reca")
length(unique(rec1_melt$variable))
length(unique(rec2_melt$variable))
head(rec1_melt)
# need to get orig orig set
reca<-read.table("~/Desktop/big_ARG_data/reca-summary-count.tsv",header=T,sep="\t",check.names=F)
dim(reca)
reca_orig_corr<-row_column_correction(reca,"reca")
head(reca_orig_corr)
reca_orig_melt<-melt(reca_orig_corr,id="reca" )
head(reca_orig_melt)
reca_super_bind<-rbind(reca_orig_melt,rec1_melt,rec2_melt)
head(reca_super_bind)
reca_cast<-dcast(reca_super_bind, variable~reca, value.var="value", fun.aggregate=mean, fill=0)
reca_meta<-merge(reca_cast,meta[,c(2:10)], by.x="variable",by.y="id_nopt")
dim(reca_meta)
reca_meta<-reca_meta[,c(1,70687:70694,2:70686)]
# write.table(reca_meta,"~/Desktop/recA_4_cocur.txt",sep="\t",row.names=F)
head(reca_meta[,1:10])
int1<-read.table("~/Desktop/big_ARG_data/int-summary.tsv",header=T,sep="\t",check.names=F)
int2<-read.table("~/Desktop/big_ARG_data/int-summary2.tsv",header=T,sep="\t",check.names=F)

int1_corr<-row_column_correction(int1,"int")
int2_corr<-row_column_correction(int2,"int")

int1_corr<-row_column_correction(int1,"int")
int2_corr<-row_column_correction(int2,"int")

int1_melt<-melt(int1_corr, id="int")
int2_melt<-melt(int2_corr, id="int")
length(unique(int1_melt$variable))
length(unique(int2_melt$variable))
head(int1_melt)

setwd("~/Google Drive/postdoc_writing/ARG/ARG_diversity/ARG_biodiversity/")
int_orig<-read.table("./datasets/int1-summary-count.tsv",header=T,sep="\t",check.names=F,strip.white=TRUE,row.names=1)
int_orig_corr<-row_column_correction(int_orig, "int")
int_orig_melt<-melt(int_orig_corr,id="int" )

int_super_bind<-rbind(int_orig_melt,int1_melt,int2_melt)
head(int_super_bind)
int_cast<-dcast(int_super_bind, variable~int, value.var="value", fun.aggregate=mean, fill=0)
int_meta<-merge(int_cast,meta[,c(2:10)], by.x="variable",by.y="id_nopt")
dim(int_meta)
head(int_meta[,1:10])
int_meta<-int_meta[,c(1,213:219,2:212)]
head(int_meta)
int_meta<-int_meta[,-219]
all_w_args<-read.csv("~/Desktop/big_ARG_data/all_ARG_data_with_covariates_for_pcoa.csv",check.names=F,header=T)
head(all_w_args[,1:20])
as.vector(unique(all_w_args$environment))
soils<-subset(all_w_args, environment=="soil")
sample_count<-ddply(soils, .(project),summarise,
samps=length(unique(variable))
)
head(soils[,1:20])
soils_melt<-melt(soils[,-c(1,2,4,8:15)],id=c("variable","project","name","bps"))
head(soils_melt)
names(soils_melt)[5]<-"ARG"
soils_melt_arg_sub<-subset(soils_melt, ARG == "AAK64454" |
ARG == "AAP46617" |
ARG == "AAQ12900" |
ARG == "AAV74563" |
ARG == "AAY62597" |
ARG == "ABN79417" |
ARG == "ACA23192" |
ARG == "YP_001691306" |
ARG == "YP_001693856" |
ARG == "YP_002242988" |
ARG == "YP_002850805" |
ARG == "YP_407313" |
ARG == "YP_788911" |
ARG == "ZP_03029420" |
ARG == "ZP_03063464" |
ARG == "ZP_03452771" |
ARG == "ZP_04466634" |
ARG == "ZP_04557016" |
ARG == "ZP_04576183" |
ARG == "ZP_04623897" )

soils_melt_merge_sampcount<-merge(soils_melt_arg_sub, sample_count,by="project")
head(soils_melt_merge_sampcount)

soils_hotargs_subbed<-subset(soils_melt_merge_sampcount, samps >=5)

soils_hotargs_cast<-dcast(soils_hotargs_subbed, project+variable~ARG, value.var="value",fun.aggregate=mean, fill=0 )

dim(soils_hotargs_cast)
recas<-as.matrix(reca_meta[,-c(1:9)])
rownames(recas)<-reca_meta$variable
colnames(recas)[1:100]
recas_small<-recas[,colSums(decostand(recas, "pa")) >= 1018/2 ]
dim(recas_small)
recas_small_df<-data.frame(recas_small[-c(63,93),],check.names=F)
recas_small_df<-cbind(rownames(recas_small_df),recas_small_df)
head(recas_small_df[,1:10])
names(recas_small_df)[1]<-"variable"
hotargs_reca<-merge(soils_hotargs_cast, recas_small_df,by="variable")
dim(hotargs_reca)
head(hotargs_reca[,1:10])
names(int_meta[,1:10])
hotargs_reca_int<-merge(hotargs_reca, int_meta[,-c(2:9)],by="variable")
# write.table(hotargs_reca_int, "~/Desktop/hotargs_reca_int.txt",sep="\t",row.names=F)
# hotargs_reca_int<-read.table("~/Desktop/big_ARG_data/hotargs_reca_int.txt",sep="\t",header=T,check.names=F)
hotargs_reca_int<-read.table("./hotargs_reca_int.txt",sep="\t",header=T,check.names=F)
library(Hmisc)
library(fdrtool)
library(vegan)
library(plyr)
library(reshape2)
dataset<-hotargs_reca_int
head(dataset[,1:10])
treatments<-as.vector(unique(dataset$project))
treatments
final_full_results<-data.frame()
for(i in 1:length(treatments)){
	#subset the data for a particular treatment
	temp<-subset(dataset, project==treatments[i])
	# making an object that has all the results in it (both rho and P values)
	results<-rcorr(as.matrix(temp[,-c(1:2)]),type="spearman")
	#make two seperate objects for p-value and correlation coefficients
	rhos<-results$r
	ps<-results$P
	# going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
	ps_melt<-na.omit(melt(ps))
	#creating a qvalue based on FDR
	ps_melt$qval<-fdrtool(ps_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
	#making column names more relevant

	names(ps_melt)[3]<-"pval"
	# if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
	ps_sub<-subset(ps_melt)
	# now melting the rhos, note the similarity between ps_melt and rhos_melt
	rhos_melt<-na.omit(melt(rhos))
	names(rhos_melt)[3]<-"rho"
	#merging together
	merged<-merge(ps_sub,rhos_melt,by=c("Var1","Var2"))
	# removing negative rhos
	merged_sub<-subset(merged)
	merged_sub$trt<-treatments[i]
	final_full_results<-rbind(final_full_results, merged_sub)
	print(paste("finished ",treatments[i],sep=""))
}

head(final_full_results)
final_full_results<-na.omit(final_full_results)
# write.table(final_full_results, "~/Desktop/big_ARG_data/final_network_w_quarter_necessary.txt",sep="\t",row.names=F)
final_full_results<-read.table("~/Desktop/big_ARG_data/final_network_w_quarter_necessary.txt",sep="\t",header=T)
stats_final<-ddply(final_full_results, .(Var1,Var2),summarise,.progress="text",
rhos=mean(rho),
rho_lo=quantile(rho,0.025),
qvals=mean(qval),
qvals_hi=quantile(qval, 0.975,na.rm=T)
)

g_list<-(subset(stats_final, rho_lo > 0 & qvals_hi < 0.05))
dim(g_list)
write.table(g_list, "~/Desktop/network_half.txt",sep="\t",row.names=F)
g_list<-read.table("~/Desktop/network_half.txt",sep="\t",check.names=F,header=T)
dim(g_list)
g_list$var1type<-NA
g_list$var2type<-NA
head(g_list)
counter<-0
for(i in 1:dim(g_list)[1]){
	# i<-1
	thing1<-unlist(strsplit(as.vector(g_list$Var1[i]),split="|",fixed=T))[1]
	thing2<-unlist(strsplit(as.vector(g_list$Var2[i]),split="|",fixed=T))[1]
	thing1
	thing2
	if(thing1 != "tr" & thing1 != "sp"){thing1<-"arg"}
	if(thing2 != "tr" & thing2 != "sp"){thing2<-"arg"}
	g_list$var1type[i]<-thing1
	g_list$var2type[i]<-thing2
	counter<-counter+1
	if(counter > 25){
		print(i/dim(g_list)[1])
	}
}

write.table(g_list, "~/Desktop/g_list.txt",row.names=F,sep="\t")
g_list<-read.table("~/Desktop/big_ARG_data/g_list.txt",header=T,sep="\t")

library(igraph)
head(g_list)
g_list<-subset(g_list, Var1 != "sp|O94762|RECQ5_HUMAN" & Var2 != "sp|O94762|RECQ5_HUMAN" & Var1 != "YP_001691306" & Var2 !="YP_001691306")
gnet<-graph.edgelist(as.matrix(g_list[,c(1,2)]),directed=F)
V(gnet)
E(gnet)$weight<-g_list$rhos
gnet<-simplify(gnet)
plot(gnet)
head(g_list)
subset(g_list, Var1=="tr|A0A0D5ZYN3|A0A0D5ZYN3_ACIBA" | Var2=="tr|A0A0D5ZYN3|A0A0D5ZYN3_ACIBA")
gnet
?degree
df<-data.frame(as.matrix(igraph::degree(gnet)))
df$nodes<-rownames(df)
names(df)[1]<-"degree"
library(plyr)
library(reshape2)
(arrange(df,degree))
library(intergraph)
library(GGally)

gnet<-asNetwork(gnet)
df<-asDF(gnet)
vs<-df$vertexes
vs$vertex_type<-NA
for(i in 1:dim(vs)[1]){
	# i<-1
	thing1<-unlist(strsplit(as.vector(vs$vertex.names[i]),split="|",fixed=T))[1]

	if(thing1 != "tr" & thing1 != "sp"){thing1<-"arg"}
	vs$vertex_type[i]<-thing1

}
vs


write.table(vs,"~/Desktop/big_ARG_data/vertexes.txt",sep="\t",row.names=F)

head(vs)
subset(vs, vertex_type=="sp")
dim(df$edges)
vs$vertex_type<-factor(vs$vertex_type)
levels(vs$vertex_type)<-c("ARG","recA","intI1")
# write.table(vs,"~/Desktop/vertexes.txt",sep="\t",row.names=F)
vs<-read.table("~/Desktop/big_ARG_data/vertexes.txt",sep="\t",header=T)
head(vs)
vs$shapes<-factor(vs$shapes,levels=c("cephalosporin", "fosmidomycin",  "macrolide", "multidrug","penicillin","tetracycline",  "vancomycin","non-ARG"))
vs$metadata<-factor(vs$metadata, levels=c("Actinobacteria","Alphaproteobacteria","Bacilli","Betaproteobacteria","Clostridia","Eukaryote","Gammaproteobacteria","Spirochaetia","ARG","intI1" ))
ggnet(gnet, size=0,method="kamadakawaii",segment.size=df$edges$weight/2)+geom_point(aes(colour=vs$metadata,shape=vs$shapes),size=4,alpha=0.75)+scale_colour_manual(name="Gene Info",values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","black","grey40"))+scale_shape_manual(name="Antibiotic",values=c(17,5,3,7,8,15,11,16))
gnet
g_list
