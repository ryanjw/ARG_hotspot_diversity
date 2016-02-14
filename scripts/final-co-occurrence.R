#combining all args together
meta<-read.csv("./data/mgrast-allwgs.csv")
m_reca<-merge(meta,reca_summary, by.x="id",by.y="variable")
m_arg<-merge(meta,args_summary, by.x="id",by.y="dup")
m_int<-merge(meta,int_summary_slim,by.x="id",by.y="variable")
env<-read.csv("./data/project_environment_key.csv")

m_arg<-merge(subset(env,environment=="soil"),m_arg,by="project")


# subbing out indicators
m_arg_sub<-subset(m_arg, ARG == "AAK64454" |
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

sample_count<-ddply(m_arg_sub, .(project), summarise,samples=length(unique(name)))

m_arg_sub<-merge(m_arg_sub, sample_count,by="project")
m_arg_repped<-subset(m_arg_sub, samples >=5)
dim(m_arg_repped)
head(m_arg_repped)

soils_hotargs_cast<-dcast(m_arg_repped, project+name~ARG, value.var="count_arg",fun.aggregate=mean, fill=0 )
dim(soils_hotargs_cast)

m_int_slim<-m_int[ m_int$project %in% m_arg_repped$project,]
head(m_int_slim)
soils_int_cast<-dcast(subset(m_int_slim, count > 0), project+name~IntI1, value.var="count",fun.aggregate=mean, fill=0 )
head(soils_int_cast)

args_ints<-merge(soils_hotargs_cast,soils_int_cast,by=c("project","name"))
library(Hmisc)
library(fdrtool)
library(vegan)
library(plyr)
library(reshape2)

head(m_reca)
m_reca_slim<-m_reca[ m_reca$project %in% m_arg_repped$project,]
length(unique(m_reca_slim$project))

# for loop for co-occurrence

treatments<-as.vector(unique(args_ints$project))
treatments
final_full_results<-data.frame()

for(i in 4:length(treatments)){
	#subset the data for a particular treatment
	# i<-1
	temp<-subset(args_ints, project==treatments[i])
	recas_temp<-dcast(subset(m_reca_slim, project==treatments[i]), project+name~reca, value.var="count",fun.aggregate=mean, fill=0 )
	temp_all<-merge(temp,recas_temp,by=c("project","name"))
	# making an object that has all the results in it (both rho and P values)
	results<-rcorr(as.matrix(temp_all[,-c(1:2)]),type="spearman")
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
	merged<-merge(ps_sub,rhos_melt,by=c("X1","X2"))
	# removing negative rhos
	merged_sub<-subset(merged)
	merged_sub$trt<-treatments[i]
	final_full_results<-rbind(final_full_results, merged_sub)
	print(paste("finished ",treatments[i],sep=""))
}

final_full_results<-na.omit(final_full_results)

stats_final<-ddply(final_full_results, .(Var1,Var2),summarise,.progress="text",
rhos=mean(rho),
rho_lo=quantile(rho,0.025),
qvals=mean(qval),
qvals_hi=quantile(qval, 0.975,na.rm=T)
)

g_list<-(subset(stats_final, rho_lo > 0 & qvals_hi < 0.05))


# making lists of variable types
recas<-as.vector(unique(m_reca$reca))
ints<-as.vector(unique(m_int$IntI1))
args<-as.vector(unique(m_arg_repped$ARG))
recas<-data.frame(recas)
ints<-data.frame(ints)
args<-data.frame(args)
recas$var_type<-"reca"
ints$var_type<-"int"
args$var_type<-"arg"
names(recas)<-names(ints)<-names(args)
var_types<-rbind(recas,ints,args)
head(var_types)
g_list<-merge(g_list,var_types,by.x="Var1",by.y="args")
g_list<-merge(g_list,var_types,by.x="Var2",by.y="args")

head(g_list)
 # all recas were actually miss labelled ints so will remove manually
g_list_sub<-subset(g_list, var_type.x != "reca" & var_type.y !="reca")
unique(g_list_sub$var_type.y)


# making the network
gnet<-igraph::graph.edgelist(as.matrix(g_list_sub[,c(1,2)]),directed=F)
E(gnet)$weight<-g_list_sub$rhos
gnet<-simplify(gnet)
data.frame(igraph::degree(gnet))
vs<-df$vertexes
vs
vs_m<-merge(vs, var_types, by.x="vertex.names",by.y="args")
dim(vs_m)
head(vs_m)
vs_m<-subset(vs_m, var_type != "reca")

# for manual curation
write.table(vs_m, "./data/vertex_list_for_network.txt",sep="\t",row.names=FALSE)
vs_m<-read.table("./data/vertex_list_for_network.txt",sep="\t",header=TRUE)
vs_m<-arrange(vs_m, intergraph_id)
head(vs_m)
as.vector(unique(vs_m$gene_name))
gnet<-asNetwork(gnet)
df<-asDF(gnet)

ggnet(gnet, size=0,method="kamadakawaii",segment.size=df$edges$weight/2)+geom_point(aes(colour=vs_m$gene_name,shape=vs_m$antibiotic),size=4,alpha=0.75)


+geom_point(aes(colour=vs$metadata,shape=vs$shapes),size=4,alpha=0.75)+scale_colour_manual(name="Gene Info",values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","black","grey40"))+scale_shape_manual(name="Antibiotic",values=c(17,5,3,7,8,15,11,16))
