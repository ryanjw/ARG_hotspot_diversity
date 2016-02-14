#combining args and recas for dissimilarity comparisons
meta<-read.csv("./data/mgrast-allwgs.csv")
m_reca<-merge(meta,reca_summary, by.x="id",by.y="variable")
m_arg<-merge(meta,args_summary, by.x="id",by.y="dup")
names(m_arg)[18]<-"count_arg"
names(m_reca)[18]<-"count_reca"
final_results<-data.frame()
groups<-as.vector(intersect(m_arg$project,m_reca$project))
for(i in 1:length(groups)){
	
	temp_arg<-subset(m_arg, project==groups[i])
	temp_reca<-subset(m_reca, project==groups[i])
	
	arg_cast<-dcast(temp_arg, name~ARG, value.var="count_arg",fun.aggregate=mean,fill=0)
	# dim(arg_cast)
	reca_cast<-dcast(temp_reca, name~reca, value.var="count_reca",fun.aggregate=mean,fill=0)
	# dim(reca_cast)
	
	if(dim(arg_cast)[1] > 1 & dim(reca_cast)[1] > 1){
		m_temp<-merge(arg_cast,reca_cast,by="name")
		arg_mat<-vegdist(decostand(m_temp[,c(2:dim(arg_cast)[2])],"pa"),"jaccard")
		reca_mat<-vegdist(decostand(m_temp[,-c(1:dim(arg_cast)[2])],"pa"),"jaccard")
		rez<-data.frame(groups[i],as.vector(arg_mat),as.vector(reca_mat))
		final_results<-rbind(final_results,rez)
	}
	print(i/length(groups))
}

head(final_results)
names(final_results)<-c("project","arg","reca")
bps<-ddply(meta, .(project),summarise, bps=mean(bps))
final_results_m<-merge(final_results,bps,by="project")
dim(final_results_m)
write.table(final_results_m, "./data/arg_reca_jaccard_distances.txt",sep="\t",row.names=F)

# doing tests

library(lme4)
library(lmerTest)
library(bbmle)
source("./scripts/mixed_model_rsq.R")
test<-lmer(arg~reca*log(bps)+(1+log(bps)|project)+(1+reca|project),data=final_results_m,REML=F)
test1<-lmer(arg~reca+log(bps)+(1+log(bps)|project)+(1+reca|project),data=final_results_m,REML=F)
test2<-lmer(arg~reca+(1+reca|project),data=final_results_m,REML=F)
test3<-lmer(arg~log(bps)+(1+log(bps)|project),data=final_results_m,REML=F)
anova(test2)
anova(test3)
r.squared.merMod(test2)
r.squared.merMod(test3)


anova(test)
summary(test)
AICctab(test,test2,test3)

r.squared.merMod(test)
r.squared.merMod(test2)
r.squared.merMod(test3)


ggplot(final_results)+geom_point(aes(x=reca,y=arg))



# making distance summaries
library(vegan)
head(reca_summary)
meta<-read.csv("./data/mgrast-allwgs.csv")
m<-merge(meta,reca_summary, by.x="id",by.y="variable")
head(meta)
head(m)
head(reca_summary)
count(meta$pt)
groups<-as.vector(unique(m$project))
length(groups)
head(count(reca_summary$variable))
head(reca_summary)
reca_dist_summary<-data.frame()
for(i in 1:length(groups)){
	
	temp<-subset(m, project==groups[i])

	reca_cast<-dcast(temp, id~reca, value.var="count",fun.aggregate=mean,fill=0)
	
	if(dim(reca_cast)[1] > 1){
		mat<-vegdist(decostand(reca_cast[,-1],"pa"),"jaccard")
		means<-mean(as.vector(mat))
		quants<-quantile(as.vector(mat),c(0.025,0.975))
		
		new_row<-data.frame(groups[i],means,quants[[1]],quants[[2]])
		reca_dist_summary<-rbind(reca_dist_summary,new_row)		
	}
}

arg_dist_summary<-data.frame()
for(i in 1:length(groups)){
	
	temp<-subset(m_arg, project==groups[i])
	
	arg_cast<-dcast(temp, id~ARG, value.var="count_arg",fun.aggregate=mean,fill=0)
	
	if(dim(arg_cast)[1] > 1){
		mat<-vegdist(decostand(arg_cast[,-1],"pa"),"jaccard")
		means<-mean(as.vector(mat))
		quants<-quantile(as.vector(mat),c(0.025,0.975))
		
		new_row<-data.frame(groups[i],means,quants[[1]],quants[[2]])
		arg_dist_summary<-rbind(arg_dist_summary,new_row)		
	}
}

m<-merge(reca_dist_summary,arg_dist_summary,by="groups.i.")
dim(m)
head(m)
as.vector(unique(m$project))
names(m)<-c("project","reca","reca_lo","reca_hi","arg","arg_lo","arg_hi")
env<-read.csv("./data/project_environment_key.csv")
head(m)
head(env)
as.matrix(env)

library(RColorBrewer)
length(unique(m$environment))
m<-merge(m,env,by="project")
m<-merge(m,bps, by="project")
display.brewer.all(8)
ggplot(subset(m, reca > 0))+geom_errorbar(aes(x=reca,y=arg,ymin=arg_lo,ymax=arg_hi),colour="grey",width=0)+geom_errorbarh(aes(x=reca,y=arg,xmin=reca_lo,xmax=reca_hi),colour="grey",height=0)+geom_point(aes(x=reca,y=arg,colour=log(bps/1000000)))+geom_abline(slope=1,intercept=0,linetype=2)+theme_bw(base_size=17)+theme(aspect.ratio=1)+scale_colour_gradientn(name="ln(Mbp Seq.)",colours=brewer.pal(8,"Spectral"))+labs(x="Microbiome Dissimilarity",y="Resistome Dissimilarity")


# checking with mantel test
reca_cast<-dcast(m_reca, project+name~reca, value.var="count_reca", fill=0,fun.aggregate=mean)
arg_cast<-dcast(m_arg, project+name~ARG, value.var="count_arg", fill=0,fun.aggregate=mean)
m_temp<-merge(reca_cast,arg_cast,by=c("project","name"))
?mantel
mantel(vegdist(decostand(m_temp[,c(3:70126)],"pa"),"jaccard"),vegdist(decostand(m_temp[,-c(1:70126)],"pa"),"jaccard"),method="spearman")