# here I am merging large matrices with metadata describing the samples from which they originate
args<-read.table("./full_args_cast.txt",header=T,sep="\t",check.names=F)
head(args[,1:10])
dim(args)
meta<-read.csv("./mgrast-allwgs.csv")
head(meta)
as.vector(unique(meta$pt))
str(meta)
args$dup<-round(as.numeric(as.character(args$dup)))

args_melt<-melt(args, id="dup")
str(args_melt)
str(meta)
args_melt_meta<-merge(args_melt,meta, by.x="dup",by.y="id_nopt")
head(args_melt_meta)

length(as.vector(unique(args_melt_meta$dup)))
# 1005 metagenomes included here

args_cast<-dcast(args_melt_meta, dup+project+bps~variable, value.var="value",fun.aggregate=sum, fill=0)
env_key<-read.csv("./project_environment_key.csv")
args_env_cast<-merge(args_cast, env_key, by="project",all.x=T)
args_no_env<-(subset(args_env_cast, is.na(environment)==T))
args_env<-(subset(args_env_cast, is.na(environment)==F))
args_no_env$environment<-"air"

args_env_cast<-rbind(args_env,args_no_env)
length(unique(args_env_cast$dup))
args_env_cast<-args_env_cast[,c(1:3,2981,4:2980)]
write.table(args_env_cast,"./args_environment_1007_metas.txt",sep="\t",row.names=F)
dim(args_env_cast)
head(args_env_cast[,1:10])


