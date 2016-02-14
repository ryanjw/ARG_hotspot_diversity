library(indicspecies)
args_env_cast<-read.table("./args_environment_1007_metas.txt",sep="\t",header=T,check.names=F)
head(args_env_cast[,1:10])

soil_hotspots<-(subset(args_env_cast, dup==4537095 |
dup==4521317 |
dup==4573678 |
dup==4554153 |
dup==4573679 |
dup==4573682 |
dup==4575385 |
dup==4554152 |
dup==4573680 |
dup==4554159 |
dup==4520326 |
dup==4573683 |
dup==4520302 |
dup==4520316 |
dup==4520310 ))

soil_notspots<-subset(args_env_cast, environment=="soil" & dup!=4537095 &
dup!=4521317 &
dup!=4573678 &
dup!=4554153 &
dup!=4573679 &
dup!=4573682 &
dup!=4575385 &
dup!=4554152 &
dup!=4573680 &
dup!=4554159 &
dup!=4520326 &
dup!=4573683 &
dup!=4520302 &
dup!=4520316 &
dup!=4520310 )
dim(soil_notspots)

category<-c(rep("hot",15),rep("not",496))

soils_cats<-cbind(category,rbind(soil_hotspots, soil_notspots))
head(soils_cats[,1:10])

# removing columns with all 0's and singletons
mat<-soils_cats[,-c(1:5)]
mat<-mat[,colSums(decostand(mat, "pa")) > 0]
soils_cats_slim<-data.frame(soils_cats[,1:5],mat)
dim(soils_cats_slim)

#now performing indicator species analysis
results<-multipatt(soils_cats_slim[,-c(1:5)],soils_cats_slim$category,print.perm=T)
write.table(soils_cats_slim, "./soils_categories_for_indicator_species_analysis_nov.txt",sep="\t",row.names=F)

# reading back in file (it was run remotely)
results<-read.table("./soils_indicators.txt",sep="\t",check.names=F,header=T)
head(results)
results<-na.omit(results)
results$qval<-fdrtool(results$p.value,statistic="pvalue",plot=F,verbose=F)$qval
results$ARG<-rownames(results)
dim(subset(results, s.hot==1 & p.value==0.005))