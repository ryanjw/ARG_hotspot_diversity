cocur_data<-read.table("./co-occurrence_soil_hotspots_nov.txt",sep="\t",header=T)
head(cocur_data)

length(unique(c(cocur_data$Var1,cocur_data$Var2)))
cocur_data_pval<-subset(cocur_data, pval < 0.05)
hist(cocur_data$rho)
length(unique(c(cocur_data_pval$Var1,cocur_data_pval$Var2)))
var_key<-read.table("./cocur_variable_type_key.txt",sep="\t",header=T)
dim(var_key)

g<-graph.edgelist(as.matrix(subset(cocur_data, rho > 0.75)[,1:2]), directed=F)
# E(g)$weight<-subset(cocur_data, qval < 0.001)$rho
g<-simplify(g)
g
plot(g)
gnet<-asNetwork(g)
df<-asDF(g)
gnet