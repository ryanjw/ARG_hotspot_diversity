# source these functions first
source("../scripts/mixed_model_rsq.R")
library(data.table)
args_env_cast<-fread("./data/merged_arg_with_meta_final.txt",sep="\t",header=T,data.table=F)
head(args_env_cast)
args_env_cast<-dcast(args_env_cast, name+project+bps+environment~ARG,value.var="count",fill=0,fun.aggregate=mean)
totals<-rowSums(args_env_cast[,-c(1:4)])
alpha<-specnumber(args_env_cast[,-c(1:4)])

totals_norm<-rowSums(args_env_cast[,-c(1:4)])/args_env_cast$bps
alpha_norm<-specnumber(args_env_cast[,-c(1:4)])/args_env_cast$bps

args_stats<-data.frame(args_env_cast[,c(1:4)],totals,alpha,totals_norm,alpha_norm)
head(args_stats)

# Does diversity of ARGs change as abundance increases?
test1<-lmer(log1p(alpha)~log1p(totals)*log(bps)+(1+log1p(totals)|project), REML=F, data=args_stats)
test2<-lmer(log1p(alpha)~log1p(totals)+log(bps)+(1+log1p(totals)|project), REML=F, data=args_stats)
test3<-lmer(log1p(alpha)~log1p(totals)+(1+log1p(totals)|project), REML=F, data=args_stats)
AICctab(test1,test2,test3)
      # dAICc df
# test1  0.0  8 
# test2 64.6  7 
# test3 67.8  6 


anova(test1)
r.squared.merMod(test1)
# Analysis of Variance Table of type III  with  Satterthwaite 
# approximation for degrees of freedom
                        # Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
# log1p(totals)          15.4952 15.4952     1 527.39  326.23 < 2.2e-16 ***
# log(bps)                3.6674  3.6674     1 734.60   77.21 < 2.2e-16 ***
# log1p(totals):log(bps)  3.6841  3.6841     1 632.38   77.56 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > r.squared.merMod(test1)
           # Class   Family     Link  Marginal Conditional      AIC
# 1 merModLmerTest gaussian identity 0.8960004   0.9792609 200.9491

anova(test3)
r.squared.merMod(test3)
# > anova(test3)
# Analysis of Variance Table of type III  with  Satterthwaite 
# approximation for degrees of freedom
              # Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
# log1p(totals) 51.355  51.355     1 99.841  1044.7 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > r.squared.merMod(test3)
           # Class   Family     Link  Marginal Conditional      AIC
# 1 merModLmerTest gaussian identity 0.8680107    0.979593 268.6565

# The additional effect is minimal, so we normalize by bps

head(args_stats)
ggplot(args_stats)+geom_point(aes(x=totals,y=alpha,colour=log(bps)))+scale_colour_gradientn(name="ln(Seq. Depth)",colours=c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#74add1","#4575b4","#313695"))+theme_bw(base_size=15)+theme(aspect.ratio=1)+labs(x="ARG Abundance",y="ARG Diversity")+scale_x_log10()+scale_y_log10()

ggplot(args_stats)+geom_point(aes(x=totals_norm,y=alpha_norm,colour=log(bps)))+scale_colour_gradientn(name="ln(Seq. Depth)",colours=c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#74add1","#4575b4","#313695"))+theme_bw(base_size=15)+theme(aspect.ratio=1)+labs(x="Normalized ARG Abundance",y="Normalized ARG Diversity")+scale_x_log10()+scale_y_log10()

# all of these are correlated, so we will normalize by seq depth

# Does diversity of ARGs change as abundance increases after normalization?
test1<-lmer(log1p(alpha_norm)~log1p(totals_norm)+(1+log1p(totals_norm)|project), REML=F, data=args_stats)
test2<-lmer(scale(alpha_norm)~scale(totals_norm)+(1+scale(totals_norm)|project), REML=F, data=args_stats)
test3<-lmer((alpha_norm)~(totals_norm)+(1+(totals_norm)|project), REML=F, data=args_stats)
test1
AICctab(test1,test2)
      # dAICc   df
# test1     0.0 6 
# test2 21136.5 6 

# model had issues converging for test1, so had to go with test2
anova(test2)
r.squared.merMod(test2)

# Analysis of Variance Table of type III  with  Satterthwaite 
# approximation for degrees of freedom
                   # Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
# scale(totals_norm) 1.0497  1.0497     1 84.822  23.497 5.597e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > r.squared.merMod(test2)
           # Class   Family     Link Marginal Conditional       AIC
# 1 merModLmerTest gaussian identity 0.340537   0.9970885 -60.90137

ggplot(args_stats)+geom_point(aes(x=totals_norm,y=alpha_norm,colour=log(bps)))+scale_colour_gradientn(name="ln(Seq. Depth)",colours=c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#74add1","#4575b4","#313695"))+theme_bw(base_size=15)+theme(aspect.ratio=1)+labs(x="ARG Abundance",y="ARG Diversity")+scale_x_log10()+scale_y_log10()

# normalization shows that the relationship between ARG diversity and abundance are largely decoupled when accounting for bps sequenced.  We will now incorporate intI1

ints_cast<-fread("./data/merged_int_with_meta_final.txt",sep="\t",header=T,data.table=F)
# removing bad ints
funky_ints<-read.csv("./data/funky-ints.csv",header=F)
ints_cast<-ints_cast[! ints_cast$IntI1 %in% funky_ints$V1,]
ints_cast<-dcast(ints_cast, name+project+bps+environment~IntI1,value.var="count",fill=0,fun.aggregate=mean)

ints_totals<-rowSums(ints_cast[,-c(1:4)])
ints_alpha<-specnumber(ints_cast[,-c(1:4)])

ints_stats<-data.frame(ints_cast$name, ints_totals,ints_alpha)
names(ints_stats)[1]<-"name"

str(ints_stats)
str(args_stats)
arg_ints_m<-merge(args_stats, ints_stats, by="name")
arg_ints_m$ints_totals_norm<-arg_ints_m$ints_totals/arg_ints_m$bps
arg_ints_m$ints_alpha_norm<-arg_ints_m$ints_alpha/arg_ints_m$bps
head(arg_ints_m)
tail(arg_ints_m)
# making the triangle plot based on quantiles

arg_ints_m<-arrange(arg_ints_m,totals_norm)
totals_q<-quantile(arg_ints_m$totals_norm,probs=seq(0,1,length.out=dim(arg_ints_m)[1]))
attr<-attributes(totals_q)
arg_ints_m$totals_q<-as.numeric(substr(attr$names,0,nchar(attr$names)-1))

arg_ints_m<-arrange(arg_ints_m,alpha_norm)
alpha_q<-quantile(arg_ints_m$alpha_norm,probs=seq(0,1,length.out=dim(arg_ints_m)[1]))
attr<-attributes(alpha_q)
arg_ints_m$alpha_q<-as.numeric(substr(attr$names,0,nchar(attr$names)-1))

arg_ints_m<-arrange(arg_ints_m,ints_totals_norm)
ints_totals_norm_q<-quantile(arg_ints_m$ints_totals_norm,probs=seq(0,1,length.out=dim(arg_ints_m)[1]))
attr<-attributes(ints_totals_norm_q)
arg_ints_m$ints_totals_norm_q<-as.numeric(substr(attr$names,0,nchar(attr$names)-1))

#writing to table to predict errors with environment labels manually
write.table(arg_ints_m,"./data/arg_ints_m_final.txt",sep="\t",row.names=F)
arg_ints_m<-read.table("./data/arg_ints_m_final_id.txt",sep="\t",header=TRUE)
tail(arg_ints_m)
ggtern::ggtern(arg_ints_m)+geom_point(aes(x=totals_q,y=alpha_q,z=ints_totals_norm_q,colour=environment),alpha=0.75)+labs(x="AMR Abundance",y="AMR Diversity",z="intI1 Abundance")+theme_bw(base_size=14)+scale_colour_manual(name="Environment",values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#f781bf","#a65628","#999999"))

# What are the summary statistics for these different hotspot types
head(arg_ints_m)
hotspot_stats<-ddply(arg_ints_m, .(environment),summarise,
q_total=mean(totals_q),
total=mean(totals_norm)*1000000,
total_sd=sd(totals_norm)*1000000,
q_alpha=mean(alpha_q),
alpha=mean(alpha_norm)*1000000,
alpha_sd=sd(alpha_norm)*1000000,
q_int=mean(ints_totals_norm_q),
qs=mean(ints_totals_norm)*1000000,
int_sd=sd(ints_totals_norm)*1000000
)
head(arrange(hotspot_stats, -q_total))
head(arrange(hotspot_stats, -q_alpha))
arrange(hotspot_stats, -q_int)






# making df based on hot or not
head(df)
df$ab_hot<-0
df$div_hot<-0
df$int_hot<-0

x<-arg_ints_m$totals_norm
psx<-pnorm(log(x), mean=mean(log(x)),sd=sd(log(x)))
y<-arg_ints_m$alpha_norm
psy<-pnorm(log(y), mean=mean(log(y)),sd=sd(log(y)))
z<-arg_ints_m$ints_totals_norm
psz<-pnorm(log1p(z), mean=mean(log1p(z)),sd=sd(log1p(z)))
df<-data.frame(arg_ints_m,x,y,z,psx,psy,psz)

ddply(df, .(environment),summarise, count=length(unique(dup)))
df_sub<-subset(df, environment=="freshwater habitat" | environment=="marine habitat" | environment=="soil" | environment=="waste water/sludge")
envs<-as.vector(unique(df$environment))
head(df_sub)
df_final<-data.frame()
for(i in 1:length(envs)){
	
	temp<-subset(df, environment==envs[i])
	x<-temp$totals_norm
	psx<-pnorm((x), mean=mean((x)),sd=sd((x)))
	y<-temp$alpha_norm
	psy<-pnorm((y), mean=mean((y)),sd=sd((y)))
	z<-temp$ints_totals_norm
	psz<-pnorm((z), mean=mean((z)),sd=sd((z)))
	df_temp<-data.frame(temp[,1:4],x,y,z,psx,psy,psz)
	df_temp$ab_hot<-0
	df_temp$div_hot<-0
	df_temp$int_hot<-0
	for(i in 1:dim(df_temp)[1]){
		if(df_temp$psx[i] > 0.95){df_temp$ab_hot[i]<-1}
		if(df_temp$psy[i] > 0.95){df_temp$div_hot[i]<-1}
		if(df_temp$psz[i] > 0.95){df_temp$int_hot[i]<-1}
		
	}
	df_final<-rbind(df_final, df_temp)
	
}

# supplemental table 1 info along with above stats (using x,y,z)
head(df_final)
final_stats<-ddply(df_final, .(environment), summarise,
m=mean(z)*1000000,
sd=sd(z)*1000000,
lo=quantile(z, 0.025)*1000000,
hi=quantile(z, 0.975)*1000000,
min=min(z)*1000000,
max=max(z)*1000000,
cv=sd/m
)
final_stats

head(df_sub)


# looking at abundance of indicator ARGs across hotspots and not hotspots
args_env_cast<-read.table("./args_environment_1007_metas.txt",sep="\t",header=T)
head(df) #from before
head(args_env_cast[,1:10])
args_env_cast_melt<-melt(args_env_cast, id=c("project","dup","bps","environment"))
summary(df_sub)
big_melt<-merge(df_final,args_env_cast_melt,by="dup")
head(big_melt)


big_melt_sub<-subset(big_melt, variable=="AAK64454"|
variable=="AAP46617"|
variable=="AAQ12900"|
variable=="AAV74563"|
variable=="AAY62597"|
variable=="ABN79417"|
variable=="ACA23192"|
variable=="YP_002242988"|
variable=="YP_002850805"|
variable=="YP_407313"|
variable=="YP_788911"|
variable=="ZP_03029420"|
variable=="ZP_03063464"|
variable=="ZP_03452771"|
variable=="ZP_04466634"|
variable=="ZP_04557016"|
variable=="ZP_04576183"|
variable=="ZP_04623897")

hot_key<-read.csv("./soils_hot_key.csv")
big_melt_sub_m<-merge(big_melt_sub, hot_key,by="variable")
head(big_melt_sub_m)
summary(big_melt_sub_m$ab_ho)
hot_arg_stats<-ddply(big_melt_sub_m, .(antibiotic, resistance_type,ab_hot),summarise,
m=(mean(value)),
SE=(sd(value)/sqrt(length(value)))
)

hot_arg_stats
ggplot(subset(hot_arg_stats))+geom_pointrange(aes(x=reorder(resistance_type,m),y=m,ymin=m-SE,ymax=m+SE,colour=as.factor(ab_hot)),position=position_dodge(width=1),shape=32)+geom_bar(aes(x=reorder(resistance_type,m),y=m,fill=as.factor(ab_hot)),stat="identity",position=position_dodge(width=1))+coord_flip()+labs(x="Resistance Genes",y="Abundance")+theme_bw(base_size=14)+theme(aspect.ratio=1,axis.text.y=element_text(face="italic"))+scale_fill_manual(name="Hotspot\nIdentity",values=c("black","red"))+scale_colour_manual(values=c("black","red"))+guides(colour=F)
comp_stats<-dcast(hot_arg_stats, antibiotic+resistance_type~ab_hot, value.var="m",fun.aggregate=mean)
comp_stats$diff<-comp_stats[,4]/comp_stats[,3]
arrange(comp_stats, -diff)