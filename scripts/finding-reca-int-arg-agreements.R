# HMP analysis was not included in this paper

args<-read.table("~/Desktop/big_ARG_data/refgenome_data/hmp-arg.summary",sep="\t",header=T,check.names=F)
head(args[,1:10])
args<-row_column_correction(args,"ref")
args_melt<-melt(args, id="ref")
args_melt$var_type<-"arg"
args_melt<-subset(args_melt, variable=="AAK64454"|
variable=="AAP46617"|
variable=="AAQ12900"|
variable=="AAV74563"|
variable=="AAY62597"|
variable=="ABN79417"|
variable=="ACA23192"|
variable=="YP_001691306"|
variable=="YP_001693856"|
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

ints<-read.table("~/Desktop/big_ARG_data/refgenome_data/hmp-int.summary",sep="\t",header=T,check.names=F)
head(ints[,1:10])
ints<-row_column_correction(ints,"ref")
ints_melt<-melt(ints, id="ref")
ints_melt$var_type<-"int"
head(ints_melt)
recs<-read.table("~/Desktop/big_ARG_data/refgenome_data/hmp-rec.summary",sep="\t",header=T,check.names=F)
head(recs[,1:10])
recs<-row_column_correction(recs,"ref")
recs_melt<-melt(recs, id="ref")
recs_melt$var_type<-"rec"
head(recs_melt)
m<-rbind(args_melt,ints_melt,recs_melt)
head(subset(m, var_type=="int"))

# readin in the network
g_list<-read.table("~/Desktop/big_ARG_data/g_list.txt",header=T,sep="\t")
head(g_list)

# these are the ints that are important
as.matrix(unique(subset(g_list, var1type=="tr" & var2type=="sp")$Var1))


#subset these out of the ints list

ints_sub<-subset(ints_melt, variable== "tr|A0A0D5ZYN3|A0A0D5ZYN3_ACIBA"|
variable== "tr|A3RLT1|A3RLT1_SALET"        |
variable== "tr|B2BE20|B2BE20_9BACT"        |
variable== "tr|C6GC53|C6GC53_9BACT"        |
variable== "tr|C6GC77|C6GC77_9BACT"        |
variable== "tr|C6GC91|C6GC91_9BACT"        |
variable== "tr|W0RY67|W0RY67_KLEPN"
& value > 0)
length(unique(ints_melt$ref))
# makes no difference here


# doing the same thing with reca since this should slim it down
as.matrix(unique(subset(g_list, var1type=="tr" & var2type=="sp")$Var2))

recs_sub<-subset(recs_melt,variable== "tr|Q8DNI9|SCPB_STRR6"  |
variable== "tr|O83964|FTSK_TREPA"  |
variable== "tr|P0C6C7|FLAE_VIBCH"  |
variable== "tr|P38086|RDH54_YEAST" |
variable== "tr|P39920|FTSK_COXBU"  |
variable== "tr|Q51564|DAPF_PSEAE"  |
variable== "tr|Q87ZS5|FTSK_PSESM"  |
variable== "tr|Q88FS8|FTSK_PSEPK"  |
variable== "tr|Q895I8|FTSK_CLOTE"  |
variable== "tr|Q8R5S4|FTSK_CALS4"  |
variable== "tr|Q8XJS8|FTSK_CLOPE"  |
variable== "tr|Q8YJB8|FTSK_BRUME"  |
variable== "tr|Q92L89|FTSK_RHIME"  |
variable== "tr|Q97I41|FTSK_CLOAB"  |
variable== "tr|Q98EH3|FTSK_RHILO"  |
variable== "tr|Q9JZG4|FTSK2_NEIMB" |
variable== "tr|Q9RNV1|FTSK_SPOUR"  |
variable== "tr|Q9Z429|NAHY_PSEPU"  |
variable== "tr|O05321|DCDA_PSEFL"  |
variable== "tr|O83045|FTSK_AZOBR"  |
variable== "tr|O86810|FTSK_STRCO"  |
variable== "tr|P0A6H1|CLPX_ECOLI"  |
variable== "tr|Q48735|HSLU_LACLE"  |
variable== "tr|Q706S1|MHPD3_PSEPU" |
variable== "tr|Q82K93|FTSK_STRAW"  |
variable== "tr|Q89WR2|FTSK_BRADU"  |
variable== "tr|Q8G4H3|FTSK_BIFLO"  |
variable== "tr|Q8U526|FTSK_AGRT5"  |
variable== "tr|Q9A262|FTSK_CAUCR"  |
variable== "tr|Q9EUQ7|SCPA_STRR6"  |
variable== "tr|Q9EUR1|SCPB_STRMT"  |
variable== "tr|O05322|DAPF_PSEFL"  |
variable== "tr|O94762|RECQ5_HUMAN" |
variable== "tr|Q81A03|FTSK_BACCR"  |
variable== "tr|Q8EQS7|FTSK_OCEIH"  |
variable== "tr|Q8PL00|FTSK_XANAC"  |
variable== "tr|Q8XRH0|FTSK1_RALSO" |
variable== "tr|Q8XWX9|FTSK2_RALSO" |
variable== "tr|P0C6C6|FLAD_VIBCH"  |
variable== "tr|Q8ZGC7|FTSK_YERPE"  |
variable== "tr|Q9I0M3|FTSK_PSEAE" & value > 0)

dim(recs_sub)

# none of these recas showed up in HMP

# attempting to merge directly into g_list based on int and arg
m2<-merge(subset(g_list, var1type=="tr" & var2type=="arg"),subset(args_melt, value > 0),by.x="Var2",by.y="variable")
has_args<-dcast(m2, ref+Var1~Var2, value.var="value",fun.aggregate=sum, fill=0)
head(has_args)
has_args_melt<-subset(melt(has_args, id=c("ref","Var1")),value > 0)
names(has_args_melt)[3]<-"Var2"

has_args_agreement<-merge(has_args_melt, g_list, by=c("Var1","Var2") )
df<-data.frame(unique(has_args_agreement$ref))

names(df)<-"HMP_references"
head(df)
write.table(df, "~/Desktop/HMP_references_with_args_and_ints.txt",row.names=F,sep="\t")

# now doing the same thing with refsoil

args<-read.table("~/Desktop/big_ARG_data/refgenome_data/soil-arg.summary",header=T, sep="\t",check.names=F)
args<-row_column_correction(args,"ref")
dim(args)
head(args[,1:10])
args_melt<-melt(args, id="ref")
args_melt$var_type<-"arg"
head(args_melt)
subset(args_melt, ref=="NC_011884" & value > 0)
length(unique(args_melt$ref))
args_soil_melt<-subset(args_melt, variable =="AAK64454"|
variable =="AAP46617"|
variable =="AAQ12900"|
variable =="AAV74563"|
variable =="AAY62597"|
variable =="ABN79417"|
variable =="ACA23192"|
variable =="YP_002242988"|
variable =="YP_002850805"|
variable =="YP_407313"|
variable =="YP_788911"|
variable =="ZP_03029420"|
variable =="ZP_03063464"|
variable =="ZP_03452771"|
variable =="ZP_04466634"|
variable =="ZP_04557016"|
variable =="ZP_04576183"|
variable =="ZP_04623897")
ints<-read.table("~/Desktop/big_ARG_data/refgenome_data/soil-int.summary",sep="\t",header=T,check.names=F)
head(ints[,1:10])
ints<-row_column_correction(ints,"ref")
ints_melt<-melt(ints, id="ref")
ints_melt$var_type<-"int"
head(ints_melt)
recs<-read.table("~/Desktop/big_ARG_data/refgenome_data/soil-rec.summary",sep="\t",header=T,check.names=F)
head(recs[,1:10])
recs<-row_column_correction(recs,"ref")
recs_melt<-melt(recs, id="ref")
recs_melt$var_type<-"rec"
m<-rbind(args_melt,ints_melt,recs_melt)

head(m)
head(g_list)

head(args_melt)
g_list_refsoil_m1<-merge(subset(g_list, var1type=="arg" & var2type=="tr"), subset(args_melt, value > 0),by.x="Var1",by.y="variable")
head(g_list_refsoil_m1)
head(ints_melt)
names(ints_melt)[2]<-"Var2"
g_list_refsoil_m2<-merge(g_list_refsoil_m1, subset(ints_melt, value > 0),by=c("ref","Var2"))
head(g_list_refsoil_m2)
taxonomy<-read.table("~/Desktop/big_ARG_data/refgenome_data/taxa2",sep="\t",header=F)
head(taxonomy)
df.ref_soil<-data.frame(unique(g_list_refsoil_m2$ref))
names(df.ref_soil)[1]<-"ref"
df_m<-merge(df.ref_soil,taxonomy, by.x="ref",by.y="V1",all.x=TRUE)
df_m
names_too<-merge(g_list_refsoil_m2,taxonomy,by.x="ref",by.y="V1",all.x=TRUE)
write.table(names_too,"./data/confirmed_refsoil_co-occurrence_relationships.txt",sep="\t",row.names=FALSE)
# P. aeroginosa showed up again

head(recs_melt)
subset(recs_melt, ref=="NC_013891" & value > 0)
subset(g_list_refsoil_m2, ref=="NC_009656" )

g_list_refsoil_m2$pairs<-paste(g_list_refsoil_m2$Var2,g_list_refsoil_m2$Var1)
arrange(ddply(g_list_refsoil_m2, .(pairs),summarise, count=length(unique(ref))),-count)
g_list$pairs<-paste(g_list$Var2,g_list$Var1)
length(unique(subset(g_list, var1type=="arg" & var2type=="tr")$pairs))
length(unique(g_list_refsoil_m2$pairs))
subset(g_list_refsoil_m2, pairs=="tr|C6GC53|C6GC53_9BACT YP_002242988")$ref
