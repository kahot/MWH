library(zoo)
library(plotrix)
library(DescTools)
library(dplyr)


# Read the dataset file
cnsl<-read.csv("Data/cnsl_data.csv", stringsAsFactors = F)
cnsl$X<-paste0("No.",cnsl$X)
colnames(cnsl)[1]<-"Sample"

#calculate the average, average/2, and average*1.5 to determine copy nubmers
cnsl["Mean",3:ncol(cnsl)]<-colMeans(cnsl[,3:ncol(cnsl)])
cnsl["deletion",3:ncol(cnsl)]<-cnsl["Mean",3:ncol(cnsl)]/2
cnsl["duplication",3:ncol(cnsl)]<-cnsl["Mean",3:ncol(cnsl)]*1.5


#determine copy numbers at each probe for each sample
copyNo<-cnsl[1:10000,1:2]
copyNo[,3:ncol(cnsl)]<-""
colnames(copyNo)[3:ncol(cnsl)]<-colnames(cnsl)[3:ncol(cnsl)]

for (i in 1:(nrow(cnsl)-3)){
    for (k in 3:102){
        no<-k-2
        copyNo[i,k]<-which.min(c(abs(cnsl[i,k]-cnsl["deletion",k]), abs(cnsl[i,k]-cnsl["Mean",k]),
                                 abs(cnsl[i,k]-cnsl["duplication",k])))
    }
}
write.csv(copyNo,"Data/cnsl.CopyNo.csv")


#find CNVs based on the known breakpoints
bk<-data.frame(index=c(1:4),breakpt5=c(32,27,20,10),breakpt3=c(38,34,40,40))

bk.summary<-data.frame(Sample=copyNo$Sample[],ethnicity=copyNo$ethnicity[], Index1="", Index2="",Index3="",Index4="")
dt<-copyNo[,3:102]
for (k in 1:nrow(bk)){ 
    bk.summary[,(k+2)]<-rowMeans(dt[bk$breakpt5[k]:bk$breakpt3[k]])
}

#Create a summay table of CNVs by ethnicity
sum<-data.frame()
ethnic<-sort(unique(copyNo$ethnicity))  

for (i in 1:length(ethnic)) { 
    df<-bk.summary[bk.summary$ethnicity==ethnic[i],]
    tb<-data.frame(matrix(ncol=2, nrow=4))
    colnames(tb)<-c("deletion", "duplication")
    
    #proportion of deletion/duplication
    for (k in 1:4){
        tb$deletion[k]<-nrow(df[df[,paste0("Index",k)]==1,])/nrow(df)
        tb$duplication[k]<-nrow(df[df[,paste0("Index",k)]==3,])/nrow(df)
        tb$del.count[k]<-nrow(df[df[,paste0("Index",k)]==1,])
        tb$dup.count[k]<-nrow(df[df[,paste0("Index",k)]==3,])
        tb$total[k]<-nrow(df)
        tb$normal.count[k]<-tb$total[]
    }
    tb$normal.count<-tb$total-tb$del.count-tb$dup.count
    tb$index<-c("Index1","Index2","Index3","Index4")
    tb$ethnicity<-ethnic[i]
    sum<-rbind(sum,tb)
}   
write.csv(sum, "Output/Index1-4.CNV_freq_summary.csv")


# Test if the observed proportions of del/dup/normal differ by ethnicity for each index
Comb<-t(combn(c("A","B","C"), 2))
G.results<-data.frame()
for (k in 1:3){
    e1<-Comb[k,1]
    e2<-Comb[k,2]
    for (i in 1:4 ){
        dt<-data.frame(ethnic1=c(sum[sum$ethnicity==e1&sum$index==paste0("Index",i),3], 
                                 sum[sum$ethnicity==e1&sum$index==paste0("Index",i),4],
                                 sum[sum$ethnicity==e1&sum$index==paste0("Index",i),6]),                       
                       ethnic2=c(sum[sum$ethnicity==e2&sum$index==paste0("Index",i),3], 
                                 sum[sum$ethnicity==e2&sum$index==paste0("Index",i),4],
                                 sum[sum$ethnicity==e2&sum$index==paste0("Index",i),6]))
        result<-GTest(dt)
        G.results[paste0("Index",i),paste0(e1,"-",e2)]<-result[[3]][1]
    }
}
write.csv(G.results,"Output/Gtest.index.csv")


#calculate SD/SE for each probes
probe<-data.frame(probe=colnames(cnsl[3:102]))
for (i in 1: 100){
    no<-i+2
    probe$de.mean[i]<-mean(cnsl[copyNo[,no]==1,no])
    probe$de.se[i]<-std.error(cnsl[copyNo[,no]==1,no])
    probe$de.sd[i]<-sd(cnsl[copyNo[,no]==1,no])
    probe$du.mean[i]<-mean(cnsl[copyNo[,no]==3,no])
    probe$du.se[i]<-std.error(cnsl[copyNo[,no]==3,no])
    probe$du.sd[i]<-sd(cnsl[copyNo[,no]==3,no])
    probe$nor.mean[i]<-mean(cnsl[copyNo[,no]==2,no])
    probe$nor.se[i]<-std.error(cnsl[copyNo[,no]==2,no])
    probe$nor.sd[i]<-sd(cnsl[copyNo[,no]==2,no])
}

#normalize SE/SD per probe for comparison
probeS<-probe[,3:4]/probe[,2]
probeS[,3:4]<-probe[,6:7]/probe[,5]
probeS[,5:6]<-probe[,9:10]/probe[,8]
probeS$probe<-probe$probe

write.csv(probeS, "Output/Normalized.se.sd.per.probe.csv")

### Calculate the rolling average (4) of copy nubmers to determine the existence of dup/del at each probe
dt<-copyNo[,3:102]
rollm<-data.frame(matrix(ncol=47))
colnames(rollm)=colnames(dt)[1:47]
for (i in 1:nrow(dt)){ 
    for (j in 1: ncol(rollm)){
        vec<-t(dt[i,])
        rollm[i,]<-rollmean(t(dt[i,]), k=4, align = "left") 
    }
}
write.csv(rollm,"Output/rollmean.cnsl.csv")

#Deterine the frequency of del/dup at each probe per ethnicity
rollm$ethnicity<-copyNo$ethnicity
Summary<-data.frame()
ethnic<-sort(unique(copyNo$ethnicity)) 
for (j in 1:3){
    df<-rollm[rollm$ethnicity==ethnic[j],1:47]
    n<-nrow(df)
    
    tb<-data.frame(probe=colnames(rollm)[1:47])
    for (i in 1:ncol(df)){
        tb$deletion[i]   <-length(df[df[,i]==1,i])/n
        tb$duplication[i]<-length(df[df[,i]==3,i])/n
    }
    tb$ethnicity<-ethnic[j]
    Summary<-rbind(Summary,tb)
}
Summary$probe<-factor(Summary$probe, levels=c(paste0("CNSL_probe_",0:49)) )

write.csv(Summary, "Output/Summary_CNVfreq_perProbe.csv")


########### 
# Find breakpoints for each samples 
copy2<-copyNo[,3:52]
SumDF<-data.frame()

for (i in 1:nrow(copy2)){
    DF<-rle(copy2[i,]) 
    df<-as.data.frame(t(DF[[2]]))
    colnames(df)<-"CopyNo"
    df$probe<-as.integer(substr(row.names(df), start = 12, stop = 13))
    df$length<-DF[[1]]
    
    df<-df[df$CopyNo!=2 & df$length>3,]
    if (nrow(df)!=0){
        df$min<-df$probe-df$length+1
        df$max<-df$probe
        df$Sample<-i
        SumDF<-rbind(SumDF, df)
    }
}

SumDF$CopyNo<-as.character(SumDF$CopyNo)
for (i in 1:nrow(SumDF)){
    SumDF$ethnicity[i]<-copyNo$ethnicity[SumDF$Sample[i]]
}
write.csv(SumDF, "Output/Breakpoints_Summary.csv")


## Create freq summary tables for each breakpoint: 

#1. 5' breakpoints
#deletion and duplication together (for plotting)

bk5<-data.frame()
for ( i in 1:3){
    tb<-data.frame(table(SumDF$min[SumDF$ethnicity==ethnic[i]]))
    tb$freq<-tb$Freq/nrow(copyNo[copyNo$ethnicity==ethnic[i],])
    tb$ethnicity<-ethnic[i]
    bk5<-rbind(bk5,tb)
}
write.csv(bk5, "Output/bk5.csv")

#1.2. dup and del separately
bk5.2<-data.frame()
for ( i in 1:3){
    tb1<-data.frame(table(SumDF$min[SumDF$CopyNo=="1" & SumDF$ethnicity==ethnic[i]]))
    tb1$type<-"deletion"
    tb2<-data.frame(table(SumDF$min[SumDF$CopyNo=="3" & SumDF$ethnicity==ethnic[i]]))
    tb2$type<-"duplication"
    tb<-rbind(tb1,tb2)
    tb$freq<-tb$Freq/nrow(copyNo[copyNo$ethnicity==ethnic[i],])
    tb$ethnicity<-ethnic[i]
    bk5.2<-rbind(bk5.2,tb)
}
write.csv(bk5.2, "Output/bk5.2.csv")

#Test for significance
for ( i in 1:3){
    tb1<-data.frame(table(SumDF$min[SumDF$CopyNo=="1" & SumDF$ethnicity==ethnic[i]]))
    colnames(tb1)<-c("probe","deletion")
    tb2<-data.frame(table(SumDF$min[SumDF$CopyNo=="3" & SumDF$ethnicity==ethnic[i]]))
    colnames(tb2)<-c("probe","duplication")
    tb<-merge(tb1,tb2, by="probe",all = T)
    tb$normal<-nrow(copyNo[copyNo$ethnicity==ethnic[i],])-tb$deletion-tb$duplication
    tname<-paste0("bk5.",ethnic[i])
    assign(tname, tb)
}

Comb<-t(combn(c("A","B","C"), 2))
Gresults<-list()
for (k in 1:3){
    e1<-Comb[k,1]
    e2<-Comb[k,2]
    tb1<-get(paste0("bk5.",e1))
    tb2<-get(paste0("bk5.",e2))
    tb<-merge(tb1, tb2, by="probe", all=T)
    for (i in 1:nrow(tb) ){
        dt<-data.frame(ethnic1=c(tb$deletion.x[i], tb$duplication.x[i], tb$normal.x[i]),
                       ethnic2=c(tb$deletion.y[i], tb$duplication.y[i], tb$normal.y[i]))
        result<-GTest(dt)
        tb$G.results[i]<-result[[3]][1]
    }
    Gresults[[k]]<-tb
    names(Gresults)[k]<-paste0(e1, ".vs.",e2)
}

sink("Output/G-test.results.bk5.txt")
for (k in 1:3){
    dt<-Gresults[[k]]
    sig<-dt$probe[which(dt$G.results<0.05)]
    cat("Significant probe numbers ",names(Gresults)[k],"\n" )
    print(paste(sig))
}
sink(NULL)


#2. 3'breakpoints
bk3<-data.frame()
for ( i in 1:3){
    tb<-data.frame(table(SumDF$max[SumDF$ethnicity==ethnic[i]]))
    tb$freq<-tb$Freq/nrow(copyNo[copyNo$ethnicity==ethnic[i],])
    tb$ethnicity<-ethnic[i]
    bk3<-rbind(bk3,tb)
}
write.csv(bk3, "Output/bk3.csv")

#2.2 dep and del separately
bk3.2<-data.frame()
for ( i in 1:3){
    tb1<-data.frame(table(SumDF$max[SumDF$CopyNo=="1" & SumDF$ethnicity==ethnic[i]]))
    tb1$type<-"deletion"
    tb2<-data.frame(table(SumDF$max[SumDF$CopyNo=="3" & SumDF$ethnicity==ethnic[i]]))
    tb2$type<-"duplication"
    tb<-rbind(tb1,tb2)
    tb$freq<-tb$Freq/nrow(copyNo[copyNo$ethnicity==ethnic[i],])
    tb$ethnicity<-ethnic[i]
    bk3.2<-rbind(bk3.2,tb)
}
write.csv(bk3.2, "Output/bk3.2.csv")

#Test for significance: bk3
for ( i in 1:3){
    tb1<-data.frame(table(SumDF$max[SumDF$CopyNo=="1" & SumDF$ethnicity==ethnic[i]]))
    colnames(tb1)<-c("probe","deletion")
    tb2<-data.frame(table(SumDF$max[SumDF$CopyNo=="3" & SumDF$ethnicity==ethnic[i]]))
    colnames(tb2)<-c("probe","duplication")
    tb<-merge(tb1,tb2, by="probe",all = T)
    tb$normal<-nrow(copyNo[copyNo$ethnicity==ethnic[i],])-tb$deletion-tb$duplication
    tname<-paste0("bk3.",ethnic[i])
    assign(tname, tb)
}

Comb<-t(combn(c("A","B","C"), 2))
Gresults2<-list()
for (k in 1:3){
    e1<-Comb[k,1]
    e2<-Comb[k,2]
    tb1<-get(paste0("bk3.",e1))
    tb2<-get(paste0("bk3.",e2))
    tb<-merge(tb1, tb2, by="probe", all=T)
    for (i in 1:nrow(tb) ){
        dt<-data.frame(ethnic1=c(tb$deletion.x[i], tb$duplication.x[i], tb$normal.x[i]),
                       ethnic2=c(tb$deletion.y[i], tb$duplication.y[i], tb$normal.y[i]))
        result<-GTest(dt)
        tb$G.results[i]<-result[[3]][1]
    }
    Gresults2[[k]]<-tb
    names(Gresults2)[k]<-paste0(e1, ".vs.",e2)
}

sink("Output/G-test.results.bk3.txt")
for (k in 1:3){
    dt<-Gresults2[[k]]
    sig<-dt$probe[which(dt$G.results<0.05)]
    cat("Significant probe numbers ",names(Gresults)[k],"\n" )
    print(paste0(sig))
}
sink(NULL)


