library(reshape2)
library(ggplot2)
library(colorspace)
library(ggstance)

colors<-qualitative_hcl(6, palette="Dark3")

cnsl<-read.csv("Data/cnsl_data.csv", stringsAsFactors = F)
cnsl$X<-paste0("No.",cnsl$X)
colnames(cnsl)[1]<-"Sample"

copyNo<-read.csv("Data/cnsl.CopyNo.csv", stringsAsFactors = F, row.names = 1)

#Plot the CNVs based on the known breakpoints
sum<-read.csv("Output/Index1-4.CNV_freq_summary.csv", stringsAsFactors = F, row.names = 1)
sum.m<-melt(sum[,c("deletion","duplication","index","ethnicity")])

ggplot(sum.m, aes(x=index, y=value, shape=variable,color=ethnicity))+
    geom_point(position=position_dodge(width=0.3),size =3)+
    geom_path(aes(x=index, y=value,group=interaction(ethnicity,variable)),position=position_dodge(width=0.3),linetype = 1, size=0.2)+
    scale_color_manual(values=colors[c(1,3,5)])+
    scale_shape_manual(values=c(19,17), name="",labels=c("Deletion", "Duplication"))+
    xlab("Del/dup index")+ylab("Observed frequency")+
    theme_bw()+
    theme(axis.text.x = element_text(size=12),axis.title.y = element_text(size=12))+
    theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank() )+
    geom_vline(xintercept = c(1:3)+0.5, color = "gray70", size=.3)+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))
ggsave("Output/Summary_plot.pdf", width = 5.5, height = 4.5)


# Plot normalized SD/SE per probe
probeS<-read.csv("Output/Normalized.se.sd.per.probe.csv", row.names = 1)
probeS$probe<-factor(probeS$probe, levels=c(paste0("CNSL_probe_",0:49),paste0("non_CNSL_probe_",0:49)) )

probeSm<-melt(probeS[,c(2,4,6,7)])

ggplot(probeSm, aes(x=probe, y=value, color=variable))+
    geom_point(size =0.6)+
    scale_color_manual(values=colors[c(1,4,5)], name="",labels=c("Deletion", "Duplication","Normal"))+
    theme(axis.title.x=element_blank())+ylab("Normalized SD")+
    theme_bw()+
    theme(axis.text.x = element_text(size = 6, angle = 90))+
    xlab("")
ggsave(filename="Output/probe.sd.pdf",width = 11, height = 2.5)


# Plot CNV proportions per ethicity
df<-sum[,c(1,2,7,8)]
df$normal<-sum$normal.count/sum$total
df.m<-melt(df)

ggplot(df.m, aes(x=ethnicity , y=value, fill=variable))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=paste0(colors[c(1,5,4)],"CC"), name="",labels=c("Deletion", "Duplication", "Normal"))+
    facet_grid(~index )+
    theme_bw()+
    ylab("Observed proportion")
ggsave("Output/Summary_proportion.by.ethinicity.pdf", width = 6, height = 6)

#Plot del/dup ratios
df$ratio<-df$duplication/df$deletion
ggplot(df, aes(x=index , y=ratio, fill=ethnicity))+
    geom_bar(stat="identity", position=position_dodge())+
    scale_fill_manual(values=paste0(colors[c(1,3,5)],"CC"), name="Ethnicity")+
    theme_bw()+
    ylab("Duplication to deletion ratio")+
    xlab("Del/dup index")+
    theme(panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank() )+
    geom_vline(xintercept = c(1:3)+0.5, color = "gray70", size=.3)
ggsave("Output/Ratio.by.ethinicity.pdf", width = 5, height = 4.2)



### Rolling average (4) of copy nubmers to determine the existence of dup/del at each probe
#rollm<-read.csv("Output/rollmean.cnsl.csv", stringsAsFactors = F, row.names = 1)
Summary<-read.csv("Output/Summary_CNVfreq_perProbe.csv", row.names = 1)
Summary$probe<-factor(Summary$probe, levels=c(paste0("CNSL_probe_",0:49)) )
Summarym<-melt(Summary)
ggplot(Summarym, aes(x=probe, y=value, shape=variable,color=ethnicity))+
    geom_point(position=position_dodge(width=0.3),size =1)+
    geom_path(aes(x=probe, y=value,group=interaction(ethnicity,variable)),position=position_dodge(width=0.3),linetype = 1, size=0.2)+
    scale_color_manual(values=colors[c(1,3,5)], name="Ethnicity")+
    scale_shape_manual(values=c(19,17), name="",labels=c("Deletion", "Duplication"))+
    ylab("Observed frequency")+
    xlab("")+
    theme_bw()+
    theme(axis.title.y = element_text(size=12))+
    theme(panel.grid.minor.y = element_blank() )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    theme(axis.text.x = element_text(angle = 90,size=8 ))
ggsave("Output/CNVs_plot.pdf", width = 11, height = 3.5)



########### 
#Breakpoints plots
#SumDF<-read.csv("Output/Breakpoints_Summary.csv",row.names = 1)
#SumDF$CopyNo<-as.character(SumDF$CopyNo)

bk5<-read.csv("Output/bk5.csv", row.names = 1)
bk5.2<-read.csv("Output/bk5.2.csv", row.names=1)
bk5$Var1<-factor(bk5$Var1, levels=0:49)
bk5.2$Var1<-factor(bk5.2$Var1, levels=0:49)

ggplot(bk5, aes(x=Var1, y=freq, fill=ethnicity, group=ethnicity))+
    geom_bar(stat="identity", position=position_dodge(), width=0.8)+
    scale_fill_manual(values=paste0(colors[c(1,3,5)]), name="Ethnicity", labels=c("A","B","C"))+
    ylab("Observed frequency (5' breakpoint)")+
    xlab("")+
    theme_bw()+
    theme(axis.text.x = element_text(size=9),axis.title.y = element_text(size=12))+
    geom_vline(xintercept = c(1:46)+0.5,  
               color = "gray80", size=.2)+
    theme(panel.grid.major.x = element_blank())+
    theme(legend.title = element_blank())
ggsave("Output/breakpoint-5.byProbe.pdf", width = 10, height = 4)


ggplot(bk5.2, aes(x=Var1, y=freq, shape=type,color=ethnicity))+
    geom_point(position=position_dodge(width=0.3),size =1)+
    geom_path(aes(x=Var1, y=freq,group=interaction(ethnicity,type)),position=position_dodge(width=0.3),linetype = 1, size=0.2)+
    scale_color_manual(values=colors[c(1,3,5)], name="Ethnicity")+
    scale_shape_manual(values=c(19,17), name="",labels=c("Deletion", "Duplication"))+
    ylab("Observed frequency")+
    xlab("CNSL Probe")+
    theme_bw()+
    theme(axis.title.y = element_text(size=12))+
    theme(panel.grid.minor.y = element_blank() )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    theme(axis.text.x = element_text(size=9 ))
ggsave("Output/Breakpoint5_dup.del.separate_pointplot.pdf", width = 9, height = 4)


bk3<-read.csv("Output/bk3.csv", row.names = 1)
bk3.2<-read.csv("Output/bk3.2.csv", row.names=1)
bk3$Var1<-factor(bk3$Var1, levels=0:49)
bk3.2$Var1<-factor(bk3.2$Var1, levels=0:49)

ggplot(bk3, aes(x=Var1, y=freq, fill=ethnicity, group=ethnicity))+
    geom_bar(stat="identity", position=position_dodge(), width=0.8)+
    scale_fill_manual(values=paste0(colors[c(1,3,5)]), name="Ethnicity", labels=c("A","B","C"))+
    ylab("Observed frequency (3' breakpoint)")+
    xlab("CNSL Probe")+
    theme_bw()+
    theme(axis.text.x = element_text(size=10),axis.title.y = element_text(size=12),axis.title.x = element_text(size=12))+
    geom_vline(xintercept = c(1:46)+0.5,  
               color = "gray80", size=.2)+
    theme(panel.grid.major.x = element_blank())+
    theme(legend.title = element_blank())
ggsave("Output/breakpoint-3.byProbe2.pdf", width = 10, height = 4)


ggplot(bk3.2, aes(x=Var1, y=freq, shape=type,color=ethnicity))+
    geom_point(position=position_dodge(width=0.3),size =1)+
    geom_path(aes(x=Var1, y=freq,group=interaction(ethnicity,type)),position=position_dodge(width=0.3),linetype = 1, size=0.2)+
    scale_color_manual(values=colors[c(1,3,5)], name="Ethnicity")+
    scale_shape_manual(values=c(19,17), name="",labels=c("Deletion", "Duplication"))+
    ylab("Observed frequency")+
    xlab("CNSL Probe")+
    theme_bw()+
    theme(axis.title.y = element_text(size=12))+
    theme(panel.grid.minor.y = element_blank() )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    theme(axis.text.x = element_text(size=9 ))
ggsave("Output/Breakpoint3_dup.del.separate_pointplot.pdf", width = 9, height = 4)




#Devide between copy nubmer 1 and 3
SumDe<-SumDF[SumDF$CopyNo=="1",]    
SumDu<-SumDF[SumDF$CopyNo=="3",]    

ggplot(SumDe, aes(y=Sample, color=ethnicity))+
    geom_linerangeh(aes(xmin=min,xmax=max))+
    scale_color_manual(values=colors[c(1,3,5)], name="Ethnicity",labels=c("A", "B", "C"))+
    ylab("")+
    theme_bw()
ggsave("Output/Deletion.matrix.pdf", width = 8,height = 8)
ggplot(SumDu, aes(y=Sample, color=ethnicity))+
    geom_linerangeh(aes(xmin=min,xmax=max))+
    scale_color_manual(values=colors[c(1,3,5)], name="Ethnicity",labels=c("A", "B", "C"))+
    ylab("")+
    theme_bw()
ggsave("Output/Duplication.matrix.pdf", width = 8,height = 8)


# plot by sorting by startpoint/length
SumDu2<-SumDu %>% 
    arrange_at("min") %>%
    arrange_at("length", desc)


SumDu2$number<-""
SumDu2$number[1]<-1
k=2
for (i in 2: nrow(SumDu2)){
    if (SumDu2$Sample[i] %in% SumDu2$Sample[1:(i-1)]) SumDu2$number[i]<-SumDu2$number[which(SumDu2$Sample==SumDu2$Sample[i])[1]]
    else {
        SumDu2$number[i]<-k
        k<-k+1}
}
SumDu2$number<-as.integer(SumDu2$number)
ggplot(SumDu2, aes(y=number, color=ethnicity))+
    geom_linerangeh(aes(xmin=min,xmax=max))+
    scale_color_manual(values=colors[c(1,3,5)], name="Ethnicity",labels=c("A", "B", "C"))+
    ylab("")+
    theme_bw()  
ggsave("Output/Duplication.matrix_ordered.pdf", width = 8,height = 8)


#plot by ethnicity
ethnic<-sort(unique(copyNo$ethnicity)) 
col.vec<-c(1,3,5)
for (i in 1:3){
    dup<-SumDu2[SumDu2$ethnicity==ethnic[i],]
    ggplot(dup, aes(y=number, color=ethnicity))+
        geom_linerangeh(aes(xmin=min,xmax=max))+
        scale_color_manual(values=colors[col.vec[i]])+
        ylab("")+
        theme_bw()  
    ggsave(paste0("Output/Duplication.matrix_ordered.", ethnic[i],".pdf"), width = 8,height = 8)
}

#deletion
SumDe2<-SumDe %>% arrange_at("min") %>% arrange_at("length", desc)
SumDe2$number<-""
SumDe2$number[1]<-1
k=2
for (i in 2: nrow(SumDe2)){
    if (SumDe2$Sample[i] %in% SumDe2$Sample[1:(i-1)]) SumDe2$number[i]<-SumDe2$number[which(SumDe2$Sample==SumDe2$Sample[i])[1]]
    else {
        SumDe2$number[i]<-k
        k<-k+1}
}
SumDe2$number<-as.integer(SumDe2$number)
ggplot(SumDe2, aes(y=number, color=ethnicity))+
    geom_linerangeh(aes(xmin=min,xmax=max))+
    scale_color_manual(values=colors[c(1,3,5)], name="Ethnicity",labels=c("A", "B", "C"))+
    ylab("")+
    theme_bw()  
ggsave("Output/Deletion.matrix_ordered.pdf", width = 8,height = 8)


#plot by ethnicity
ethnic<-sort(unique(copyNo$ethnicity)) 
col.vec<-c(1,3,5)
for (i in 1:3){
    dup<-SumDe2[SumDe2$ethnicity==ethnic[i],]
    ggplot(dup, aes(y=number, color=ethnicity))+
        geom_linerangeh(aes(xmin=min,xmax=max))+
        scale_color_manual(values=colors[col.vec[i]])+
        ylab("")+
        theme_bw()  
    ggsave(paste0("Output/Deletion.matrix_ordered.", ethnic[i],".pdf"), width = 8,height = 8)
}

