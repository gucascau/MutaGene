

### draw some examples to show the distribution of SNV/Indel/CNV burden for overlapping genes


#### The script is for figure 3D

library("ggplot2")

library(reshape2)
library(dplyr)

getwd()
setwd("/Users/xinwang/Documents/Projects/Proposal/Figure3/MainFigureV1020_20220325V/Overlapping_Mutators_20200623/OverlappingMutator_forfigures/Example_intergrated/")

#myfiles<-list.files(pattern="*associatedgenes.txt"

#### show the low number of mutations

Mutationload_gene<-read.table("Mutationload_DecMut/CNVdecreased_Omutators.forMLF2.txt",header = T,sep = "\t", na.strings = NA, fill = T)

attach(Mutationload_gene)
levels(Type)

Figure3DCNVdata<-cbind(Mutationload_gene,log10(Mutationload_gene$Mutationload))
write.table(as.matrix(Figure3DCNVdata),file="Fig3DPart2_datasouce.xls",row.names = TRUE, sep="\t")


ggplot(Mutationload_gene,aes(x=factor(Mutationload_gene$Gene),y=log10(Mutationload_gene$Mutationload),fill=Type))+
 # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_boxplot(outlier.shape=NA,width=0.5)+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
  ) +ylim(1,4)




#### For list2 increased genes

#setwd("/Users/xinwang/Documents/Project/Proposal/Figure3/Overlapping_Mutators_20200623/OverlappingMutator_forfigures/Example_intergrated/Mutationload_IncMut/")
Mutationload_gene<-read.table("Mutationload_IncMut/CNVInc_Omutators.forMLF2.txt",header = T,sep = "\t", na.strings = NA, fill = T)

attach(Mutationload_gene)
levels(Mutationload_gene$Type)

Figure3D2CNVdata<-cbind(Mutationload_gene,log10(Mutationload_gene$Mutationload))
write.table(as.matrix(Figure3DCNVdata),file="Fig3DPart1_datasouce.xls",row.names = TRUE, sep="\t")

ggplot(Mutationload_gene,aes(x=factor(Mutationload_gene$Gene),y=log10(Mutationload_gene$Mutationload),fill=Type))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_boxplot(outlier.shape=NA,width=0.5)+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
        ) +ylim(1,4)




##### for lists indel gene increased

#setwd("/Users/xinwang/Documents/Project/Proposal/Figure3/Overlapping_Mutators_20200623/OverlappingMutator_forfigures/Example_intergrated/Indelload_IncMut")

Mutationload_gene<-read.table("Indelload_IncMut/CNVInc_IndelloadExam2.txt",header = T,sep = "\t", na.strings = NA, fill = T)

attach(Mutationload_gene)
levels(Mutationload_gene$Type)


Figure3D3CNVdata<-cbind(Mutationload_gene,log10(Mutationload_gene$Indelload))
write.table(as.matrix(Figure3D3CNVdata),file="Fig3DPart3_datasouce.xls",row.names = TRUE, sep="\t")


ggplot(Mutationload_gene,aes(x=factor(Gene),y=log10(Indelload),fill=Type))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_boxplot(outlier.shape=NA,width=0.5)+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
  ) 


##### for lists indel gene decreased

#setwd("/Users/xinwang/Documents/Project/Proposal/Figure3/Overlapping_Mutators_20200623/OverlappingMutator_forfigures/Example_intergrated/Indelload_DecMut")

  Mutationload_gene<-read.table("Indelload_DecMut/CNVdec_Indelloadexam2.txt",header = T,sep = "\t", na.strings = NA, fill = T)

attach(Mutationload_gene)
levels(Mutationload_gene$Type)

Figure3D4CNVdata<-cbind(Mutationload_gene,log10(Mutationload_gene$Indelload))
write.table(as.matrix(Figure3D4CNVdata),file="Fig3DPart4_datasouce.xls",row.names = TRUE, sep="\t")



ggplot(Mutationload_gene,aes(x=factor(Gene),y=log10(Indelload),fill=Type))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_boxplot(outlier.shape=NA,width=0.5)+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
  ) 



##### for lists CNV gene decreased

#setwd("/Users/xinwang/Documents/Project/Proposal/Figure3/Overlapping_Mutators_20200623/OverlappingMutator_forfigures/Example_intergrated/CNVload_DecMut/")

  Mutationload_gene<-read.table("CNVload_DecMut/CNVload_decMutExam2.txt",header = T,sep = "\t", na.strings = NA, fill = T)

attach(Mutationload_gene)
levels(Mutationload_gene$Type)

Figure3D5CNVdata<-cbind(Mutationload_gene,log10(Mutationload_gene$CNV))
write.table(as.matrix(Figure3D5CNVdata),file="Fig3DPart5_datasouce.xls",row.names = TRUE, sep="\t")



ggplot(Mutationload_gene,aes(x=factor(Gene),y=log10(CNV),fill=Type))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_boxplot(outlier.shape=NA,width=0.5)+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
  ) 



##### for lists CNV gene increase

#setwd("/Users/xinwang/Documents/Project/Proposal/Figure3/Overlapping_Mutators_20200623/OverlappingMutator_forfigures/Example_intergrated/CNVload_IncMut/")

Mutationload_gene<-read.table("CNVload_IncMut/CNVload_increaseMut2.txt",header = T,sep = "\t", na.strings = NA, fill = T)

attach(Mutationload_gene)
levels(Mutationload_gene$Type)

Figure3D6CNVdata<-cbind(Mutationload_gene,log10(Mutationload_gene$CNV))
write.table(as.matrix(Figure3D6CNVdata),file="Fig3DPart6_datasouce.xls",row.names = TRUE, sep="\t")



ggplot(Mutationload_gene,aes(x=factor(Gene),y=log10(CNV),fill=Type))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_boxplot(outlier.shape=NA,width=0.5)+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8),legend.position="bottom"
  ) 
