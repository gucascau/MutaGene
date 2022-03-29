
### For list of Pathway mutations and mutation load indel load and CNV loads

library("ggplot2")

library(reshape2)
library(dplyr)
setwd("/Users/xinwang/Documents/Projects/Proposal/Figure4/V20211205/")

## add stringAsFactors true


Mutationload_gene<-read.table("ALLSample_Divided_intoCate_neanti_withDDRNM.txt",header = T,sep = "\t", na.strings = NA, fill = T,stringsAsFactors=T)

#Mutationload_gene<-Mutationload_gene2[grep("NonDDR",Mutationload_gene2$Pathway,invert = TRUE),]
attach(Mutationload_gene)



levels(Mutationload_gene$Pathway)
Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(4,5,10,11,1,9,12,8,15,16,2,3,6,7,14,13)],ordered = TRUE)

levels(Mutationload_gene$Pathway)

attach(Mutationload_gene)
head(Mutationload_gene)



Figure4B2<-cbind(Pathway,log10(Mutationload))
write.table(as.matrix(Figure4B2),file="Fig4BPartSNV2_datasouce.xls",row.names = TRUE, sep="\t")


### For Mutation load
ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(Mutationload), fill=Pathway))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  # geom_violin(fill="cyan2")+ 
  #geom_boxplot(width=0.1, fill="cyan2",coef=1e30)+
  ### indianred1, orchid1, cornflowerblue
  #stat_boxplot(geom ='errorbar')+
  geom_boxplot(outlier.shape=NA)+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
        
        ##### color: lightgoldenrod1  ML:mediumaquamarine, Indel:orchid1, CNV:lightskyblue
  )  +scale_fill_manual(values=c("brown1","wheat",
                                 "mediumaquamarine","mediumaquamarine","mediumaquamarine",
                                 "mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine",
                                 "mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine",
                                 "gray")) +
  ylim(1,5)


#### calculate the Pvalue for mutationload 4, indel load 5, CNV 6, Neoantigene load 11, BindingExpressed PMHC 13

for (i in (levels(Mutationload_gene$Pathway))){
  
  mutator<-Mutationload_gene[Mutationload_gene$Pathway==i,4]
  nonmutant<-Mutationload_gene[Mutationload_gene$Pathway=="NonDDR",4]
  pvalue<-wilcox.test(mutator, nonmutant, alternative = "greater")$p.value
  pvalue
  print(paste(i,pvalue))
}


### For indel load
ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(Indel), fill=Pathway))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  # geom_violin(fill="cyan2")+ 
  #geom_boxplot(width=0.1, fill="cyan2",coef=1e30)+
  ### indianred1, orchid1, cornflowerblue
  #stat_boxplot(geom ='errorbar')+
  geom_boxplot(outlier.shape=NA)+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
        
        ##### color: lightgoldenrod1  ML:mediumaquamarine, Indel:orchid1, CNV:lightskyblue
  )  +scale_fill_manual(values=c("brown1","wheat",
                                 "orchid1","orchid1","orchid1",
                                 "orchid1","orchid1","orchid1","orchid1","orchid1",
                                 "orchid1","orchid1","orchid1","orchid1","orchid1",
                                "gray")) +
  ylim(0,4)



#### calculate the Pvalue for mutationload 4, indel load 5, CNV 6, Neoantigene load 11, BindingExpressed PMHC 13

for (i in (levels(Mutationload_gene$Pathway))){
  
  mutator<-Mutationload_gene[Mutationload_gene$Pathway==i,6]
  nonmutant<-Mutationload_gene[Mutationload_gene$Pathway=="NonDDR",6]
  pvalue<-wilcox.test(mutator, nonmutant, alternative = "greater")$p.value
  pvalue
  print(paste(i,pvalue))
}

### For CNV load
ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(CNV), fill=Pathway))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  # geom_violin(fill="cyan2")+ 
  #geom_boxplot(width=0.1, fill="cyan2",coef=1e30)+
  ### indianred1, orchid1, cornflowerblue
  #stat_boxplot(geom ='errorbar')+
  geom_boxplot(outlier.shape=NA)+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
        
        ##### color: lightgoldenrod1  ML:mediumaquamarine, Indel:orchid1, CNV:lightskyblue
  )  +scale_fill_manual(values=c("brown1","wheat",
                                 "lightskyblue","lightskyblue","lightskyblue",
                                 "lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue",
                                 "lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue",
                                 "gray")) +
  ylim(1.5,4)



#### calculate the Pvalue for mutationload 4, indel load 5, CNV 6, Neoantigene load 11, BindingExpressed PMHC 13

for (i in (levels(Mutationload_gene$Pathway))){
  
  mutator<-Mutationload_gene[Mutationload_gene$Pathway==i,6]
  nonmutant<-Mutationload_gene[Mutationload_gene$Pathway=="NonDDR",6]
  pvalue<-wilcox.test(mutator, nonmutant, alternative = "greater")$p.value
  pvalue
  print(paste(i,pvalue))
}




#### Neoantigen load
ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(numberOfImmunogenicMutation), fill=Pathway))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_boxplot(outlier.shape=NA)+
  ### cyan2 and orchid1
  # geom_violin()+ 
   # geom_boxplot(width=0.1, coef=1e30)+
  ### indianred1, orchid1, cornflowerblue
  #geom_boxplot(outlier.shape=NA,fill="orchid1")+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
        ###mediumaquamarine and cornflowerblue
  ) +scale_fill_manual(values=c("brown1","wheat",
                                "chocolate","chocolate","chocolate","chocolate","chocolate","chocolate","chocolate","chocolate","chocolate","chocolate","chocolate","chocolate","chocolate",
                                "darkgray","gray")) +
  ylim(0,5)


for (i in (levels(Mutationload_gene$Pathway))){
  
  mutator<-Mutationload_gene[Mutationload_gene$Pathway==i,11]
  nonmutant<-Mutationload_gene[Mutationload_gene$Pathway=="DDRM",11]
  pvalue<-wilcox.test(mutator, nonmutant, alternative = "less")$p.value
  pvalue
  print(paste(i,pvalue))
}




















Mutationload_gene<-read.table("ALLSample_Divided_intoCate.txt",header = T,sep = "\t", na.strings = NA, fill = T,stringsAsFactors=T)

attach(Mutationload_gene)


factor(Mutationload_gene$Pathway)
levels(Mutationload_gene$Pathway)
Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(4,9,10,1,8,11,7,14,15,2,3,5,6,13,12)],ordered = TRUE)

levels(Mutationload_gene$Pathway)

attach(Mutationload_gene)

ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(Mutationload)))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
 # geom_violin(fill="lightblue")+ 
 # geom_boxplot(width=0.1, fill="lightblue",coef=1e30)+
  ### indianred1, orchid1, cornflowerblue
  geom_boxplot(outlier.shape=NA,fill="indianred1")+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
  ) +ylim(1,5)

#### calculation the pvalue of each pathway compared with nonDDR

tail(Mutationload_gene)

##CCNV ####
NonmutCNV<-Mutationload_gene[Mutationload_gene$Pathway=="NER",6]
MutPCNV<-Mutationload_gene[Mutationload_gene$Pathway=="DDRM",6]
MutPCNV 

wilcox.test(MutPCNV, NonmutCNV, alternative = "greater")

##ML ####
NonmutML<-Mutationload_gene[Mutationload_gene$Pathway=="Othersusp",4]
MutPML<-Mutationload_gene[Mutationload_gene$Pathway=="DDRM",4]
NonmutML
MutPML 

wilcox.test(MutPML, NonmutML, alternative = "greater")


##Indel ####
NonmutIndel<-Mutationload_gene[Mutationload_gene$Pathway=="Othersusp",5]
MutPIndel<-Mutationload_gene[Mutationload_gene$Pathway=="DDRM",5]
MutPIndel 

wilcox.test(MutPIndel, NonmutIndel, alternative = "greater")




#### neoantigen ####

setwd("/Users/xinwang/Documents/Project/Proposal/Figure3/Version_20200511/Final_All/BRCA_figure/DDR_deficiency/ImmuneDDR/")

Mutationload_gene<-read.table("ALLSample_Divided_intoCate_neanti.txt",header = T,sep = "\t", na.strings = NA, fill = T)

attach(Mutationload_gene)


levels(Pathway)
Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(4,9,10,1,8,11,7,14,15,2,3,5,6,13,12)],ordered = TRUE)

levels(Mutationload_gene$Pathway)

attach(Mutationload_gene)
head(Mutationload_gene)
ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(numberOfImmunogenicMutation)))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_violin(fill="cyan2")+ 
  geom_boxplot(width=0.1, fill="cyan2",coef=1e30)+
  ### indianred1, orchid1, cornflowerblue
  #geom_boxplot(outlier.shape=NA,fill="coral")+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
  ) +ylim(0,5)


ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(numberOfImmunogenicMutation)))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_violin(fill="cyan2")+ 
  geom_boxplot(width=0.1, fill="cyan2",coef=1e30)+
  ### indianred1, orchid1, cornflowerblue
  #geom_boxplot(outlier.shape=NA,fill="coral")+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
  ) +ylim(0,4)

####### Then neoantigen P value calculation

### number of ImmuogenicMutations
head(Mutationload_gene)
NonmutML<-Mutationload_gene[Mutationload_gene$Pathway=="EditNucl",11]
MutPML<-Mutationload_gene[Mutationload_gene$Pathway=="DDRM",11]
NonmutML
MutPML 

wilcox.test(MutPML, NonmutML, alternative = "greater")

### number of numberOfBindingExpressedPMHC
head(Mutationload_gene)
NonmutML<-Mutationload_gene[Mutationload_gene$Pathway=="Othersusp",14]
MutPML<-Mutationload_gene[Mutationload_gene$Pathway=="DDRM",14]
NonmutML
MutPML 


wilcox.test(MutPML, NonmutML, alternative = "greater")



#### Only calculated the DDRM and nonDDR

setwd("/Users/xinwang/Documents/Project/Proposal/Figure3/Version_20200511/Final_All/BRCA_figure/DDR_deficiency/ImmuneDDR/")

Mutationload_gene<-read.table("ALLSample_Divided_intoCate_neanti_DDR_DDRM_nonDDR.txt",header = T,sep = "\t", na.strings = NA, fill = T)


attach(Mutationload_gene)
head(Mutationload_gene)

Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(2,1,3)],ordered = TRUE)
Mutationload_gene$Pathway

ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(numberOfBindingExpressedPMHC), fill=factor(Pathway)))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_violin()+ 
  geom_boxplot(width=0.1,coef=1e30, fill="gray")+
  ### indianred1, orchid1, cornflowerblue
  #geom_boxplot(outlier.shape=NA,fill="coral")+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
  ) +scale_fill_brewer(palette="Dark2") +ylim(0,4)



mutant<-Mutationload_gene[Mutationload_gene$Pathway=="DDR",13]
mutant
nonmutant<-Mutationload_gene[Mutationload_gene$Pathway=="NonDDR",13]
mutator<-Mutationload_gene[Mutationload_gene$Pathway=="DDRM",13]

wilcox.test(mutant, nonmutant, alternative = "greater")




ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(numberOfImmunogenicMutation), fill=factor(Pathway)))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_violin()+ 
  geom_boxplot(width=0.1,coef=1e30, fill="gray")+
  ### indianred1, orchid1, cornflowerblue
  #geom_boxplot(outlier.shape=NA,fill="coral")+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
  ) +scale_fill_brewer(palette="Dark2") +ylim(0,4)


### caculate the pvalue for different groups


+scale_fill_manual(values=c("red", "blue","gray" ))


+scale_x_discrete(limits=c("DDRM", "NA", "NonDDR"))+ylim(0,4)



###


####################################################################
###################Update version after discussion #################
####################################################################

############# Here is to compare DDRM and DDRNM

#### Only calculated the DDRM and DDRNM

setwd("/Users/xinwang/Documents/Project/Proposal/Figure3/Version_20200511/Final_All/BRCA_figure/DDR_deficiency/Version_July/")
Mutationload_gene<-read.table("DDRM_DDRNM_Categories.txt",header = T,sep = "\t", na.strings = NA, fill = T)


attach(Mutationload_gene)
head(Mutationload_gene)

Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(1,2)],ordered = TRUE)
Mutationload_gene$Pathway

  ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(numberOfBindingExpressedPMHC), fill=factor(Pathway)))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  geom_violin()+ 
  geom_boxplot(width=0.1,coef=1e30, fill="gray")+
  ### indianred1, orchid1, cornflowerblue
  #geom_boxplot(outlier.shape=NA,fill="coral")+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
  ) +scale_fill_brewer(palette="Dark2") +ylim(0,4)



### Compare between DDRM and DDRNM


DDRNM_Binding<-Mutationload_gene[Mutationload_gene$Pathway=="DDRNM",11]
DDRM_Binding<-Mutationload_gene[Mutationload_gene$Pathway=="DDRM",11]

wilcox.test(DDRM_Binding, DDRNM_Binding, alternative = "greater")


#### Combined with DDRNM, DDRM, DDR pathway and NonDDR pathway
dev.off()

Mutationload_gene<-read.table("ALLSample_Divided_intoCate_neanti_withDDRNM.txt",header = T,sep = "\t", na.strings = NA, fill = T)

#Mutationload_gene<-Mutationload_gene2[grep("NonDDR",Mutationload_gene2$Pathway,invert = TRUE),]
attach(Mutationload_gene)



levels(Mutationload_gene$Pathway)
Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(5,6,11,12,1,10,13,9,16,17,2,3,7,8,15,4,14)],ordered = TRUE)

levels(Mutationload_gene$Pathway)

attach(Mutationload_gene)
head(Mutationload_gene)
ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(Indel), fill=Pathway))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
 # geom_violin(fill="cyan2")+ 
  #geom_boxplot(width=0.1, fill="cyan2",coef=1e30)+
  ### indianred1, orchid1, cornflowerblue
  #stat_boxplot(geom ='errorbar')+
  geom_boxplot(outlier.shape=NA)+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
        
        ##### color: lightgoldenrod1  ML:mediumaquamarine, Indel:orchid1, CNV:lightskyblue
  )  +scale_fill_manual(values=c("brown1","wheat",
                                 "orchid1","orchid1","orchid1",
                                 "orchid1","orchid1","orchid1","orchid1","orchid1",
                                 "orchid1","orchid1","orchid1","orchid1","orchid1",
                                 "darkgray","gray")) +
  ylim(0,4)


#### calculate the Pvalue for mutationload 4, indel load 5, CNV 6, Neoantigene load 11, BindingExpressed PMHC 13

for (i in (levels(Mutationload_gene$Pathway))){
  
  mutator<-Mutationload_gene[Mutationload_gene$Pathway==i,6]
  nonmutant<-Mutationload_gene[Mutationload_gene$Pathway=="NonDDR",6]
  pvalue<-wilcox.test(mutator, nonmutant, alternative = "greater")$p.value
pvalue
  print(paste(i,pvalue))
}


#### Neoantigen load


ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(numberOfBindingPMHC), fill=Pathway))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  ### cyan2 and orchid1
  geom_violin()+ 
  geom_boxplot(width=0.1, coef=1e30)+
  ### indianred1, orchid1, cornflowerblue
  #geom_boxplot(outlier.shape=NA,fill="orchid1")+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
        ###mediumaquamarine and cornflowerblue
  ) +scale_fill_manual(values=c("brown1","wheat",
                                "mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine","mediumaquamarine",
                                "gray")) +
  ylim(0,6)

## expressed neoantige

ggplot(Mutationload_gene,aes(x=factor(Pathway),y=log10(numberOfBindingExpressedPMHC), fill=Pathway))+
  # geom_violin()+ 
  #geom_boxplot(notch = TRUE,width=0.5)+
  # geom_dotplot(binaxis ='y', stackdir='center',stackratio=1, dotsize=0.5)+
  
  # geom_violin(fill="gray")+ 
  ### cyan2 and orchid1
   geom_violin()+ 
  geom_boxplot(width=0.1, coef=1e30)+
  ### indianred1, orchid1, cornflowerblue
  #geom_boxplot(outlier.shape=NA,fill="orchid1")+
  #stat_compare_means(comparisons = my_comparisons, label.y = c(29, 35, 40))+     # Add global p-value
  # stat_compare_means(label.y = 45) +
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size=8)
        ###mediumaquamarine and cornflowerblue
  ) +scale_fill_manual(values=c("brown1","wheat",
                                "lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue","lightskyblue",
                                "darkgray","gray")) +
  ylim(0,6)


