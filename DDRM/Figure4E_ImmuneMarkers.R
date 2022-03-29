
##### This script is to draw the volcano plot to compare the significant different expression genes between DDRM and nonDDR cancer types

setwd("/Users/xinwang/Documents/Projects/Proposal/Figure4/DDR_deficiency/ImmuneDDR/")
### for Amplication

library(ggplot2)
library(ggrepel)
Expression_var<-read.table("DDRM_Diff_geneExp_viocano2.txt", header = T,sep = "\t", na.strings = NA, fill = T)
attach(Expression_var)
ggplot(Expression_var, aes(x=log2(Difference_2), y=Padjust)) +
  geom_point(size=1.2,color= ifelse(Difference_2>2.5 & Padjust>=2,"red",
                                    ifelse(Difference_2<=0.4 & Padjust>=2,"blue","gray")))+
  geom_text_repel(aes(label=ifelse((Difference_2>=5| Difference_2<= 0.1 ) & Padjust>=3,as.character(Gene),'')),hjust=0.5, 
                  size=2,segment.size=0.02)+
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )+
  theme(legend.justification = c(1, 1), 
        legend.text= element_text(face="italic"),legend.position = c(1,1))

log2(1/5)

+
  scale_x_continuous(limits=c(-6,6),breaks =c(-6,-4,-2,0,2,4,6))



#### Draw the immune therapy genes expression
library(reshape2)
ImmuneGeneExp<-read.table("DDRM_ImmuneGeneExp.txt", sep = "\t",header = T)
attach(ImmuneGeneExp)
ImmuneGeneExpTranfer<-melt(ImmuneGeneExp)
head(ImmuneGeneExpTranfer)

levels(ImmuneGeneExpTranfer$variable)
ImmuneGeneExpTranfer$variable=factor(ImmuneGeneExpTranfer$variable,levels = levels(ImmuneGeneExpTranfer$variable)[c(5,1,2,3,4)],ordered = TRUE)

ggplot(ImmuneGeneExpTranfer, aes(x=variable, y=value,fill=Type)) +facet_wrap(.~variable, scale="free",nrow = 1)+ geom_boxplot(outlier.shape=NA,position=position_dodge(0.8)) + theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )+ theme(legend.position = "bottom")


dev.off()
library(reshape2)
ImmuneGeneExp<-read.table("DDRM_DDR_NonDDR_Immune_genes.txt", sep = "\t",header = T)
attach(ImmuneGeneExp)
ImmuneGeneExpTranfer<-melt(ImmuneGeneExp)
head(ImmuneGeneExpTranfer)

levels(ImmuneGeneExpTranfer$variable)
levels(ImmuneGeneExpTranfer$Type)
ImmuneGeneExpTranfer$variable=factor(ImmuneGeneExpTranfer$variable,levels = levels(ImmuneGeneExpTranfer$variable)[c(5,1,2,3,4)],ordered = TRUE)
### only draw one varaib
ImmuneGeneExpTranfer$Type=factor(ImmuneGeneExpTranfer$Type,levels = levels(ImmuneGeneExpTranfer$Type)[c(2,1,3)],ordered = TRUE)
ggplot(ImmuneGeneExpTranfer, aes(x=variable, y=log2(value),fill=Type)) +facet_wrap(.~variable, scale="free",nrow = 1)+ geom_boxplot(outlier.shape=NA,position=position_dodge(0.8)) + theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  #)+ theme(legend.position = "bottom") +scale_fill_manual(values=c("red","green","blue")) +ylim(0,12.5)

        
  )+ theme(legend.position = "bottom") +scale_fill_brewer(palette="Dark2") +ylim(0,12.5)



tail(ImmuneGeneExpTranfer)

PDL1<-ImmuneGeneExpTranfer[ImmuneGeneExpTranfer$variable=="CD8A",]
tail(PDL1)

levels(PDL1$Type)
mutator<-PDL1[PDL1$Type=="DDRM",3]
mutant<-PDL1[PDL1$Type=="DDR",3]
nonmutant<-PDL1[PDL1$Type=="NonDDR",3]

wilcox.test(mutator, nonmutant, alternative = "greater")
wilcox.test(mutant, nonmutant, alternative = "greater")
wilcox.test(mutator, mutant, alternative = "greater")




#### Pathway different Exrpression of five immunetherapy genes
#### for genes CD274 CD8A, CTLA4, PDCD1, and PDCD1LG2

setwd("/Users/xinwang/Documents/Projects/Proposal/Figure4/DDR_deficiency/ImmuneDDR/")

Mutationload_gene<-read.table("ALLSample_Divided_intoCate_immunegeneExp.txt",header = T,sep = "\t", na.strings = NA, fill = T)
attach(Mutationload_gene)

GeneExpressTr<- melt(Mutationload_gene,id.vars = c("Pathway2"),measure.vars = c("PDL1","CD8A","PD1","PDL2","CTLA4","CD8B"))
GeneExpressTr

attach(GeneExpressTr)
levels(GeneExpressTr$variable)
head(GeneExpressTr)
GeneExpressTr$variable<- factor(GeneExpressTr$variable, levels=levels(GeneExpressTr$variable)[c(3,1,4,2,6,5)],ordered=TRUE)

GeneExpressTr$Pathway2=factor(GeneExpressTr$Pathway2,levels = levels(GeneExpressTr$Pathway2)[c(2,1,3)],ordered = TRUE)

#Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(2,1,3)],ordered = TRUE)
#Mutationload_gene$Pathway
#Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(4,9,10,1,8,11,7,14,15,2,3,5,6,13,12)],ordered = TRUE)

ggplot(GeneExpressTr, aes(x=variable, y=log2(value),fill=Pathway2)) +facet_wrap(.~variable, scale="free",nrow = 1)+ geom_boxplot(outlier.shape=NA,position=position_dodge(0.8)) + theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        #)+ theme(legend.position = "bottom") +scale_fill_manual(values=c("red","green","blue")) +ylim(0,12.5)
        
        
  )+ theme(legend.position = "bottom") +scale_fill_brewer(palette="Dark2") +ylim(0,14.5)



PDL1<-GeneExpressTr[GeneExpressTr$variable=="CTLA4",]
tail(PDL1)

levels(PDL1$Pathway2)
mutator<-PDL1[PDL1$Pathway2=="DDRM",3]
mutant<-PDL1[PDL1$Pathway2=="DDR",3]
nonmutant<-PDL1[PDL1$Pathway2=="NonDDR",3]


wilcox.test(mutator, nonmutant, alternative = "greater")
wilcox.test(mutant, nonmutant, alternative = "greater")
wilcox.test(mutator, mutant, alternative = "greater")



#### Draw the different RSEM across diffferent DDR pathway

GeneExpressTr2<- melt(Mutationload_gene,id.vars = c("Pathway"),measure.vars = c("PDL1","CD8A","PD1","PDL2","CTLA4","CD8B"))
GeneExpressTr2

attach(GeneExpressTr2)

head(GeneExpressTr2)
levels(GeneExpressTr2$variable)
GeneExpressTr2$variable<- factor(GeneExpressTr2$variable, levels=levels(GeneExpressTr2$variable)[c(3,1,4,2,6,5)],ordered=TRUE)
GeneExpressTr2$Pathway<-factor(GeneExpressTr2$Pathway,levels = levels(GeneExpressTr2$Pathway)[c(12,13,6,5,3,2,15,14,7,11,8,1,10,9,4)],ordered = TRUE)
  #4,9,10,1,8,11,7,14,15,2,3,5,6,13,12)],ordered = TRUE)




ggplot(GeneExpressTr2, aes(x=Pathway, y=log2(value),fill=variable)) +facet_wrap(.~variable, scale="free",nrow = 1)+ 
  geom_violin()+ 
  geom_boxplot(width=0.1,coef=1e30, fill="gray")+
  #geom_boxplot(outlier.shape=NA,position=position_dodge(0.8)) + 
  
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=10),
        axis.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        #)+ theme(legend.position = "bottom") +scale_fill_manual(values=c("red","green","blue")) +ylim(0,12.5)
        
        
  )+ coord_flip()+theme(legend.position = "bottom") +scale_fill_brewer(palette="Dark2") +ylim(0,14.5)


### calculate the Pvalue of each genes


PDL1<-GeneExpressTr2[GeneExpressTr2$variable=="CTLA4",]
tail(PDL1)

levels(PDL1$Pathway)

for (i in (rev(levels(PDL1$Pathway)))){
  
  mutator<-PDL1[PDL1$Pathway=="DDRM",3]
  nonmutant<-PDL1[PDL1$Pathway==i,3]
  pvalue<-wilcox.test(mutator, nonmutant, alternative = "greater")$p.value
  pvalue
  print(paste(i,pvalue))
}
mutator<-PDL1[PDL1$Pathway==i,3]
#mutant<-PDL1[PDL1$Pathway2=="DDR",3]
nonmutant<-PDL1[PDL1$Pathway=="NonDDR",3]


wilcox.test(mutator, nonmutant, alternative = "greater")



#############################
#### Update version #########
#############################

##### This script is to draw the volcano plot to compare the significant different expression genes between DDRM and nonDDR cancer types

setwd("/Users/xinwang/Documents/Project/Proposal/Figure3/Version_20200511/Final_All/BRCA_figure/DDR_deficiency/Version_July/")
### for Amplication

library(ggplot2)
library(ggrepel)
Expression_var<-read.table("DDRM_DDRNM_DEG_pvalue.txt", header = T,sep = "\t", na.strings = NA, fill = T)
attach(Expression_var)
ggplot(Expression_var, aes(x=log2(Difference_2), y=Padjust)) +
  geom_point(size=1,color= ifelse(Difference_2>2 & Padjust>=1.3,"red",
                                    ifelse(Difference_2<=0.5 & Padjust>=1.3,"blue","gray")))+
  geom_text_repel(aes(label=ifelse((Difference_2>=3| Difference_2<= 0.1 ) & Padjust>=2,as.character(Gene),'')),hjust=0.5, 
                  size=2,segment.size=0.02)+
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )+
  theme(legend.justification = c(1, 1), 
        legend.text= element_text(face="italic"),legend.position = c(1,1))



#### updated new version
setwd("/Users/xinwang/Documents/Projects/Proposal/Figure4/OriginalDataForFig4/")

Mutationload_gene<-read.table("ALLSample_Divided_intoCate_final_immnueExp_forfigure.txt",header = T,sep = "\t", na.strings = NA, fill = T)

GeneExpressTr<- melt(Mutationload_gene,id.vars = c("Pathway2"),measure.vars = c("PDL1","CD8A","PD1","PDL2","CTLA4","CD8B"))
GeneExpressTr$Pathway2

attach(GeneExpressTr)
levels(factor(GeneExpressTr$Pathway2))
head(GeneExpressTr)
GeneExpressTr$variable<- factor(GeneExpressTr$variable, levels=levels(GeneExpressTr$variable)[c(3,1,4,2,6,5)],ordered=TRUE)
GeneExpressTr$Pathway2=factor(GeneExpressTr$Pathway2,levels = levels(factor(GeneExpressTr$Pathway2))[c(2,3,1,4)],ordered = TRUE)
GeneExpressTr$Pathway2
#Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(2,1,3)],ordered = TRUE)
#Mutationload_gene$Pathway
#Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(4,9,10,1,8,11,7,14,15,2,3,5,6,13,12)],ordered = TRUE)

ggplot(GeneExpressTr, aes(x=variable, y=log2(value),fill=GeneExpressTr$Pathway2)) +facet_wrap(.~variable, scale="free",nrow = 2)+ geom_boxplot(outlier.shape=NA,position=position_dodge(0.9)) + theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        #)+ theme(legend.position = "bottom") +scale_fill_manual(values=c("red","green","blue")) +ylim(0,12.5)
        
        
  )+ theme(legend.position = "bottom") +scale_fill_manual(values=c("brown1","wheat","darkgray","gray")) +ylim(0,14.5)



PDL1<-GeneExpressTr[GeneExpressTr$variable=="PD1",]
tail(PDL1)

levels(PDL1$Pathway2)

mutator<-PDL1[PDL1$Pathway2=="DDRM",3]
Nmutator<-PDL1[PDL1$Pathway2=="DDRNM",3]
mutant<-PDL1[PDL1$Pathway2=="DDR",3]
nonmutant<-PDL1[PDL1$Pathway2=="NonDDR",3]




##### Draw alll the genes under different Pathway

Mutationload_gene<-read.table("ALLSample_Divided_intoCate_final_immnueExp.txt",header = T,sep = "\t", na.strings = NA, fill = T)

GeneExpressTr<- melt(Mutationload_gene,id.vars = c("Pathway"),measure.vars = c("PDL1","CD8A","PD1","PDL2","CTLA4","CD8B"))
GeneExpressTr

attach(GeneExpressTr)

head(GeneExpressTr)
GeneExpressTr$variable<- factor(GeneExpressTr$variable, levels=levels(GeneExpressTr$variable)[c(3,1,4,2,6,5)],ordered=TRUE)
levels(Pathway)
GeneExpressTr$Pathway=factor(GeneExpressTr$Pathway,levels = levels(GeneExpressTr$Pathway)[c(5,6,11,12,1,10,13,9,16,17,2,3,7,8,15,4,14)],ordered = TRUE)

#Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(2,1,3)],ordered = TRUE)
#Mutationload_gene$Pathway
#Mutationload_gene$Pathway<-factor(Mutationload_gene$Pathway,levels = levels(Mutationload_gene$Pathway)[c(4,9,10,1,8,11,7,14,15,2,3,5,6,13,12)],ordered = TRUE)


ggplot(GeneExpressTr, aes(x=Pathway, y=log2(value),fill=variable)) +facet_wrap(.~variable, scale="free",nrow = 1)+ 
  geom_violin()+ 
  geom_boxplot(width=0.1,coef=1e30, fill="gray")+
  #geom_boxplot(outlier.shape=NA,position=position_dodge(0.8)) + 
  
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=10),
        axis.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        #)+ theme(legend.position = "bottom") +scale_fill_manual(values=c("red","green","blue")) +ylim(0,12.5)
        
        
  )+ coord_flip()+theme(legend.position = "bottom") +scale_fill_brewer(palette="Dark2") +ylim(0,14.5)

####
+scale_fill_manual(values=c("brown1","wheat", "orchid1","orchid1","orchid1",
                                                                  "orchid1","orchid1","orchid1","orchid1","orchid1",
                                                                  "orchid1","orchid1","orchid1","orchid1","orchid1",
                                                                  "darkgray","gray")) +ylim(0,14.5)




wilcox.test(mutator, nonmutant, alternative = "greater")
wilcox.test(mutant, nonmutant, alternative = "greater")
wilcox.test(mutator, mutant, alternative = "greater")
wilcox.test(mutator, Nmutator, alternative = "greater")



#### Draw the different RSEM across diffferent DDR pathway

GeneExpressTr2<- melt(Mutationload_gene,id.vars = c("Pathway"),measure.vars = c("PDL1","CD8A","PD1","PDL2","CTLA4","CD8B"))
GeneExpressTr2

attach(GeneExpressTr2)

head(GeneExpressTr2)
levels(GeneExpressTr2$variable)
GeneExpressTr2$variable<- factor(GeneExpressTr2$variable, levels=levels(GeneExpressTr2$variable)[c(3,1,4,2,6,5)],ordered=TRUE)
GeneExpressTr2$Pathway<-factor(GeneExpressTr2$Pathway,levels = levels(GeneExpressTr2$Pathway)[c(12,13,6,5,3,2,15,14,7,11,8,1,10,9,4)],ordered = TRUE)
#4,9,10,1,8,11,7,14,15,2,3,5,6,13,12)],ordered = TRUE)


### calculate the Pvalue of each genes


PDL1<-GeneExpressTr2[GeneExpressTr2$variable=="CTLA4",]
tail(PDL1)

levels(PDL1$Pathway)

for (i in (rev(levels(PDL1$Pathway)))){
  
  mutator<-PDL1[PDL1$Pathway=="DDRM",3]
  nonmutant<-PDL1[PDL1$Pathway==i,3]
  pvalue<-wilcox.test(mutator, nonmutant, alternative = "greater")$p.value
  pvalue
  print(paste(i,pvalue))
}



