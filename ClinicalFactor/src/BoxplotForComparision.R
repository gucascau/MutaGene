
### Here is to draw figure 1c of gender difference for STAD
library(ggplot2)
library(reshape2)



setwd("/Users/xinwang/Documents/Projects/Proposal/Figure1/")

Mutation<-read.table("TCGA_mutation_cna_clinical_combined3.txt",header = T,sep = "\t", na.strings = NA, fill = T)
attach(Mutation)


#### draw the figure of mutation load in SKCM under Gender ###

Test<- Mutation[Sample_Status== "hypermutated sample" ||Sample_Status== "PASS" ,]

SingleTumor<-Test[type=="SKCM",]

Mutationload =(SingleTumor$Flank3+SingleTumor$Flank5+SingleTumor$UTR5 + SingleTumor$UTR3 +SingleTumor$Missense_Mutation +SingleTumor$Nonsense_Mutation +SingleTumor$Nonstop_Mutation +SingleTumor$Silent +SingleTumor$Splice_Site +SingleTumor$Translation_Start_Site)

result_gender<-wilcox.test(log10(Mutationload+0.1)~SingleTumor$gender, data=SingleTumor)
result_gender
p_gender

Figure1Ddata<-cbind(SingleTumor$ID, SingleTumor$gender,log10(Mutationload+0.1))
write.table(as.matrix(Figure1Ddata),file="Fig1D_datasouce.xls",row.names = TRUE, sep="\t")

boxplot(log10(Mutationload+0.1)~SingleTumor$gender, data=SingleTumor,  range =0.5,boxwex = 0.5, outline=FALSE,col="white", cex=1,ylim=c(1.5,4)) 

#### draw the figure of mutation load in THCA under different tumor stages ###

Test<- Mutation[Sample_Status== "hypermutated sample" ||Sample_Status== "PASS" ,]

SingleTumor<-Test[type=="THCA",]
SingleTumor$ID
SingleTumor$ajcc_pathologic_tumor_stage
SingleTumor$type
Mutationload =(SingleTumor$Flank3+SingleTumor$Flank5+SingleTumor$UTR5 + SingleTumor$UTR3 +SingleTumor$Missense_Mutation +SingleTumor$Nonsense_Mutation +SingleTumor$Nonstop_Mutation +SingleTumor$Silent +SingleTumor$Splice_Site +SingleTumor$Translation_Start_Site)

result_gender<-aov(log10(Mutationload+0.1)~SingleTumor$ajcc_pathologic_tumor_stage, data=SingleTumor)
p_gender<-summary(result_gender)[[1]]$'Pr(>F)'[1]
p_gender

Figure1Fdata<-cbind(SingleTumor$ID,SingleTumor$ajcc_pathologic_tumor_stage,log10(Mutationload+0.1))
write.table(as.matrix(Figure1Fdata),file="Fig1F_datasouce.xls",row.names = TRUE, sep="\t")

boxplot(log10(Mutationload+0.1)~SingleTumor$ajcc_pathologic_tumor_stage, data=SingleTumor,  range =0.5,boxwex = 0.5, outline=FALSE,col="white", cex=1,ylim=c(0.5,2)) 

#### draw the figure of mutation load in BLCA under different tumor stages ###

Test<- Mutation[Sample_Status== "hypermutated sample" ||Sample_Status== "PASS" ,]

SingleTumor<-Test[type=="BLCA" & (race =="ASIAN" || race =="WHITE" || race == "BLACK OR AFRICAN AMERICAN"),]


SingleTumor<-read.table("Tumor_Race_forFigure_BLCA.txt",header = T,sep = "\t", na.strings = NA, fill = T)
SingleTumor$ID
Mutationload =(SingleTumor$Flank3+SingleTumor$Flank5+SingleTumor$UTR5 + SingleTumor$UTR3 +SingleTumor$Missense_Mutation +SingleTumor$Nonsense_Mutation +SingleTumor$Nonstop_Mutation +SingleTumor$Silent +SingleTumor$Splice_Site +SingleTumor$Translation_Start_Site)

result_gender<-aov(log10(Mutationload+0.1)~SingleTumor$race, data=SingleTumor)
p_gender<-summary(result_gender)[[1]]$'Pr(>F)'[1]
p_gender
Figure1Edata<-cbind(SingleTumor$ID,SingleTumor$race,log10(Mutationload+0.1))
write.table(as.matrix(Figure1Edata),file="Fig1E_datasouce.xls",row.names = TRUE, sep="\t")


boxplot(log10(Mutationload+0.1)~SingleTumor$race, data=SingleTumor,  range =0.5,boxwex = 0.5, outline=FALSE,col="white", cex=1, ylim=c(1,3)) 

dev.off()



### draw the figure of muation load in LGG under the radiation
setwd("/Users/xinwang/Documents/Database/Publish_Data/Published_TCGA/ATLAS_version1011/modified_clinical_version1016/Radiation_usedForPvalue")

SingleTumor<-read.table("Clinical_radiation.lgg.RadiationALL.txt",header = T,sep = "\t", na.strings = NA, fill = T)

SingleTumor$CNV_ID
### anova test
result_radiation<-aov(log10(Mutationload+0.1)~SingleTumor$Status, data=SingleTumor)

## wilcox test
result_radiation2<-wilcox.test(log10(Mutationload+0.1)~SingleTumor$Status, data=SingleTumor)
result_radiation2$p.value
p_radition<-summary(result_radiation)[[1]]$'Pr(>F)'[1]
p_radition

Figure1Gdata<-cbind(SingleTumor$ID,SingleTumor$Status,log10(Mutationload+0.1) )
write.table(as.matrix(Figure1Gdata),file="Fig1G_datasouce.xls",row.names = TRUE, sep="\t")

boxplot(log10(Mutationload+0.1)~SingleTumor$Status, data=SingleTumor,  range =0.5,boxwex = 0.5, outline=FALSE,col="white", cex=1, ylim=c(1,2)) 


### draw the figure of muation load in LUAD under the smoking
setwd("/Users/xinwang/Documents/Database/Publish_Data/Published_TCGA/ATLAS_version1011/modified_clinical_version1016/Smoking_usedForPvalue/")

SingleTumor<-read.table("TCGA_Smoking.LUAD.smokingComb.txt",header = T,sep = "\t", na.strings = NA, fill = T,row.names=NULL)
SingleTumor$Status


results_smoking<-aov(log10(Mutationload+0.1)~SingleTumor$Status, data=SingleTumor)
p_smoking<-summary(results_smoking)[[1]]$'Pr(>F)'[1]
p_smoking


result_smoking<-wilcox.test(log10(Mutationload+0.1)~SingleTumor$Status, data=SingleTumor)

result_smoking$p.value


Figure1Hdata<-cbind(SingleTumor$ID,SingleTumor$Status,log10(Mutationload+0.1) )
write.table(as.matrix(Figure1Hdata),file="Fig1H_datasouce.xls",row.names = TRUE, sep="\t")


boxplot(log10(Mutationload+0.1)~SingleTumor$Status, data=SingleTumor,  range =0.5,boxwex = 0.5, outline=FALSE,col="white", cex=1,ylim=c(1,3.5)) 




