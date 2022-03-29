#### Calculate the different factor influence ####

setwd("/Users/xinwang/Documents/Projects/Proposal/Figure1")
library(ggplot2)
library(reshape2)
Mutation<-read.table("TCGA_mutation_cna_clinical_combined3.txt",header = T,sep = "\t", na.strings = NA, fill = T)
attach(Mutation)

####### Here we need to remove some samples and modified the age 

## 
#Age_modified<-  ceiling(Mutation$age_at_initial_pathologic_diagnosis/10)*10
#levels(as.factor(Age_modified))
#density(na.omit(Age_modified))

#plot(density(na.omit(Age_modified)))


### modified according to age distribution on  https://www.cancer.gov/about-cancer/causes-prevention/risk/age
# AgeM <-c()
# for (i in Mutation$age_at_initial_pathologic_diagnosis) {
#   
#   #AgeM<-append(AgeM,ceiling(i/10)*10)
#   print (AgeM)
#  if (is.null(i) || i==''){
#     AgeM<- append(AgeM, NA)
#   } else if (i<20){
#     AgeM<-append(AgeM, 20)
#   } else if (i>=20 && i<35){
#     AgeM<-append(AgeM, 30)
#   } else if (i>=35 && i<45){
#     AgeM<-append(AgeM, 40)
#   } else if (i>=45 && i<55){
#     AgeM<-append(AgeM, 50)
#   } else if (i>=55 && i<65){
#     AgeM<-append(AgeM, 60)
#   } else if (i>=65 && i<75){
#     AgeM<-append(AgeM, 70)
#   } else if (i>=75 && i<85){
#     AgeM<-append(AgeM, 80)
#   } else if (i>=85){
#       AgeM<-append(AgeM, 85)
#    }
# }



cancer_type<- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")




#### This is to calculate the CNV load, here we used CNV burden as log10 (CNV number of segments)
#### Here we removed the samples with RNA degradation

Test<- Mutation[Sample_Status== "hypermutated sample" ||Sample_Status== "PASS" ,]


# Test data using the specific cancer type
# 
# 
# Test1<-Test[type=="BRCA",]
# # #test_age<-ceiling(Test1$age_at_initial_pathologic_diagnosis/10)*10
# # testage<-aov(x=factor(Test1$Age_modified),y=log10(Test1$CNV_nsegs+0.1),method = 'spearman')$p.value# result_age<-aov(log10(Test1$CNV_nsegs+0.1)~(ceiling(Test1$age_at_initial_pathologic_diagnosis/10)*10), data=Test1)
# # #Parameter<-c("age_at_initial_pathologic_diagnosis","gender","race","ajcc_pathologic_tumor_stage")
# result_tumorstage<-aov(log10(Test1$CNV_nsegs+0.1)~factor(Test1$Age_modified), data=Test1)
# 
# test_tumorstage<-summary(result_tumorstage)[[1]]$'Pr(>F)'[1]
# print (paste("BRCA",test_tumorstage))




for (i in cancer_type) {
 # print (i)
  SingleTumor<- Test[type==i,]
  
  ##### This part is for age
  result_age<-aov(log10(SingleTumor$CNV_nsegs+1)~factor(SingleTumor$Age_modified), data=SingleTumor)
 
   p_age<-summary(result_age)[[1]]$'Pr(>F)'[1]
  
   #### this part is for race
   result_race<-aov(log10(SingleTumor$CNV_nsegs+1)~SingleTumor$race, data=SingleTumor)
   p_race<-summary(result_race)[[1]]$'Pr(>F)'[1]
   
   
   
  #####this part is for the gender
  if (i=="CESC"|| i=="OV" ||i=="PRAD" ||i=="TGCT" || i=="UCEC" || i=="UCS"){
    p_gender<- c("NA")
  }else{
    result_gender<-aov(log10(SingleTumor$CNV_nsegs+1)~SingleTumor$gender, data=SingleTumor)
    p_gender<-summary(result_gender)[[1]]$'Pr(>F)'[1]
  }


  ###### this part is for the tumor stage because some tumor they don't have stage information
  if (i=="CESC" || i== "DLBC" || i=="OV" || i=="THYM" || i=="UCEC" || i=="UCS"){
    result_tumorstage<-aov(log10(SingleTumor$CNV_nsegs+1)~SingleTumor$clinical_stage, data=SingleTumor)
    p_tumorstage<-summary(result_tumorstage)[[1]]$'Pr(>F)'[1]
    
  }else if (i=="GBM" || i=="LAML" || i=="LGG" || i=="PCPG" ||i=="PRAD" || i=="SARC"){
    p_tumorstage<- c("NA")
  }else{
    result_tumorstage<-aov(log10(SingleTumor$CNV_nsegs+1)~SingleTumor$ajcc_pathologic_tumor_stage, data=SingleTumor)
    p_tumorstage<-summary(result_tumorstage)[[1]]$'Pr(>F)'[1]
    
  }
#print (paste(i,p_gender,p_tumorstage))
# result<-rbind (result,test)
  print (paste(i,p_age,p_gender,p_race,p_tumorstage))
}

Pvalue_mutation<-data.frame(Tumor_type,Age_Pvalue,Gender_Pvalue,Race_Pvalue,Tumor_Pvalue)
write.table(Pvalue_mutation, file = "mutation_test_pvalue.txt",sep = "\t")


######################################Figure 1e
#### draw the figure of CNV load in different cancer type
boxplot(log10(Test$CNV_nseg+1)~Test$type, data=Test, range =0.5, outline=FALSE,boxwex = 0.6,col="gray", whisklty = 1, cex=1) 
par(cex.axis=0.6)
help(boxplot )

ggplot(Test, aes(x=Test$type, y = log10(Test$CNV_nseg+1)))+ 
  #geom_violin(color="blue",fill="gray")+ 
  #geom_boxplot(outlier.shape=NA, fill="blue")+
  geom_boxplot(width=0.5, fill="blue",coef=1e30)+
 # geom_boxplot(width=0.1, fill="white",coef=1e30)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

#############

test<-read.table("CNVload_Pvalue_Allfactors.txt",header = T, sep = "\t")

test_melted<-melt(test,id="Type")
head(test_melted)

p<-ggplot(test_melted, aes(x=variable, y=Type))

p+geom_point( aes(size=value),colour= ifelse(test_melted$value>1.301029996,"blue",ifelse(test_melted$value==0,"white","gray")), fill= ifelse(test_melted$value>1.301029996,"red",ifelse(test_melted$value==0,"white","gray"))) +
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        
        
  ) + 
  #### reverse the order of y axis
  scale_y_discrete(limits = rev(levels(test_melted$Type)))+
  scale_size_area(max_size=5)
#Add labels to axes



###### Calculate the Indel load ###

Test<- Mutation[Sample_Status== "hypermutated sample" ||Sample_Status== "PASS" ,]

#### modify age #### 


INDEL<- log10(Test$Frame_Shift_Del+Test$Frame_Shift_Ins+Test$In_Frame_Del +Test$In_Frame_Ins +0.1)
Mutation2<-cbind(INDEL,Test)

######################################Figure 1d
#### draw the Indel load boxplot across different cancer type: 

ggplot(Mutation2, aes(x=Mutation2$type, y = Mutation2$INDEL))+ 
 # geom_violin(fill="gray",color="purple")+ 
  #geom_boxplot(width=0.1, fill="white",coef=1e30)+
  geom_boxplot(width=0.5, fill="purple",coef=1e30)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

#### calcualte the Pvalue of different factors across samples for Indel load
#result<-c()

help(anova)
Tumor_type<-c()
Age_Pvalue<-c()
Race_Pvalue<-c()
Tumor_Pvalue<-c()
Gender_Pvalue<-c()
Tumor_size<-c()

for (i in cancer_type) {
  # print (i)
  SingleTumor<- Mutation2[type==i,]
  
  Rnumber<- nrow(SingleTumor)
  Tumor_size<-append(Tumor_size,Rnumber)
  
  
  
  #####this part is for the age
  result_age<-aov(SingleTumor$INDEL~factor(SingleTumor$Age_modified), data=SingleTumor)
  
  p_age<-summary(result_age)[[1]]$'Pr(>F)'[1]
  
  result_age<-cor.test(x=SingleTumor$age_at_initial_pathologic_diagnosis,y=SingleTumor$INDEL,method = 'spearman')
  p_age<-result_age$p.value 
  
  #####this part is for the gender
  if (i=="CESC"|| i=="OV" ||i=="PRAD" ||i=="TGCT" || i=="UCEC" || i=="UCS"){
    p_gender<- c("NA")
  }else{
    result_gender<-aov(SingleTumor$INDEL~SingleTumor$gender, data=SingleTumor)
    p_gender<-summary(result_gender)[[1]]$'Pr(>F)'[1]
  }
  

  #####this part is for the race 
  
  result_race<-aov(SingleTumor$INDEL~SingleTumor$race, data=SingleTumor)
  p_race<-summary(result_race)[[1]]$'Pr(>F)'[1]
  
 
  
  ###### this part is for the tumor stage because some tumor they don't have stage information
  if (i=="CESC" || i== "DLBC" || i=="OV" || i=="THYM" || i=="UCEC" || i=="UCS"){
    result_tumorstage<-aov(SingleTumor$INDEL~SingleTumor$clinical_stage, data=SingleTumor)
    p_tumorstage<-summary(result_tumorstage)[[1]]$'Pr(>F)'[1]
    
  }else if (i=="GBM" || i=="LAML" || i=="LGG" || i=="PCPG" ||i=="PRAD" || i=="SARC"){
    p_tumorstage<- c("NA")
  }else{
    result_tumorstage<-aov(SingleTumor$INDEL~SingleTumor$ajcc_pathologic_tumor_stage, data=SingleTumor)
    p_tumorstage<-summary(result_tumorstage)[[1]]$'Pr(>F)'[1]
    
  }
  #print (paste(i,p_gender,p_tumorstage))
  # result<-rbind (result,test)
  print (paste(i,p_age,p_gender,p_race,p_tumorstage))
  
  Tumor_type<-append(Tumor_type,i)
  Age_Pvalue<-append(Age_Pvalue,p_age)
  Gender_Pvalue<-append(Gender_Pvalue,p_gender)
  Race_Pvalue<-append(Race_Pvalue,p_race)
  Tumor_Pvalue<-append(Tumor_Pvalue,p_tumorstage)
  # test<-paste(i,p_age,p_gender,p_race,p_tumorstage)
  # result<-rbind (result,test)
}

Pvalue_Indel<-cbind(Tumor_type,Age_Pvalue,Gender_Pvalue,Race_Pvalue,Tumor_Pvalue)


###############
test<-read.table("Indelload_Pvalue_v0815.txt",header = T, sep = "\t")

test_melted<-melt(test,id="Type")
head(test_melted)

p<-ggplot(test_melted, aes(x=variable, y=Type))

p+geom_point( aes(size=value),colour= ifelse(test_melted$value>1.301029996,"purple",ifelse(test_melted$value==0,"white","gray")), fill= ifelse(test_melted$value>1.301029996,"red",ifelse(test_melted$value==0,"white","gray"))) +
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        
        
  ) + 
  #### reverse the order of y axis
  scale_y_discrete(limits = rev(levels(test_melted$Type)))+
  scale_size(range = c(0,8), breaks = c(0,2,4,6,8))


#scale_size_area(max_size=5)+
 # scale_fill_gradientn(limits = c(0,9), breaks = c(0, 3, 6, 9),labels=c(0, 3, 6, 9), colors = black)   
#Add labels to axes

help()






Tumor_information<-cbind(Tumor_type,Tumor_size)
write.table(as.matrix(Tumor_information),file="Fig1C_SampleSize_datasouce.xls",row.names = TRUE, sep="\t")

#### plot sample size within the cancer type
################################ Figure 1b ###########
barplot(Tumor_size, names.arg = Tumor_type,width=0.5)
help("barplot")

######################################################
#### calculate the mutation load  #####
#######################################

###

#### 37.59 Mb
Mutationload =(Test$Flank3+Test$Flank5+Test$UTR5 + Test$UTR3 +Test$Missense_Mutation +Test$Nonsense_Mutation +Test$Nonstop_Mutation +Test$Silent +Test$Splice_Site +Test$Translation_Start_Site)


Mutationload_log = log10(Mutationload+0.1)
Mutation3<-cbind(Mutationload_log,Mutationload,Test)



############## draw figure 1b ##### mutation load ###

ggplot(Mutation3, aes(x=Mutation3$type, y = Mutation3$Mutationload_log))+ 
  geom_violin(fill="gray",color="red")+ 
 geom_boxplot(width=0.1, fill="white",coef=1e30)+
  #geom_boxplot(width=0.5, fill="red",coef=1e30)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

#########################################################



Tumor_type<-c()
Age_Pvalue<-c()
Race_Pvalue<-c()
Tumor_Pvalue<-c()
Gender_Pvalue<-c()
Tumor_size<-c()

for (i in cancer_type) {
  # print (i)
  SingleTumor<- Mutation3[type==i,]
  
  Rnumber<- nrow(SingleTumor)
  Tumor_size<-append(Tumor_size,Rnumber)
  
  
  
  #####this part is for the age
  result_age<-aov(SingleTumor$Mutationload_log~factor(SingleTumor$Age_modified), data=SingleTumor)
  
  p_age<-summary(result_age)[[1]]$'Pr(>F)'[1]
  
 # result_age<-cor.test(x=SingleTumor$age_at_initial_pathologic_diagnosis,y=SingleTumor$INDEL,method = 'spearman')
 # p_age<-result_age$p.value 
  
  #####this part is for the gender
  if (i=="CESC"|| i=="OV" ||i=="PRAD" ||i=="TGCT" || i=="UCEC" || i=="UCS"){
    p_gender<- c("NA")
  }else{
    result_gender<-aov(SingleTumor$Mutationload_log~SingleTumor$gender, data=SingleTumor)
    p_gender<-summary(result_gender)[[1]]$'Pr(>F)'[1]
  }
  
  
  #####this part is for the race 
  
  result_race<-aov(SingleTumor$Mutationload_log~SingleTumor$race, data=SingleTumor)
  p_race<-summary(result_race)[[1]]$'Pr(>F)'[1]
  
  
  
  ###### this part is for the tumor stage because some tumor they don't have stage information
  if (i=="CESC" || i== "DLBC" || i=="OV" || i=="THYM" || i=="UCEC" || i=="UCS"){
    result_tumorstage<-aov(SingleTumor$Mutationload_log~SingleTumor$clinical_stage, data=SingleTumor)
    p_tumorstage<-summary(result_tumorstage)[[1]]$'Pr(>F)'[1]
    
  }else if (i=="GBM" || i=="LAML" || i=="LGG" || i=="PCPG" ||i=="PRAD" || i=="SARC"){
    p_tumorstage<- c("NA")
  }else{
    result_tumorstage<-aov(SingleTumor$Mutationload_log~SingleTumor$ajcc_pathologic_tumor_stage, data=SingleTumor)
    p_tumorstage<-summary(result_tumorstage)[[1]]$'Pr(>F)'[1]
    
  }
  #print (paste(i,p_gender,p_tumorstage))
  # result<-rbind (result,test)
  print (paste(i,p_age,p_gender,p_race,p_tumorstage))
  
  Tumor_type<-append(Tumor_type,i)
  Age_Pvalue<-append(Age_Pvalue,p_age)
  Gender_Pvalue<-append(Gender_Pvalue,p_gender)
  Race_Pvalue<-append(Race_Pvalue,p_race)
  Tumor_Pvalue<-append(Tumor_Pvalue,p_tumorstage)
  # test<-paste(i,p_age,p_gender,p_race,p_tumorstage)
  # result<-rbind (result,test)
}


###### print Tumor mutation load Pvalue

# ######create a pvalue matrix and put into file
# Pvalue_mutation<-data.frame(Tumor_type,Age_Pvalue,Gender_Pvalue,Race_Pvalue,Tumor_Pvalue)
# #col.names(Pvalue_mutation) <-Tumor_type
# melter_Pvalue<- na.omit(melt(Pvalue_mutation,id = "Tumor_type"))
# 
# ggplot(melter_Pvalue, aes(x=variable, y=Tumor_type))+ geom_point(aes(size=value,colour= ifelse(melter_Pvalue$value<0.05,"red","gray"))) +
#   theme_bw()+
#   theme(plot.background = element_blank(),
#         text = element_text(size=8),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank()
#   ) +
#   scale_size_area(max_size=5)
# #Add labels to axes

##############################
#########figure 1a ###########
##############################
############# print Pvalue of Mutation load and associated factor into file
Pvalue_mutation<-data.frame(Tumor_type,Age_Pvalue,Gender_Pvalue,Race_Pvalue,Tumor_Pvalue)
write.table(Pvalue_mutation, file = "mutation_test_pvalue.txt",sep = "\t")

#### modified the file because of NA information -log10(..) NA =0



test<-read.table("Mutationload_Pvalue_v0815.txt",header = T, sep = "\t")

test_melted<-melt(test,id="Type")
head(test_melted)

p<-ggplot(test_melted, aes(x=variable, y=Type))

p+geom_point( aes(size=value),colour= ifelse(test_melted$value>1.301029996,"red",ifelse(test_melted$value==0,"white","gray")), fill= ifelse(test_melted$value>1.301029996,"red",ifelse(test_melted$value==0,"white","gray"))) +
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      
     
  ) + 
  #### reverse the order of y axis
  scale_y_discrete(limits = rev(levels(test_melted$Type)))+
  scale_size(range = c(0,8), breaks = c(0,2,4,6,8))
  
 
+
  scale_fill_continuous(limits = c(0,9), breaks = c(0, 3, 6, 9))         #Add labels to axes


 
 #### CNV 
 
 test<-read.table("CNVload_Pvalue_Allfactors_v0815.txt",header = T, sep = "\t")
 
 test_melted<-melt(test,id="Type")
 head(test_melted)
 
 p<-ggplot(test_melted, aes(x=variable, y=Type))
 
 p+geom_point( aes(size=value),colour= ifelse(test_melted$value>1.301029996,"blue",ifelse(test_melted$value==0,"white","gray")), fill= ifelse(test_melted$value>1.301029996,"red",ifelse(test_melted$value==0,"white","gray"))) +
   theme_bw()+
   theme(plot.background = element_blank(),
         text = element_text(size=8),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank()
         
         
   ) + 
   #### reverse the order of y axis
   scale_y_discrete(limits = rev(levels(test_melted$Type)))+
   scale_size(range = c(0,8), breaks = c(0,2,4,6,8))


 
 
 #### Some examples for the clinical risk factors
 
 setwd("/Users/xinwang/Documents/Database/Published_TCGA/ATLAS_version1011/modified_clinical_version1016/")
 #myfiles<-list.files(pattern="TCGA_mutation_clinic_combined.*.txt.modified.txt")
 #myfiles
 
