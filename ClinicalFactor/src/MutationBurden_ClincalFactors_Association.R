setwd("/Users/xinwang/Documents/Projects/Proposal/Figure1/")
library(ggplot2)
library(reshape2)
Mutation<-read.table("TCGA_mutation_cna_clinical_combined3.txt",header = T,sep = "\t", na.strings = NA, fill = T)
attach(Mutation)

##############################
#########figure 1a ###########
##############################


############# print Pvalue of Mutation burden and associated factor into file

test<-read.table("SNVburden_Pvalue.txt",header = T, sep = "\t")

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
  scale_y_discrete(limits = rev(unique(test_melted$Type)))+
  scale_size(range = c(0,8), breaks = c(0,2,4,6,8))

###
###############
test<-read.table("Indelburden_Pvalue.txt",header = T, sep = "\t")

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
  scale_y_discrete(limits = rev(unique(test_melted$Type)))+
  scale_size(range = c(0,8), breaks = c(0,2,4,6,8))



### Draw the figure for Figure 1a CNV
test<-read.table("CNVburden_Pvalue.txt",header = T, sep = "\t")

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
  scale_y_discrete(limits = rev(unique(test_melted$Type)))+
  scale_size(range = c(0,8), breaks = c(0,2,4,6,8))



###### Calculate the SNV Burden ###

#### 37.59 Mb
Test<- Mutation[Sample_Status== "hypermutated sample" ||Sample_Status== "PASS" ,]
Mutationload =(Test$Flank3+Test$Flank5+Test$UTR5 + Test$UTR3 +Test$Missense_Mutation +Test$Nonsense_Mutation +Test$Nonstop_Mutation +Test$Silent +Test$Splice_Site +Test$Translation_Start_Site)


Mutationload_log = log10(Mutationload+0.1)
Mutation3<-cbind(Mutationload_log,Mutationload,Test)

Figure1BSNVdata<-cbind(Mutation3$ID, Mutation3$type, Mutationload_log)
write.table(as.matrix(Figure1BSNVdata),file="Fig1BSNV_datasouce.xls",row.names = TRUE, sep="\t")


############## draw figure 1b ##### mutation load ###

ggplot(Mutation3, aes(x=Mutation3$type, y = Mutation3$Mutationload_log))+ 
  geom_violin(fill="gray",size=0.1)+ 
  geom_boxplot(width=0.1, fill="white",coef=1e30)+
  #geom_boxplot(width=0.5, fill="red",coef=1e30)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + scale_y_continuous(position = "right")



###### Calculate the Indel burden ###

Test<- Mutation[Sample_Status== "hypermutated sample" ||Sample_Status== "PASS" ,]

#### modify age #### 


INDEL<- log10(Test$Frame_Shift_Del+Test$Frame_Shift_Ins+Test$In_Frame_Del +Test$In_Frame_Ins +0.1)
Mutation2<-cbind(INDEL,Test)

### print into a table for indel burden

Figure1BIndeldata<-cbind(Mutation2$ID, Mutation2$type, INDEL)
write.table(as.matrix(Figure1BIndeldata),file="Fig1BIndel_datasouce.xls",row.names = TRUE, sep="\t")


######################################Figure 1d
#### draw the Indel load boxplot across different cancer type: 

ggplot(Mutation2, aes(x=Mutation2$type, y = Mutation2$INDEL))+ 
  geom_violin(fill="gray",size=0.1)+ 
  geom_boxplot(width=0.1, fill="white",coef=1e30)+
 # geom_boxplot(width=0.5, fill="purple",coef=1e30)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +scale_y_continuous(position = "right")

######################################Figure 1e
#### draw the figure of CNV load in different cancer type
boxplot(log10(Test$CNV_nseg+0.1)~Test$type, data=Test, range =0.5, outline=FALSE,boxwex = 0.6,col="gray", whisklty = 1, cex=1) 
par(cex.axis=0.6)
help(boxplot )

### print into a table for indel burden

Figure1BCNVdata<-cbind(Test$ID, Test$type,log10(Test$CNV_nseg+0.1))
write.table(as.matrix(Figure1BCNVdata),file="Fig1BCNV_datasouce.xls",row.names = TRUE, sep="\t")


ggplot(Test, aes(x=Test$type, y = log10(Test$CNV_nseg+0.1)))+ 
  geom_violin(fill="gray",size=0.1)+ 
  #geom_boxplot(outlier.shape=NA, fill="blue")+
  #geom_boxplot(width=0.5, fill="blue",coef=1e30)+
 geom_boxplot(width=0.1, fill="white",coef=1e30)+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+scale_y_continuous(position = "right")


dev.off()






