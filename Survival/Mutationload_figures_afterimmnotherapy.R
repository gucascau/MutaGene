


library(ggplot2)
setwd("/Users/xinwang/Documents/Project/Proposal/Figure3/Version_1103/Survival_analysis_afterimmnuotherapy/BRCA_mutationload_output/")

myfiles<-list.files(pattern="*.txt")


length(myfiles)


for (i in 1:40) {
  ###i=8, 37
i=4
  myfiles[i]
  Mutationload_genes<-read.table(myfiles[i],header = T,sep = "\t", na.strings = NA, fill = T)
  mutant<-Mutationload_genes[Mutationload_genes$x=="Mutant",4]
  nonmutant<-Mutationload_genes[Mutationload_genes$x=="NonMutant",4]
  wilcox.test(mutant,nonmutant,alternative = "greater")
  
  
 pdf(paste(myfiles[i],".mutationload.pdf",sep=""))
 #attach(Mutationload_genes)
  ggplot(Mutationload_genes, aes(x=Mutationload_genes$x, y =log10(Mutationload_genes$Mutationload),fill=x ))+ 
    #geom_violin(fill="gray")+ 
    #geom_jitter(shape=16, position=position_jitter(0.2),col="gray",size=1)+
   # geom_boxplot(width=0.1, fill="white",coef=1e30)+
    
   geom_violin(trim=F)+ geom_boxplot(width=0.2, fill="gray",coef=1e30)+
    
    #geom_boxplot(outlier.shape=NA, fill=c("red","blue"),alpha = 0.8,width=0.8)+ 
    theme_bw()+ scale_fill_manual(values=c("#FF0000", "#0000FF"))+
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  dev.off()

  
  
  
  
  
  }
  
  
  

  
  
