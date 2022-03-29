#### This script is to draw the mutation distribution of DDR genes

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")

BiocManager::install("rngtools")
library("maftools")


BiocManager::install("maftools")

#BiocManager::install("PoisonAlien/maftools")

library(maftools)

??oncoplots

help("oncoplot")
### for 
setwd("/Users/xinwang/Documents/Projects/Proposal/Figure4/V20211205/")

#brac_maf=system.file('./','data_mutations_extended.txt',package = 'maftools')


#brac_clinc=system.file('extdata','data_bcr_clinical_data_patient.txt',package = 'maftools')

Brca_maf=read.maf(maf = "BRCA_req.allSampRemoveExtr.txt", clinicalData = "BRCA_req.allclinic2.txt")


getClinicalData(x=Brca_maf)

getSampleSummary(Brca_maf)

help(plotmafSummary)

plotmafSummary(maf = Brca_maf, rmOutlier = TRUE, addStat='median', dashboard=TRUE, titvRaw=FALSE)

oncoplot(maf=Brca_maf, top=10)

##### Draw the picture of all DDR genes mutation in breast cancer

DDR_ALL<-c("MSH2","MSH3","MSH6","MLH1","PMS2","MSH4","MSH5","MLH3","PMS1","PMS2P3","HFM1","XPC","RAD23B","CETN2","RAD23A","XPA","DDB1","DDB2","RPA1","RPA2","RPA3","TFIIH","ERCC3","ERCC2","GTF2H1","GTF2H2","GTF2H3","GTF2H4","GTF2H5","GTF2E2","CDK7","CCNH","MNAT1","ERCC5","ERCC1","ERCC4","LIG1","ERCC8","ERCC6","UVSSA","XAB2","MMS19","UNG","SMUG1","MBD4","TDG","OGG1","MUTYH","NTHL1","MPG","NEIL1","NEIL2","NEIL3","APEX1","APEX2","LIG3","XRCC1","PNKP","APLF","HMCES","RAD51","RAD51B","RAD51D","HELQ","SWI5","SWSAP1","ZSWIM7","SPIDR","PDS5B","DMC1","XRCC2","XRCC3","RAD52","RAD54L","RAD54B","BRCA1","BARD1","ABRAXAS1","PAXIP1","SMC5","SMC6","SHLD1","SHLD2","SHLD3","SEM1","RAD50","MRE11A","NBN","RBBP8","MUS81","EME1","EME2","SLX1A","SLX1B","GEN1","XRCC6","XRCC5","PRKDC","LIG4","XRCC4","DCLRE1C","NHEJ1","FANCA","FANCB","FANCC","BRCA2","FANCD2","FANCE","FANCF","FANCG","FANCI","BRIP1","FANCL","FANCM","PALB2","RAD51C","SLX4(FANCP)","FAAP20","FAAP24","FAAP100","UBE2T","POLA1","POLB","POLD1","POLD2","POLD3","POLD4","POLE","POLE2","POLE3","POLE4","REV3L","MAD2L2","REV1","POLG","POLH","POLI","POLQ","POLK","POLL","POLM","POLN","PRIMPOL","DNTT","PARP1","PARP2","PARP3","PARG","PARPBP","MGMT","ALKBH2","ALKBH3","TDP1","TDP2","SPRTN","NUDT1","DUT","RRM2B","PARK7","DNPH1","NUDT15","NUDT18","FEN1","FAN1","TREX1","TREX2","EXO1","APTX","SPO11","ENDOV","DNA2","DCLRE1A","DCLRE1B","EXO5","UBE2A","UBE2B","RAD18","SHPRH","HLTF","RNF168","RNF8","RNF4","UBE2V2","UBE2N","USP1","WDR48","HERC2","H2AX","CHAF1A","SETMAR","ATRX","BLM","RMI1","TOP3A","WRN","RECQL4","ATM","MPLKIP","RPA4","PRPF19","RECQL","RECQL5","RDM1","NABP2","ATR","ATRIP","MDC1","PCNA","RAD1","RAD9A","HUS1","RAD17","CHEK1","CHEK2","TP53BP1","RIF1","TOPBP1","CLK2","PER1")

oncoplot(maf=Brca_maf, genes = DDR_ALL, draw_titv = TRUE,keepGeneOrder =TRUE)

plotmafSummary(maf = Brca_maf, rmOutlier = TRUE, addStat='median', dashboard=TRUE, titvRaw=FALSE)

help(oncoplot)


#### Selected DDR driver genes



DDR_mutators_BRCA<-c ("RNF168","TOP3A","ALKBH1","DCLRE1B","RAD18","PMS2","FANCI","POLN","XPC","REC8","MSH4"
                      
                      ,"MSH2","MSH6","LIG1","USP7","RFC1","UBC","ERCC2","MMS19","ERCC6","RAD51","ZFYVE26","NIPBL","TONSL","SMC5","SPIDR","CDK9","EYA1","UBR5","BRCA2","BRCA1")



DDR_mutators_BRC_Tp53<-c ("RNF168","DCLRE1B","FANCI","POLN","ALKBH1","TOP3A","PMS2","XPC","RAD18","REC8","MSH4","MSH2","MSH6","LIG1","USP7","RFC1","UBC","ERCC2","MMS19","ERCC6","RAD51","ZFYVE26","NIPBL","TONSL","SMC5","SPIDR","CDK9","EYA1","UBR5","BRCA2","BRCA1","TP53")

Testpathway = data.frame(
  Genes = c(
    "RNF168","DCLRE1B","FANCI","POLN","ALKBH1","TOP3A","PMS2","XPC","RAD18","REC8","MSH2","MSH6","PMS2","LIG1","XPC","USP7","RFC1","UBC","ERCC2","MMS19","ERCC6","RAD51","ZFYVE26","NIPBL","TONSL","SMC5","SPIDR","CDK9","EYA1","UBR5","BRCA2","BRCA1"
  ),
  Pathway = rep(c(
    "Mutator", "OtherMutator"
  ), c(10,22)),
  stringsAsFactors = FALSE
)



DDR_mutators_TP53<-c ("MSH2","MSH6","PMS2","LIG1","XPC","USP7","RFC1","UBC","ERCC2","MMS19","ERCC6","DCLRE1B","RNF168","FANCI","RAD51","POLN","ZFYVE26","NIPBL","TONSL","SMC5","SPIDR","ALKBH1","TOP3A","CDK9","EYA1","UBR5","BRCA2","BRCA1","TP53")


#aml_genes_vaf = subsetMaf(maf = Brca_maf, genes = aml_genes, fields = "i_TumorVAF_WU", mafObj = FALSE)[,mean(i_TumorVAF_WU, na.rm = TRUE), Hugo_Symbol]
Brca_maf@variants.per.sample

sample_order = Brca_maf@variants.per.sample$Tumor_Sample_Barcode

help(plotmafSummary)
#Brca_maf@variants.per.sample$Tumor_Sample_Barcode

oncoplot(maf=Brca_maf, genes = DDR_mutators_BRCA, draw_titv = TRUE,keepGeneOrder =TRUE, sampleOrder = sample_order,logColBar=TRUE )

help("oncoplot")

oncoplot(maf=Brca_maf, genes = DDR_mutators, draw_titv = TRUE,sortByAnnotation=TRUE,keepGeneOrder =TRUE )


help("oncoplot")
dev.off()
somaticInteractions(maf = Brca_maf, genes = DDR_mutators_BRC_Tp53, pvalue=c(0.05,0.1))

help(somaticInteractions)


aml_genes=c("BRCA2","ERCC5","MSH2","MSH6","TERT","TOP3A")

DDR_genes=c("APLF","APTX","ASCC3","DNTT","LIG1","LIG3","LIG4","MRE11A","NBN","NHEJ1","PARG","PARP1","PARP3","PARPBP","PNKP","POLB","POLL","POLM","PRKDC","RAD50","RNF168","RNF8","TP53BP1","XRCC1","XRCC2","XRCC3","XRCC4","XRCC5","XRCC6","UBE2A","EXO1","HMGB1","MLH1","MLH3","MSH2","MSH3","MSH6","PCNA","PMS1","PMS2","POLD1","POLD2","POLD3","POLD4","RFC1","RFC2","RFC3","RFC4","RFC5","RPA1","RPA2","RPA3","RPA4","ALKBH1","ALKBH2","ALKBH3","APEX1","APEX2","APITD1","ATM","ATR","ATRIP","ATRX","BARD1","BLM","BRE","BRIP1","CCNH","CDK7","CETN2","CHAF1A","CHEK1","CHEK2","CLK2","CUL3","CUL4A","CUL5","DCLRE1A","DCLRE1B","DCLRE1C","DDB1","DDB2","DMC1","DNA2","DUT","EID3","EME1","EME2","ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC6","ERCC8","FAAP100","FAAP24","FAAP20","FAM175A","FAN1","FANCA","FANCB","FANCC","FANCD2","FANCE","FANCF","FANCG","FANCI","FANCL","FANCM","FEN1","GADD45A","GADD45G","GEN1","GTF2H1","GTF2H2","GTF2H3","GTF2H4","GTF2H5","H2AFX","HELQ","HES1","HFM1","HLTF","HMGB2","HUS1","INO80","KAT5","MAD2L2","MBD4","MDC1","MGMT","MMS19","MNAT1","MPG","MPLKIP","MRPL40","MUS81","MUTYH","NABP2","NEIL1","NEIL2","NEIL3","NFATC2IP","NSMCE1","NSMCE2","NSMCE3","NSMCE4A","NTHL1","NUDT1","NUDT15","NUDT18","RRM1","RRM2","OGG1","PALB2","PARP2","PARP4","PAXIP1","PER1","POLA1","POLE","POLE2","POLE3","POLE4","POLG","POLH","POLI","POLK","POLN","POLQ","PPP4C","PPP4R1","PPP4R2","PPP4R4","PRPF19","RAD1","RAD17","RAD18","RAD23A","RAD23B","RAD51","RAD51B","RAD51C","RAD51D","RAD52","RAD54B","RAD54L","RAD9A","RBBP8","RBX1","RDM1","RECQL","RECQL4","RECQL5","REV1","REV3L","RIF1","RMI1","RMI2","RNMT","RRM2B","RTEL1","SETMAR","SHFM1","SHPRH","SLX1A","SLX1B","SLX4","SMARCAD1","SMC5","SMC6","SMUG1","SPO11","STRA13","SWSAP1","TCEA1","TCEB1","TCEB2","TCEB3","TDG","TDP1","TELO2","TOP3A","TOP3B","TOPBP1","TP53","TREX1","TREX2","TYMS","UBE2B","UBE2N","UBE2T","UBE2V2","UIMC1","UNG","USP1","UVSSA","WDR48","WRN","XAB2","XPA","XPC","ZSWIM7","PTEN","TDP2","ENDOV","SPRTN","RNF4","SMARCA4","IDH1","SOX4","WEE1","RAD9B","AEN","PLK3","EXO5","CDC5L","BCAS2","PLRG1","YWHAB","YWHAG","YWHAE","CDC25A","CDC25B","CDC25C","BABAM1","BRCC3","TTK","SMARCC1","SWI5","MORF4L1","RNF169","HERC2")
somaticInteractions(maf = Brca_maf, genes = DDR_genes, pvalue=c(0.05,0.1))
DDR_Mutator=c("BRCA2","RNF168","DCLRE1B","FANCI","ALKBH1","POLN","TOP3A","PMS2","XPC","MSH2","MSH4","RAD18","REC8","RFC1","USP7","ATM","CDK12","ERCC4","TOX3","UPF1","LIG1","WDR33","MSH6","ERCC6","UBC")

MDDR_Mutator2<-c("DCLRE1B","SMG1","DOT1L","XPC","ALKBH1","UBC","TNKS1BP1","PMS2","FANCI","UPF1","PARP4","LIG1","EYA1","TOP3A","REC8","WDR33","RNF168","MSH6","MSH2","PSME4","ERCC2","ERCC6","MMS19","RAD51","NIPBL","POLN")

IDDR_Mutator<-c ("DCLRE1B","FANCM","EYA1","POLN","CHRNA4","DOT1L","USP45","HUWE1","BRCA2","SMC1A","MTOR","WDR33","RNF168","ALKBH1","UBC","TIMELESS","MMS19")
#pathway=rep(c("NHEJ","Checkpoint","dealkylation","HR","MMR","Postrepair"),c(2,1,1,2,2,2), stringsAsFactors=FALSE)

DDR_mutators<-c ("MSH2","MSH6","PMS2","LIG1","XPC","USP7","RFC1","UBC","ERCC2","MMS19","ERCC6","DCLRE1B","RNF168","FANCI","RAD51","POLN","ZFYVE26","NIPBL","TONSL","SMC5","SPIDR","ALKBH1","TOP3A","CDK9","EYA1","UBR5","BRCA2","BRCA1")

oncoplot(maf=Brca_maf, genes = IDDR_Mutator, draw_titv = TRUE)

help("oncoplot")

somaticInteractions(maf = Brca_maf, genes = DDR_Mutator, pvalue=c(0.05,0.1))


dev.off()
Test.titv=titv(maf=Brca_maf, plot=FALSE, useSyn = TRUE, pathway=pathway)

plotTiTv(res = Test.titv)
somaticInteractions(maf = Brca_maf, top = 25, pvalue=c(0.05, 0.1))




############
Brca_maf@variants.per.sample$Tumor_Sample_Barcode
Brca_maf=read.maf(maf = "BRCA_req.allSampRemoveExtr.txt", clinicalData = "BRCA_req.allclinic2.txt")
sample_order=sample(Brca_maf@variants.per.sample$Tumor_Sample_Barcode)
help(sample)


### test
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = laml.maf)
annotateTable <- data.frame(Tumor_Sample_Barcode = levels(laml@data[["Tumor_Sample_Barcode"]]),
                            Age = rpois(193, 50))

anno.spl.sort = anno.spl.sort[names(sort(unlist(lapply(anno.spl.sort, ncol)), decreasing = TRUE))]

oncoplot(maf = laml, top = 3, annotationDat = annotateTable, clinicalFeatures = "Age",
         groupAnnotationBySize = FALSE)






##### Draw the co-occurence of muations 
DDR_mutators_BRCACancer<-c ("RNF168","TOP3A","ALKBH1","DCLRE1B","RAD18","PMS2","FANCI","POLN","XPC","REC8","MSH4")

fabcolors = RColorBrewer::brewer.pal(n = 8,name = 'Spectral')
names(fabcolors) = c("M0", "M1", "M2", "M3", "M4", "M5", "M6", "M7")
fabcolors = list(FAB_classification = fabcolors)

somaticInteractions(maf = Brca_maf, genes = DDR_mutators_BRCACancer, pvalue=c(0.05,0.1))
