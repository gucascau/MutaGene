# MutaGene
MutaGene
Genomic instability is a hallmark of most cancers and presumably results from dysfunctions of mutator genes. Yet a catalogue of mutator genes and their roles in cancer is lacking. In this work, we developed a computation method to analyze the TCGA dataset and identified 1534 mutator genes.

The pipeline is to identify potential mutators that can drive the changes of tumor mutation burden.

## Install Prerequesites
GNU utilities
Perl -- 
    module requirement:
    Statistics::Test ; 
    Statistics::R ;  
    List::Util
R

## Tutorial 

### Data collection and quality control

MAF file was collected from Multi-Cancer Mutation calling in the TCGA cohort. 
TCGA clinical (demographic, treatment, pathologic, and survival data) and molecular characterization data (mRNA, protein expression) from 33 cancer types were collected from the cBioPortal (https://www.cbioportal.org/datasets). 
Copy number variation previously defined using ABSOLUTE-estimated purity and ploidy values of each sample were obtained from the PancanAtlas publication page (https://gdc.cancer.gov/about-data/publications/pancanatlas). 
Predicted SNV neoantigen counts and peptide binding data from NetMHCpan were also collected from the PancanAtlas publication page (https://gdc.cancer.gov/about-data/publications/panimmune). 
Cancer-related gene mutations and clinical data from 1662 metastasis patients treated with immune checkpoint inhibitors (ICI) were collected from MSKCC (https://www.cbioportal.org/study/summary?id=tmb_mskcc_2018). 

### Clinical association with mutation load, indel load and CNV burdens

1. Calculate the P value for age, tumor stage, race using a two-way ANOVA
Mutiple_Clinical_status_association_Pvalue.R 
2. Calculate the P value for gender, radiation, smoking using one side Mann-Whitney U 
Two_Clinical_status_association_Pvalue.R

### Sample filteration

1. Remove patients with radiation
2. Remove patients with smoking 
3. Remove patients with extremely high mutation, indel load and CNV burdens.

### Mutator identification

1. calculate P1 and P2 value for each genes based on their existence matrix and mutaiton load information using Mann-Whitney U test

perl MutatorIdentification/MutatorIdentification.pl 
            
           		-i: query files of nonsilent mutations statistics (Mutation_information.txt)
			-g: index files of sample mutation load (mc3)
			-q: request gene lists
			-s: cancer type
           
 2. Extract the mutators from different cancer types : require both P1 and P2 <=0 .05 and larger than 2.
 
### DDRM Analysis

### Immunotherapy Treatment



