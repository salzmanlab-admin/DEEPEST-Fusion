require(data.table)
require(stringr)
require(compare)
library(dplyr)

#All input files are available in DEEPEST-Fusion Github repository, custom scripts folder at: https://github.com/salzmanlab/DEEPEST-Fusion/tree/master/custom_scripts/files 


########### Input files ##############
DEEPEST_fusions_with_ID_file = "DEEPEST-Fusion/custom_scripts/files/DEEPEST_fusions_with_ID.txt"
TCGA_sample_run_file = "DEEPEST-Fusion/custom_scripts/files/TCGA_sample_run_stats.txt"
mean_fusion_number_per_cancer_file = "DEEPEST-Fusion/custom_scripts/files/mean_detected_fusion_per_cancer.txt"
kinase_genes_file = "DEEPEST-Fusion/custom_scripts/files/kinase_genes_new.csv"
COSMIC_genes_file = "DEEPEST-Fusion/custom_scripts/files/cancer_gene_census.csv"
#############################################


###### read in input files ####################
DEEPEST_fusions_with_ID = fread(DEEPEST_fusions_with_ID_file,sep="\t",header = TRUE)
TCGA_sample_run = fread(TCGA_sample_run_file,header = TRUE,sep = "\t",fill = TRUE)
TCGA_sample_run = TCGA_sample_run[!(Project%like%"GTEx")]
TCGA_sample_run = TCGA_sample_run[!(Project%like%"ALL")]
mean_fusion_number_per_cancer = fread(mean_fusion_number_per_cancer_file,header = TRUE,sep = "\t")
kinase_genes = fread(kinase_genes_file,sep = ",",header = TRUE)
cosmic_genes = fread(COSMIC_genes_file,sep = ",",header = TRUE)
#############################################


#### finding kinase fusions  ##############
DEEPEST_fusions_with_ID = unique(DEEPEST_fusions_with_ID)
DEEPEST_fusions_with_ID[,is.kinase:=0]
DEEPEST_fusions_with_ID[gene1%in%kinase_genes$Name | gene1%in%kinase_genes$Old_Name| gene1%in%kinase_genes$Entrez_Symbol |gene1%in%kinase_genes$Entrez_Synonyms | gene2%in%kinase_genes$Name | gene2%in%kinase_genes$Old_Name| gene2%in%kinase_genes$Entrez_Symbol |gene2%in%kinase_genes$Entrez_Synonyms,is.kinase:=1]

##### finding cosmic fusion s###############
DEEPEST_fusions_with_ID[,cosmic:=0]
all_cosmic_synonyms=unlist((strsplit(cosmic_genes$Synonyms,split = ",")))
DEEPEST_fusions_with_ID[gene1%in%cosmic_genes$`Gene Symbol` | gene1%in%all_cosmic_synonyms |gene2%in%cosmic_genes$`Gene Symbol` | gene2%in%all_cosmic_synonyms,cosmic:=1]


## ratio of samples with kinase fusions
num_kinase_sample_per_cancer = DEEPEST_fusions_with_ID[is.kinase == 1,length(unique(sample_name)),by = TCGA_Project]
num_kinase_sample_per_cancer = merge(num_kinase_sample_per_cancer,TCGA_sample_run,all.x = TRUE,all.y = TRUE,by.x = "TCGA_Project",by.y = "Project")
num_kinase_sample_per_cancer[,kinase_ratio := V1/machete_tumore_samples,by = 1:nrow(num_kinase_sample_per_cancer)]
num_kinase_sample_per_cancer[is.na(num_kinase_sample_per_cancer)] = 0
num_kinase_sample_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(num_kinase_sample_per_cancer)]
num_kinase_sample_per_cancer = setorder(num_kinase_sample_per_cancer, -kinase_ratio)

# computing the p-value for the fraction of samples with kinase fusions for each TCGA tumor type
p_kinase_fusion = 0.05556942  #the probability of having kinase fusions based on random pairing of genes: 1-(1-620/22000)^2  
num_kinase_sample_per_cancer = merge(num_kinase_sample_per_cancer,mean_fusion_number_per_cancer,by.x = "TCGA_Project",by.y="TCGA_Project")
num_kinase_sample_per_cancer[,null_frac_kinase_samples:=1-(1-p_kinase_fusion)^mean,by=1:nrow(num_kinase_sample_per_cancer)]
num_kinase_sample_per_cancer = setorder(num_kinase_sample_per_cancer, -kinase_ratio)
num_kinase_sample_per_cancer[,p_value:= 1-pbinom(q=V1,size=machete_tumore_samples,prob=null_frac_kinase_samples),by=1:nrow(num_kinase_sample_per_cancer)]

#### plotting kinase statistical analysis results
plot(num_kinase_sample_per_cancer$null_frac_kinase_samples, num_kinase_sample_per_cancer$kinase_ratio, xlim=c(0.0105,0.4),ylim=c(0.001,0.4), xlab="Expected fraction of samples with kinase fusions" ,ylab="Fraction of samples with detected kinase fusions", pch=19,cex=1.4,cex.lab=1.1,col = "blue")
abline(a=0,b=1,cex=20,col=1,lwd=2)
text(0.35,0.37, "y=x", col = 1,cex=1.6)
text(x = x_pos,y = y_pos,labels= num_kinase_sample_per_cancer$TCGA_Project, cex = 0.7,pos = 1)
################


## Now ratio of samples with cosmic fusions
num_cosmic_sample_per_cancer = DEEPEST_fusions_with_ID[cosmic == 1,length(unique(sample_name)),by = TCGA_Project]
num_cosmic_sample_per_cancer = merge(num_cosmic_sample_per_cancer,TCGA_sample_run,all.x = TRUE,all.y = TRUE,by.x = "TCGA_Project",by.y = "Project")
num_cosmic_sample_per_cancer[,cosmic_ratio := V1/machete_tumore_samples,by = 1:nrow(num_cosmic_sample_per_cancer)]
num_cosmic_sample_per_cancer[is.na(num_cosmic_sample_per_cancer)] = 0
num_cosmic_sample_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(num_cosmic_sample_per_cancer)]

# computing the p-value for the fraction of samples with cosmic fusions for each TCGA tumor type
p_cosmic_fusion = 0.0643   #the probability of having cosmic fusions based on random pairing of genes: 1-(1-719/22000)^2
num_cosmic_sample_per_cancer = merge(num_cosmic_sample_per_cancer,mean_fusion_number_per_cancer,by.x = "TCGA_Project",by.y="TCGA_Project")
num_cosmic_sample_per_cancer[,null_frac_cosmic_samples:=1-(1-p_cosmic_fusion)^mean,by=1:nrow(num_cosmic_sample_per_cancer)]
num_cosmic_sample_per_cancer[,p_value:= pbinom(q=V1-1,size=machete_tumore_samples,prob=null_frac_cosmic_samples,lower.tail = FALSE),by=1:nrow(num_cosmic_sample_per_cancer)]
num_cosmic_sample_per_cancer = setorder(num_cosmic_sample_per_cancer, -cosmic_ratio)

par(mar = c(10, 5, 5, 1))
points(x = plt, y = num_cosmic_sample_per_cancer$null_frac_cosmic_samples,lwd=2,cex=1.3, pch = 19,col = "indianred1")
points(x = plt, y = num_cosmic_sample_per_cancer$cosmic_ratio,lwd=2,pch = 19,cex=1.3,col = "red3")
legend(x = 30, y = 0.6, legend = c("Detected", "Expected"), pch = 19:19, col = c("red3", "indianred1"),cex = 1.2)
