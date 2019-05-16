require(data.table)
require(stringr)
require(compare)
library(dplyr)


#plotting the bar plots for the ratio of samples with druggable fusions
rotate_x <- function(data, column_to_plot, labels_vec, rot_angle,y_label,main_label,y_limit) {
  par(mar = c(10, 5, 5, 1))
  colfunc <- colorRampPalette(c("red", "white"))
  cols = colfunc(33)
  plt <- barplot(data[[column_to_plot]], col=cols, xaxt="n",ylim=y_limit,ylab=y_label,las=2,cex.lab=1.5,main=main_label)
  text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.3,0.5), xpd = TRUE, cex=1.1 )
  
}

########### Inputs ##############
# smachete_fusions_file = "/scratch/PI/horence/Roozbeh/TCGA_Project/consolidated_files/smachete_fusions_filtered.txt"     #the list of fusions called by smachete for all TCGA cancers
# knife_fusions_file = "/scratch/PI/horence/Roozbeh/TCGA_Project/consolidated_files/knife_fusions_filtered.txt"           #the list of fusions called by knife for all TCGA cancers
# TCGA_sample_run_file = "/scratch/PI/horence/Roozbeh/TCGA_Project/utility_files/TCGA_sample_run_stats.txt"
# mean_fusion_number_per_cancer_file = "/scratch/PI/horence/Roozbeh/TCGA_Project/consolidated_files/mean_detected_fusion_per_cancer.txt"
# MSK_onco_tsg_file = "/scratch/PI/horence/Roozbeh/TCGA_Project/systems_biology_files/MSK_oncogenes_TSG_List.txt"         #the list of MSK oncogenes and tumor suppressor genes
# gene_drug_interaction_file = "/scratch/PI/horence/Roozbeh/TCGA_Project/systems_biology_files/interactions.tsv"
# dGene_file = "/scratch/PI/horence/Roozbeh/TCGA_Project/systems_biology_files/dGene.tsv"
# msk_file = "/scratch/PI/horence/Roozbeh/TCGA_Project/systems_biology_files/allActionableVariants.txt"
# kinase_genes_file = "/scratch/PI/horence/Roozbeh/TCGA_Project/systems_biology_files/kinase_genes_new.csv"
# COSMIC_genes_file = "/scratch/PI/horence/Roozbeh/TCGA_Project/systems_biology_files/cancer_gene_census.csv"
################################


########### Inputs (local PC) ##############
smachete_fusions_file = "G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\consolidated_files\\smachete_fusions_filtered.txt"     #the list of fusions called by smachete for all TCGA cancers
knife_fusions_file = "G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\consolidated_files\\knife_fusions_filtered.txt"           #the list of fusions called by knife for all TCGA cancers
TCGA_sample_run_file = "G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\utility_files\\TCGA_sample_run_stats.txt"
mean_fusion_number_per_cancer_file = "G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\consolidated_files\\mean_detected_fusion_per_cancer.txt"
MSK_onco_tsg_file = "G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\systems_biology_files\\MSK_oncogenes_TSG_List.txt"         #the list of MSK oncogenes and tumor suppressor genes
gene_drug_interaction_file = "G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\systems_biology_files\\interactions.tsv"
dGene_file = "G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\systems_biology_files\\dGene.tsv"
msk_file = "G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\systems_biology_files\\allActionableVariants.txt"
kinase_genes_file = "G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\systems_biology_files\\kinase_genes_new.csv"
COSMIC_genes_file = "G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\systems_biology_files\\cancer_gene_census.csv"
#############################################


######read in input files####################
smachete_fusions = fread(smachete_fusions_file,sep = "\t",header = TRUE)
knife_fusions = fread(knife_fusions_file,sep="\t",header = TRUE)
TCGA_sample_run = fread(TCGA_sample_run_file,header = TRUE,sep = "\t",fill = TRUE)
TCGA_sample_run = TCGA_sample_run[!(Project%like%"GTEx")]
TCGA_sample_run = TCGA_sample_run[!(Project%like%"ALL")]
mean_fusion_number_per_cancer = fread(mean_fusion_number_per_cancer_file,header = TRUE,sep = "\t")
MSK_onco_tsg = fread(MSK_onco_tsg_file,sep = "\t",header = TRUE)
gene_drug_interaction = fread(gene_drug_interaction_file,sep="\t",header = TRUE)
gene_drug_interaction = gene_drug_interaction[interaction_claim_source%like%"Cancer"|interaction_types=="inhibitor"]
dGene = fread(dGene_file,sep = "\t",header = TRUE)
msk = fread(msk_file,sep = "\t",header = TRUE)
kinase_genes = fread(kinase_genes_file,sep = ",",header = TRUE)
cosmic_genes = fread(COSMIC_genes_file,sep = ",",header = TRUE)
#############################################



knife_fusions[,gene1:=strsplit(junction,split = ":")[[1]][2],by = 1:nrow(knife_fusions)]
knife_fusions[,gene2:=strsplit(junction,split = ":")[[1]][4],by = 1:nrow(knife_fusions)]

setnames(smachete_fusions,old = "X.Junction",new = "junction")
#all_fusions = smachete_fusions[,list(TCGA_Project,sample_name,fusion,gene1,gene2)]
all_fusions = rbind(knife_fusions[,list(TCGA_Project,sample_name,fusion,gene1,gene2)],smachete_fusions[,list(TCGA_Project,sample_name,fusion,gene1,gene2)])
all_fusions = unique(all_fusions)   #all fusions called by knife and smachete for all TCGA cancers

#########genes found in the gene-drug interaction list############
all_fusions[,is.drug_interaction:=0]
all_fusions[gene1%in%gene_drug_interaction$gene_name | gene2%in%gene_drug_interaction$gene_name | gene1%in%gene_drug_interaction$gene_claim_name | gene2%in%gene_drug_interaction$gene_claim_name,is.drug_interaction:=1]


########genes found in the dGene list as potentially druggable fusions#####################3
all_fusions[,is.druggable:=0]
all_fusions[gene1%in%dGene$Symbol | gene2%in%dGene$Symbol | gene1%in%dGene$Synonyms | gene2%in%dGene$Synonyms,is.druggable:=1]

##############msk fusions######
all_fusions[,is.msk:=0]
all_fusions[gene1%in%msk$Gene | gene2%in%msk$Gene,is.msk:=1]

####kinase fusions##############
all_fusions[,is.kinase:=0]
all_fusions[gene1%in%kinase_genes$Name | gene1%in%kinase_genes$Old_Name| gene1%in%kinase_genes$Entrez_Symbol |gene1%in%kinase_genes$Entrez_Synonyms | gene2%in%kinase_genes$Name | gene2%in%kinase_genes$Old_Name| gene2%in%kinase_genes$Entrez_Symbol |gene2%in%kinase_genes$Entrez_Synonyms,is.kinase:=1]

#####cosmic fusions###############
all_fusions[,cosmic:=0]
all_cosmic_synonyms=unlist((strsplit(cosmic_genes$Synonyms,split = ",")))
all_fusions[gene1%in%cosmic_genes$`Gene Symbol` | gene1%in%all_cosmic_synonyms |gene2%in%cosmic_genes$`Gene Symbol` | gene2%in%all_cosmic_synonyms,cosmic:=1]


################################
##########Below we compute the number of samples within each cancer that have the specific kind of fusion being considered 

##ratio of samples with msk fusions
num_msk_sample_per_cancer=all_fusions[is.msk==1,length(unique(sample_name)),by=TCGA_Project] 
num_msk_sample_per_cancer=merge(num_msk_sample_per_cancer,TCGA_sample_run,all.x = TRUE,all.y=TRUE,by.x = "TCGA_Project",by.y = "Project")
num_msk_sample_per_cancer[,msk_ratio:=V1/machete_tumore_samples,by=1:nrow(num_msk_sample_per_cancer)]
num_msk_sample_per_cancer[is.na(num_msk_sample_per_cancer)]=0


##ratio of samples with potentially druggable fusions
num_druggable_sample_per_cancer=all_fusions[is.druggable==1,length(unique(sample_name)),by=TCGA_Project]
num_druggable_sample_per_cancer=merge(num_druggable_sample_per_cancer,TCGA_sample_run,all.x = TRUE,all.y=TRUE,by.x = "TCGA_Project",by.y = "Project")
num_druggable_sample_per_cancer[,druggable_ratio:=V1/machete_tumore_samples,by=1:nrow(num_druggable_sample_per_cancer)]
num_druggable_sample_per_cancer[is.na(num_druggable_sample_per_cancer)]=0

##ratio of samples with drug_inetraction fusions
num_drug_interaction_sample_per_cancer=all_fusions[is.drug_interaction==1,length(unique(sample_name)),by=TCGA_Project]
num_drug_interaction_sample_per_cancer=merge(num_drug_interaction_sample_per_cancer,TCGA_sample_run,all.x = TRUE,all.y=TRUE,by.x = "TCGA_Project",by.y = "Project")
num_drug_interaction_sample_per_cancer[,drug_interaction_ratio:=V1/machete_tumore_samples,by=1:nrow(num_drug_interaction_sample_per_cancer)]
num_drug_interaction_sample_per_cancer[is.na(num_drug_interaction_sample_per_cancer)]=0

##ratio of samples with kinase fusions
num_kinase_sample_per_cancer = all_fusions[is.kinase == 1,length(unique(sample_name)),by = TCGA_Project]
num_kinase_sample_per_cancer = merge(num_kinase_sample_per_cancer,TCGA_sample_run,all.x = TRUE,all.y = TRUE,by.x = "TCGA_Project",by.y = "Project")
num_kinase_sample_per_cancer[,kinase_ratio := V1/machete_tumore_samples,by = 1:nrow(num_kinase_sample_per_cancer)]
num_kinase_sample_per_cancer[is.na(num_kinase_sample_per_cancer)] = 0
num_kinase_sample_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(num_kinase_sample_per_cancer)]
num_kinase_sample_per_cancer = setorder(num_kinase_sample_per_cancer, -kinase_ratio)

p_kinase_fusion = 0.05556942  #the probability of having kinase fusions based on random pairing of genes: 1-(1-620/22000)^2  
num_kinase_sample_per_cancer = merge(num_kinase_sample_per_cancer,mean_fusion_number_per_cancer,by.x = "TCGA_Project",by.y="TCGA_Project")
num_kinase_sample_per_cancer[,null_frac_kinase_samples:=1-(1-p_kinase_fusion)^mean,by=1:nrow(num_kinase_sample_per_cancer)]
num_kinase_sample_per_cancer = setorder(num_kinase_sample_per_cancer, -kinase_ratio)
num_kinase_sample_per_cancer[,p_value:= 1-pbinom(q=V1,size=machete_tumore_samples,prob=null_frac_kinase_samples),by=1:nrow(num_kinase_sample_per_cancer)]

par(mar = c(10, 5, 5, 1))
colfunc <- colorRampPalette(c("red", "white"))
cols = colfunc(33)
plt <- barplot(num_kinase_sample_per_cancer[['kinase_ratio']], col=cols, xaxt="n",ylim=c(0,0.45),ylab="Fraction of samples with kinase fusions",las=2,cex.lab=1.5,main="")
text(plt, par("usr")[3], labels = num_kinase_sample_per_cancer$TCGA_Project, srt = 90, adj = c(1.3,0.5), xpd = TRUE, cex=1.1 )
points(x = plt, y = num_kinase_sample_per_cancer$null_frac_kinase_samples,lwd=2,pch = 19)
points(x = plt, y = num_kinase_sample_per_cancer$kinase_ratio,lwd=2,pch = 19)


###here I adjust the positions for adding the project names on the graph
x_pos= num_kinase_sample_per_cancer$null_frac_kinase_samples
y_pos=num_kinase_sample_per_cancer$kinase_ratio

y_pos[4]=y_pos[4]+0.038       #STAD

y_pos[5]=y_pos[5]+0.038       #UCS

x_pos[6]=x_pos[6]-0.012       #BRCA
y_pos[6]=y_pos[6]+0.018

y_pos[8]=y_pos[8]+0.04        #LUAD

x_pos[9]=x_pos[9]+0.013       #BLCA
y_pos[9]=y_pos[9]+0.021 

y_pos[10]=y_pos[10]+0.038     #LUSC

y_pos[11]=y_pos[11]+0.038     #CHOL

x_pos[13]=x_pos[13]+0.013     #PRAD
y_pos[13]=y_pos[13]+0.02 

y_pos[15]=y_pos[15]+0.038     #HNSC

x_pos[18]=x_pos[18]+0.012     #ACC
y_pos[18]=y_pos[18]+0.014

x_pos[19]=x_pos[19]+0.012     #LIHC
y_pos[19]=y_pos[19]+0.014

x_pos[20]=x_pos[20]+0.012     #MESO
y_pos[20]=y_pos[20]+0.039

x_pos[21]=x_pos[21]+0.012     #LGG
y_pos[21]=y_pos[21]+0.014

x_pos[22]=x_pos[22]+0.012     #READ
y_pos[22]=y_pos[22]+0.014

y_pos[23]=y_pos[23]+0.038     #PAAD

x_pos[24]=x_pos[24]+0.001     #COAD
y_pos[24]=y_pos[24]+0.03      

#x_pos[25]=x_pos[25]+0.013     #LAML
y_pos[25]=y_pos[25]+0.03

x_pos[26]=x_pos[26]+0.012     #KIRP
y_pos[26]=y_pos[26]+0.021

x_pos[27]=x_pos[27]-0.002    #DLBC
y_pos[27]=y_pos[27]+0.03     


x_pos[28]=x_pos[28]-0.002   
y_pos[28]=y_pos[28]-0.005   #UVM

x_pos[29]=x_pos[29]+0.023     #KIRC
y_pos[29]=y_pos[29]+0.017

x_pos[30]=x_pos[30]-0.01     #PCPG
y_pos[30]=y_pos[30]+0.01

x_pos[31]=x_pos[31]-0.012     #TGCT
y_pos[31]=y_pos[31]+0.02

x_pos[32]=x_pos[32]+0.013       #THYM
y_pos[32]=y_pos[32]+0.019

x_pos[33]=x_pos[33]+0.013       #KICH
y_pos[33]=y_pos[33]+0.021 


plot(num_kinase_sample_per_cancer$null_frac_kinase_samples, num_kinase_sample_per_cancer$kinase_ratio, xlim=c(0.0105,0.4),ylim=c(0.001,0.4), xlab="Expected fraction of samples with kinase fusions" ,ylab="Fraction of samples with detected kinase fusions", pch=19,cex=1.4,cex.lab=1.1,col = "blue")
abline(a=0,b=1,cex=20,col=1,lwd=2)
text(0.35,0.37, "y=x", col = 1,cex=1.6)

text(x = x_pos,y = y_pos,labels= num_kinase_sample_per_cancer$TCGA_Project, cex = 0.7,pos = 1)
##ratio of samples with cosmic fusions


num_cosmic_sample_per_cancer = all_fusions[cosmic == 1,length(unique(sample_name)),by = TCGA_Project]
num_cosmic_sample_per_cancer = merge(num_cosmic_sample_per_cancer,TCGA_sample_run,all.x = TRUE,all.y = TRUE,by.x = "TCGA_Project",by.y = "Project")
num_cosmic_sample_per_cancer[,cosmic_ratio := V1/machete_tumore_samples,by = 1:nrow(num_cosmic_sample_per_cancer)]
num_cosmic_sample_per_cancer[is.na(num_cosmic_sample_per_cancer)] = 0
num_cosmic_sample_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(num_cosmic_sample_per_cancer)]


p_cosmic_fusion = 0.0643   #the probability of having cosmic fusions based on random pairing of genes: 1-(1-719/22000)^2
num_cosmic_sample_per_cancer = merge(num_cosmic_sample_per_cancer,mean_fusion_number_per_cancer,by.x = "TCGA_Project",by.y="TCGA_Project")
num_cosmic_sample_per_cancer[,null_frac_cosmic_samples:=1-(1-p_cosmic_fusion)^mean,by=1:nrow(num_cosmic_sample_per_cancer)]
num_cosmic_sample_per_cancer[,p_value:= pbinom(q=V1-1,size=machete_tumore_samples,prob=null_frac_cosmic_samples,lower.tail = FALSE),by=1:nrow(num_cosmic_sample_per_cancer)]
num_cosmic_sample_per_cancer = setorder(num_cosmic_sample_per_cancer, -cosmic_ratio)

par(mar = c(10, 5, 5, 1))
plt <- barplot(num_cosmic_sample_per_cancer[['cosmic_ratio']], xaxt="n",ylim=c(0,0.7),col="white",border=NA,ylab="Fraction of samples with COSMIC fusions",las=2,cex.lab=1.3,main="")
text(plt, par("usr")[3], labels = num_cosmic_sample_per_cancer$TCGA_Project, srt = 90, adj = c(1.3,0.5), xpd = TRUE, cex=1.1 )
points(x = plt, y = num_cosmic_sample_per_cancer$null_frac_cosmic_samples,lwd=2,cex=1.3, pch = 19,col = "indianred1")
points(x = plt, y = num_cosmic_sample_per_cancer$cosmic_ratio,lwd=2,pch = 19,cex=1.3,col = "red3")
legend(x = 30, y = 0.6, legend = c("Detected", "Expected"), pch = 19:19, col = c("red3", "indianred1"),cex = 1.2)


##ratio of samples with detected fusions
ratio_samples_with_fusions = all_fusions[,length(unique(sample_name)),by=TCGA_Project]
ratio_samples_with_fusions = merge(ratio_samples_with_fusions,TCGA_sample_run,all.x = TRUE,all.y = TRUE,by.x = "TCGA_Project",by.y = "Project")
ratio_samples_with_fusions[,sample_with_fusion_ratio := V1/machete_tumore_samples,by = 1:nrow(ratio_samples_with_fusions)]


########################################################################
############the fraction of fusions with specific category for each caner

total_number_detected_fusions_per_cancer = all_fusions[,.N,by = TCGA_Project]
setnames(total_number_detected_fusions_per_cancer,"N","total_number_detected_fusions")
total_number_kinase_fusions_per_cancer = all_fusions[is.kinase==1,.N,by = TCGA_Project]
setnames(total_number_kinase_fusions_per_cancer,"N","total_number_kinase_fusions")
total_number_kinase_fusions_per_cancer = merge(total_number_kinase_fusions_per_cancer,total_number_detected_fusions_per_cancer,by.x="TCGA_Project",by.y="TCGA_Project")
total_number_kinase_fusions_per_cancer[,ratio_kinase_fusions:=total_number_kinase_fusions/total_number_detected_fusions,by=1:nrow(total_number_kinase_fusions_per_cancer)]
abline(h = 0.0439, col="black", lwd=3, lty=2)

total_number_kinase_fusions_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(total_number_kinase_fusions_per_cancer)]
total_number_kinase_fusions_per_cancer = setorder(total_number_kinase_fusions_per_cancer, -ratio_kinase_fusions)
rotate_x(total_number_kinase_fusions_per_cancer, 'ratio_kinase_fusions', total_number_kinase_fusions_per_cancer$TCGA_Project, 90,"Fraction of kinase fusions","",c(0,0.26))


######################################################################


total_number_kinase_cosmic_fusions_per_cancer = all_fusions[is.kinase==1|cosmic==1,.N,by = TCGA_Project]
setnames(total_number_kinase_cosmic_fusions_per_cancer,"N","total_number_kinase/cosmic_fusions")
total_number_kinase_cosmic_fusions_per_cancer = merge(total_number_kinase_cosmic_fusions_per_cancer,total_number_detected_fusions_per_cancer,by.x="TCGA_Project",by.y="TCGA_Project")
total_number_kinase_cosmic_fusions_per_cancer[,ratio_kinase_cosmic_fusions:=`total_number_kinase/cosmic_fusions`/total_number_detected_fusions,by=1:nrow(total_number_kinase_cosmic_fusions_per_cancer)]


total_number_kinase_cosmic_fusions_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(total_number_kinase_cosmic_fusions_per_cancer)]
total_number_kinase_cosmic_fusions_per_cancer = setorder(total_number_kinase_cosmic_fusions_per_cancer, -ratio_kinase_cosmic_fusions)
rotate_x(total_number_kinase_cosmic_fusions_per_cancer, 'ratio_kinase_cosmic_fusions', total_number_kinase_cosmic_fusions_per_cancer$TCGA_Project, 90,"Fraction of fusions involving either a kinase or cosmic fusion","",c(0,0.6))
abline(h = 0.0887, col="black", lwd=3, lty=2)



#########################################################
#######plotting the results for each fusion type

num_druggable_sample_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(num_druggable_sample_per_cancer)]
num_druggable_sample_per_cancer=setorder(num_druggable_sample_per_cancer, -druggable_ratio)
rotate_x(num_druggable_sample_per_cancer, 'druggable_ratio', num_druggable_sample_per_cancer$TCGA_Project, 90,"Fraction of samples with druggable fusions","Druggable Fusions Across 33 TCGA Cancers",c(0,0.6))

num_msk_sample_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(num_msk_sample_per_cancer)]
num_msk_sample_per_cancer=setorder(num_msk_sample_per_cancer, -msk_ratio)
rotate_x(num_msk_sample_per_cancer, 'msk_ratio', num_msk_sample_per_cancer$TCGA_Project, 90,"Fraction of samples with druggable fusions","Druggable Fusions Across 33 TCGA Cancers",c(0,0.2))

num_drug_interaction_sample_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(num_drug_interaction_sample_per_cancer)]
num_drug_interaction_sample_per_cancer=setorder(num_drug_interaction_sample_per_cancer, -drug_interaction_ratio)
rotate_x(num_drug_interaction_sample_per_cancer, 'drug_interaction_ratio', num_drug_interaction_sample_per_cancer$TCGA_Project, 90,"Fraction of samples with druggable fusions","Druggable Fusions Across 33 TCGA Cancers",c(0,0.6))

num_kinase_sample_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(num_kinase_sample_per_cancer)]
num_kinase_sample_per_cancer=setorder(num_kinase_sample_per_cancer, -kinase_ratio)
rotate_x(num_kinase_sample_per_cancer, 'kinase_ratio', num_kinase_sample_per_cancer$TCGA_Project, 90,"Fraction of samples with kinase fusions","Kinase Fusions Across 33 TCGA Cancers",c(0,0.4))

num_cosmic_sample_per_cancer[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by=1:nrow(num_cosmic_sample_per_cancer)]
num_cosmic_sample_per_cancer=setorder(num_cosmic_sample_per_cancer, -cosmic_ratio)
rotate_x(num_cosmic_sample_per_cancer, 'cosmic_ratio', num_cosmic_sample_per_cancer$TCGA_Project, 90,"Fraction of samples with fusions involving COSMIC genes","Ratio of Samples with detected Fusions involving COSMIC genes",c(0,0.6))


ratio_samples_with_fusions[,TCGA_Project:=gsub("TCGA_","",TCGA_Project),by = 1:nrow(ratio_samples_with_fusions)]
ratio_samples_with_fusions = setorder(ratio_samples_with_fusions, -sample_with_fusion_ratio)
rotate_x(ratio_samples_with_fusions, 'sample_with_fusion_ratio', ratio_samples_with_fusions$TCGA_Project, 90,"Fraction of samples containing fusions","",c(0,1))
#######################################################################

####################################
#############identifying fusions involving oncogenes and tumor suppressors
onco_genes = MSK_onco_tsg[`OncoKB Oncogene` == "Yes"]$`Hugo Symbol`  #genes identified as oncogenes in the msk file
tumor_suppressor_genes = MSK_onco_tsg[`OncoKB TSG` == "Yes"]$`Hugo Symbol`  #genes identified as tumor suppressor genes in the msk file

onco_genes_dt = data.table(onco_genes)
tumor_suppressor_genes_dt = data.table(tumor_suppressor_genes)

all_fusions[,is.five_prime_oncogene:=0]           #5' gene in oncoegene
all_fusions[,is.three_prime_oncogene:=0]          #3' gene is oncogene 
all_fusions[,is.five_prime_tumor_suppressor:=0]   #5' gene is tumor suppressor
all_fusions[,is.three_prime_tumor_suppressor:=0]  #3' gene is tumor suppressor
all_fusions[gene1%in%onco_genes_dt$onco_genes,is.five_prime_oncogene:=1]
all_fusions[gene2%in%onco_genes_dt$onco_genes,is.three_prime_oncogene:=1]
all_fusions[gene1%in%tumor_suppressor_genes_dt$tumor_suppressor_genes,is.five_prime_tumor_suppressor:=1]
all_fusions[gene2%in%tumor_suppressor_genes_dt$tumor_suppressor_genes,is.three_prime_tumor_suppressor:=1]
##################################################################################



write.table(all_fusions,"G:\\My Drive\\Postdoc_Research\\Projects\\TCGA_sMACHETE\\consolidated_files\\oncogenic_TSG_fusions.txt", row.names = FALSE, quote = FALSE, sep = "\t")