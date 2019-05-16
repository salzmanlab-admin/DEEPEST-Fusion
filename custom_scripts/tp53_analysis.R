require(data.table)
require(stringr)
require(compare)
library(dplyr)
library(ggplot2)


#All input files are available in DEEPEST-Fusion Github repository, custom scripts folder at: https://github.com/salzmanlab/DEEPEST-Fusion/tree/master/custom_scripts/files 

################### Inputs files #####################
cell_fusions_file = "DEEPEST-Fusion/custom_scripts/files/cell_fusions.csv"
P53_freq_file = "DEEPEST-Fusion/custom_scripts/files/TP53_mutation.txt"
average_num_deepest_fusions_file = "DEEPEST-Fusion/custom_scripts/files/mean_detected_fusion_per_cancer.txt"
cell_paper_samples_file = "DEEPEST-Fusion/custom_scripts/files/samples_used_in_cell_paper.csv"
DEEPEST_fusions_with_ID_file = "DEEPEST-Fusion/custom_scripts/files/DEEPEST_fusions_with_ID.txt"
DEEPEST_samples_file = "DEEPEST-Fusion/custom_scripts/files/samples_used_in_DEEPEST.txt"
P53_mutated_samples_file = "DEEPEST-Fusion/custom_scripts/files/TCGA_samples_with_tp53_mutation.tsv"
tumorfusions_file = "DEEPEST-Fusion/custom_scripts/files/tumorfusions.csv"
tumorfusions_samples_file = "DEEPEST-Fusion/custom_scripts/files/samples_used_in_PRADA.csv"
###########################################################



#######  read in input files  #########################
cell_fusions = fread(cell_fusions_knife_file,sep=",",header = TRUE)
samples_used_in_cell_paper = fread(cell_paper_samples_file, header = TRUE, sep = ",")
P53_freq = fread(P53_freq_file,header = TRUE,sep="\t")
average_num_DEEPEST_fusions = fread(average_num_DEEPEST_fusions_file,header = TRUE,sep="\t")
samples_used_in_cell_paper = fread(cell_paper_samples_file, header = TRUE, sep = ",")
samples_used_in_DEEPEST = fread(DEEPEST_samples_file,header = TRUE,sep = "\t" )
DEEPEST_fusions_with_ID = fread(DEEPEST_fusions_with_ID_file,sep="\t",header = TRUE)
P53_mutated_samples = fread(P53_mutated_samples_file,sep="\t",header=TRUE)
tumorfusions = fread(tumorfusions_knife_file, header = TRUE, sep = ",") 
samples_used_in_tumorfusions = fread(tumorfusions_samples_file, header = TRUE, sep = ",")
##############################################


###compute correlations for DEEPEST
average_num_DEEPEST_fusions_vs_p53_freq = merge(average_num_DEEPEST_fusions,P53_freq,by.x = "TCGA_Project",by.y = "cancer")
DEEPEST_pearson = cor(average_num_DEEPEST_fusions_vs_p53_freq$mean, average_num_DEEPEST_fusions_vs_p53_freq$mutation_frequency, method = "pearson")
DEEPEST_spearman = cor(average_num_DEEPEST_fusions_vs_p53_freq$mean, average_num_DEEPEST_fusions_vs_p53_freq$mutation_frequency, method = "spearman")


#####compute correlations for cell
tumor_samples_used_in_cell_paper = samples_used_in_cell_paper[!((Sample%like%"-11A-") | (Sample%like%"-11B-") | (Sample%like%"-11C-"))]
num_cell_samples = tumor_samples_used_in_cell_paper[,length(Sample),by = Cancer]
num_tot_cell_fusions_per_cancer = cell_fusions[,length(Fusion),by=Cancer]
num_tot_cell_fusions_per_cancer = merge(num_tot_cell_fusions_per_cancer,num_cell_samples,by.x="Cancer",by.y = "Cancer",all.x = TRUE,all.y=TRUE)
average_num_cell_fusions = num_tot_cell_fusions_per_cancer[,V1.x/V1.y,by=Cancer]   #average number of fusions reported in the cell paper per cancer
average_num_cell_fusions_vs_p53_freq = merge(average_num_cell_fusions,P53_freq,by.x = "Cancer",by.y = "cancer")

cell_pearson = cor(average_num_cell_fusions_vs_p53_freq$V1,  average_num_cell_fusions_vs_p53_freq$mutation_frequency, method = "pearson")
cell_spearman = cor(average_num_cell_fusions_vs_p53_freq$V1,  average_num_cell_fusions_vs_p53_freq$mutation_frequency, method = "spearman")


#####compute correlations for prada
tumor_samples_used_in_tumorfusions = samples_used_in_tumorfusions[!((barcode%like%"-11A-") | (barcode%like%"-11B-") | (barcode%like%"-11C-"))]
num_tumorfusions_samples = tumor_samples_used_in_tumorfusions[,.N,by = Disease]
num_tot_tumorfusions_per_cancer = tumorfusions[,.N,by=Cancer]
num_tot_tumorfusions_per_cancer = merge(num_tot_tumorfusions_per_cancer,num_tumorfusions_samples,by.x="Cancer",by.y = "Disease",all.x = TRUE,all.y=TRUE)
average_num_tumorfusions = num_tot_tumorfusions_per_cancer[,N.x/N.y,by=Cancer]   #average number of fusions reported in the cell paper per cancer
average_num_tumorfusions_vs_p53_freq = merge(average_num_tumorfusions,P53_freq,by.x = "Cancer",by.y = "cancer")

tumorfusions_pearson = cor(average_num_tumorfusions_vs_p53_freq$V1,  average_num_tumorfusions_vs_p53_freq$mutation_frequency, method = "pearson")
tumorfusions_spearman = cor(average_num_tumorfusions_vs_p53_freq$V1,  average_num_tumorfusions_vs_p53_freq$mutation_frequency, method = "spearman")

################################
#comarping the abudnance of fusions in samples w and w/o TP53  mutation


samples_used_in_DEEPEST[,sample_name_new:=substr(sample_name,1,nchar(sample_name)-1),by=1:nrow(samples_used_in_DEEPEST)]
samples_used_in_DEEPEST[,is.TP53_mutated:=0]
samples_used_in_DEEPEST[sample_name_new%in%P53_mutated_samples$`Sample ID`,is.TP53_mutated:=1]
samples_used_in_DEEPEST = samples_used_in_DEEPEST[!((sample_name%like%"11A") | (sample_name%like%"11B") | (sample_name%like%"11C"))]

DEEPEST_fusions_with_ID[,sample_fusion:=paste(sample_name,fusion_unified,sep="--"),by=1:nrow(DEEPEST_fusions_with_ID)]
DEEPEST_fusions_with_ID = DEEPEST_fusions_with_ID[!(duplicated(sample_fusion))]
num_fusion_per_sample = DEEPEST_fusions_with_ID[,.N,by=sample_name]
samples_num_fusion_tp53 = merge(samples_used_in_DEEPEST,num_fusion_per_sample,by.x = "sample_name",by.y = "sample_name",all.x=TRUE,all.y = FALSE)
samples_num_fusion_tp53[is.na(N)]$N=0

samples_num_fusion_tp53$is.TP53_mutated=as.factor(samples_num_fusion_tp53$is.TP53_mutated)
samples_num_fusion_tp53[is.TP53_mutated==0,is.TP53_mutated:=as.factor("No")]
samples_num_fusion_tp53[is.TP53_mutated==1,is.TP53_mutated:=as.factor("Yes")]
samples_num_fusion_tp53 = samples_num_fusion_tp53[!(cancer%in%c("UVM","THYM","THCA","PCPG","KIRC","KIRP","TGCT","KICH"))] #remove those cancer with very low TP53 mutation rates
TCGA_Cancers=sort(unique(samples_num_fusion_tp53$cancer))
a=c()
b=c()
c=c()
d=c()
for (counter in 1:25){
  samples_num_fusion_tp53_per_cancer = samples_num_fusion_tp53[cancer%like%TCGA_Cancers[counter]]  
  g <- ggplot(samples_num_fusion_tp53_per_cancer, aes(is.TP53_mutated, N))
  sts_no_TP53 <- boxplot.stats( samples_num_fusion_tp53_per_cancer[is.TP53_mutated=="No"]$N)$stats 
  sts_TP53 <- boxplot.stats( samples_num_fusion_tp53_per_cancer[is.TP53_mutated=="Yes"]$N)$stats 
  a[counter]=t.test(samples_num_fusion_tp53_per_cancer[is.TP53_mutated=="Yes"]$N,samples_num_fusion_tp53_per_cancer[is.TP53_mutated=="No"]$N,alternative = "less")$p.value
  b[counter]=wilcox.test(samples_num_fusion_tp53_per_cancer[is.TP53_mutated=="Yes"]$N,samples_num_fusion_tp53_per_cancer[is.TP53_mutated=="No"]$N,alternative = "less")$p.value
 # print(g + geom_boxplot(outlier.colour=NA,width = 0.5,lwd =0.7,aes(fill=factor(is.TP53_mutated)))+coord_cartesian(ylim = c(0 ,max(sts_no_TP53[5],sts_TP53[5])*1.05))+ labs(title=TCGA_Cancers[counter])+ stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="black", fill="black")+xlab("TP53 mutation") + ylab("Number of fusions")+theme(legend.position="none",axis.text.x = element_text(size=14),axis.text.y = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16)))
#  print(g +  geom_violin(width = 0.5,lwd =0.7,aes(fill=factor(is.TP53_mutated)))+coord_cartesian(ylim = c(0 ,max(samples_num_fusion_tp53_per_cancer[is.TP53_mutated=="Yes"]$N,samples_num_fusion_tp53_per_cancer[is.TP53_mutated=="No"]$N)*1.01))+ labs(title=TCGA_Cancers[counter])+ stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="black", fill="black")+xlab("TP53 mutation") + ylab("Number of fusions")+theme(legend.position="none",axis.text.x = element_text(size=14),axis.text.y = element_text(size=14),axis.title.x = element_text(size=16),axis.title.y = element_text(size=16)))
}

####plot the graph for average fusion by DEEPEST vs TP53 mutation frequency
plot(average_num_DEEPEST_fusions_vs_p53_freq$mutation_frequency, average_num_DEEPEST_fusions_vs_p53_freq$mean, xlim=c(0,95),ylim=c(0,9.5), xlab="TP53 mutation frequency (%)" ,ylab="Average number of fusions", pch=19,cex=1.4,cex.lab=1.4)
text(13,9, "Pearson's r = 0.49667", col = 1,cex=1.2)
text(13,8, "Spearman's rho = 0.6370", col = 1,cex=1.2)
abline(lm(average_num_DEEPEST_fusions_vs_p53_freq$mean ~ average_num_DEEPEST_fusions_vs_p53_freq$mutation_frequency),col=1,cex=20)

x_pos= average_num_DEEPEST_fusions_vs_p53_freq$mutation_frequency-2
x_pos[17]=x_pos[17]+8
x_pos[9]=x_pos[9]+8
x_pos[24]=x_pos[24]+8
x_pos[26]=x_pos[26]+8
x_pos[5]=x_pos[5]+6
x_pos[19]=x_pos[19]+8
x_pos[7]=x_pos[7]+8
x_pos[14]=x_pos[14]+4
x_pos[12]=x_pos[12]+7
x_pos[4]=x_pos[4]+3
x_pos[21]=x_pos[21]+1
y_pos=average_num_DEEPEST_fusions_vs_p53_freq$mean+0.5
y_pos[21]=y_pos[21]-0.3
y_pos[14]=y_pos[14]-0.3
y_pos[29]=y_pos[29]-0.3
y_pos[4]=y_pos[4]+0.3
text(x = x_pos,y = y_pos,labels= average_num_DEEPEST_fusions_vs_p53_freq$TCGA_Project, cex = 0.7,pos = 1)
