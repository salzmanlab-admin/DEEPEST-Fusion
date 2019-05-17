require(data.table)
require(stringr)
require(compare)
library(dplyr)
library(DescTools)

#All input files are available in DEEPEST-Fusion Github repository, custom scripts folder at: https://github.com/salzmanlab/DEEPEST-Fusion/tree/master/custom_scripts/files 


######################## Inputs files ###################
five_prime_recurrent_file = "DEEPEST-Fusion/custom_scripts/files/five_prime_recurrent_genes.txt"
three_prime_recurrent_file = "DEEPEST-Fusion/custom_scripts/files/three_prime_recurrent_genes.txt"
DEEPEST_fusions_with_ID_file = "DEEPEST-Fusion/custom_scripts/files/DEEPEST_fusions_with_ID.txt"
COSMIC_genes_file = "DEEPEST-Fusion/custom_scripts/files/cancer_gene_census.csv"
DEEPEST_SBT_GTEx_hits_file = "DEEPEST-Fusion/custom_scripts/files/DEEPEST_SBT_GTEx_hits.txt"
#####################################################

#### read in input files ########
five_prime_genes = fread(five_prime_recurrent_file,header = TRUE,sep = "\t")   #genes that have partnerd with many genes as five_prime genes in the fusions
three_prime_genes = fread(three_prime_recurrent_file,header = TRUE,sep = "\t")
DEEPEST_fusions_with_ID = fread(DEEPEST_fusions_with_ID_file,header = TRUE,sep = "\t")
cosmic_genes = fread(COSMIC_genes_file,sep = ",",header = TRUE)
DEEPEST_SBT_GTEx_hits = fread(DEEPEST_SBT_GTEx_hits_file,sep="\t",header = TRUE)
###############################


DEEPEST_SBT_GTEx_hits[,gene1 := strsplit(fusion,split="--")[[1]][1],by=1:nrow(DEEPEST_SBT_GTEx_hits)]
DEEPEST_SBT_GTEx_hits[,gene2 := strsplit(fusion,split="--")[[1]][2],by=1:nrow(DEEPEST_SBT_GTEx_hits)]
DEEPEST_fusions_with_ID[,pos1:=as.numeric(strsplit(junction,split = ":",fixed = TRUE)[[1]][3]),by=1:nrow(DEEPEST_fusions_with_ID)]
DEEPEST_fusions_with_ID[,pos2:=as.numeric(strsplit(junction,split = ":",fixed = TRUE)[[1]][7]),by=1:nrow(DEEPEST_fusions_with_ID)]
DEEPEST_fusions_with_ID[,junc:=strsplit(junction,split = ":",fixed = TRUE)[[1]][9],by=1:nrow(DEEPEST_fusions_with_ID)]
DEEPEST_fusions_with_ID[,type:="machete"]
DEEPEST_fusions_with_ID[junc=="rev"&abs(pos1-pos2) < 1000000,type:="knife"]
TCGA_fusions_unique = unique(DEEPEST_fusions_with_ID[,list(sample_name,fusion_unified)])
num_each_fusion = TCGA_fusions_unique[,.N,by=fusion_unified]
setnames(num_each_fusion,"N","occurrence")

#recurrent partner genes
five_prime_genes_recurrent = five_prime_genes[num_three_prime_partners>1]    #only work with those genes that have appeared more than once as five prime gene
three_prime_genes_recurrent = three_prime_genes[num_five_prime_partners>1]
five_prime_genes_recurrent = setorder(five_prime_genes_recurrent,-num_three_prime_partners) 
three_prime_genes_recurrent = setorder(three_prime_genes_recurrent,-num_five_prime_partners)
dist_three_prime_partners = five_prime_genes_recurrent[,.N,by=num_three_prime_partners]  #
dist_five_prime_partners = three_prime_genes_recurrent[,.N,by=num_five_prime_partners]

fusions_recurrent = num_each_fusion[occurrence>1]
fusions_recurrent = setorder(fusions_recurrent,-occurrence)
dist_fusions = fusions_recurrent[,.N,by=occurrence]


num_fusions = nrow(unique(DEEPEST_fusions_with_ID[,list(sample_name,fusion_unified)]))   #total number of detected fusions 
num_all_five_prime_genes = length(unique(DEEPEST_fusions_with_ID$gene1))   #total number of 5' genes in the detected fusions
num_all_three_prime_genes = length(unique(DEEPEST_fusions_with_ID$gene2))
num_distinct_fusions =length(unique(TCGA_fusions_unique$fusion_unified))

#################################################################################
#################################################################################
#################### Compute p-value for the number of recurrent genes #####

#first for recurrent 5' genes
n = 20000
k_n = num_fusions
c = dist_three_prime_partners$num_three_prime_partners[1]
t = k_n/(n^(1-1/c))

prob_five_prime_genes = log10(ppois(dist_three_prime_partners$N[1]-1, t^c/factorial(c), log = FALSE,lower.tail = FALSE)) #we compute the logarithm of the p-value for all recurrent partner genes
for (counter in 2:nrow(dist_three_prime_partners)){
  
  n = 20000 - dist_three_prime_partners$N[counter-1]
  k_n = num_fusions - dist_three_prime_partners$N[counter-1]*dist_three_prime_partners$num_three_prime_partners[counter-1]
  c = dist_three_prime_partners$num_three_prime_partners[counter]
  t = k_n/(n^(1-1/c))
  prob_five_prime_genes =  prob_five_prime_genes+ log10(ppois(dist_three_prime_partners$N[counter]-1, t^c/factorial(c), log = FALSE,lower.tail = FALSE))
}

#Now  recurrent 3' genes
c = dist_five_prime_partners$num_five_prime_partners[1]
t = k_n/(n^(1-1/c))
prob_three_prime_genes = log10(ppois(dist_five_prime_partners$N[1]-1, t^c/factorial(c), log = FALSE,lower.tail = FALSE))
for (counter in 2:nrow(dist_five_prime_partners)){
  n = 20000 - dist_five_prime_partners$N[counter-1]
  k_n = num_fusions - dist_five_prime_partners$N[counter-1]*dist_five_prime_partners$num_five_prime_partners[counter-1]
  c = dist_five_prime_partners$num_five_prime_partners[counter]
  t = k_n/(n^(1-1/c))
  prob_three_prime_genes =  prob_three_prime_genes + log10(ppois(dist_five_prime_partners$N[counter]-1, t^c/factorial(c), log = FALSE,lower.tail = FALSE))
}



#########################################################################################################
#########################################################################################################
########Compute CI and expected number of recurrent genes for each number of c (minimum number of partners) based on BHY FDR control


######first five prime recurrent genes
n = 20000
k_n = num_fusions   #total number of detected fusions 
m = nrow(dist_three_prime_partners)   #the number of hypotheses needed for BHY method
c_m = sum(1/1:m)
alpha = 0.01/c_m
k = 0
all_dist_three_prime_partners = data.table(num_three_prime_partners=65:2)
all_dist_three_prime_partners = merge(all_dist_three_prime_partners,dist_three_prime_partners,by.x="num_three_prime_partners",by.y="num_three_prime_partners",all.x=TRUE,all.y=FALSE)
all_dist_three_prime_partners[is.na(N),N:=0]
all_dist_three_prime_partners = setorder(all_dist_three_prime_partners,-num_three_prime_partners)
dist_three_prime_partners = all_dist_three_prime_partners
cum_sum_five_prime_recurrent_genes = cumsum(dist_three_prime_partners$N)
p_value_five_prime_genes = c()
for (counter in 1:nrow(dist_three_prime_partners)){
  k = k + 1   #the parameter K needed for BHY FDR control
  c = dist_three_prime_partners$num_three_prime_partners[counter] 
  t = k_n/(n^(1-1/c))
  p_value_five_prime_genes[counter] = ppois(cum_sum_five_prime_recurrent_genes[counter],t^c/factorial(c),lower.tail = FALSE)
}


##Now compute expectation and its upper CI based on null
m = 65   #the number of hypotheses needed for BHY method
c_m = sum(1/1:m)
k = 0
upper_CI_num_partners = c()  #based on qpois function
expected_freq = c()
for (c in 65:2){
  k = k + 1
  k_n = num_fusions
  t = k_n/(n^(1-1/c))
  upper_CI_num_partners[c] = qpois(0.01*k/m/c_m,t^c/factorial(c), lower.tail = FALSE)
  expected_freq [c] = t^c/factorial(c)
}

expected_freq[expected_freq < 0.013] = 0.013   #i add these values as it is difficulat to plot absolute zero values in a logarithmic plot
upper_CI_num_partners[upper_CI_num_partners < 0.013] = 0.013

#now 3' recurrent genes
k = 0
m = nrow(dist_five_prime_partners)   #the number of hypotheses needed for BHY method
c_m = sum(1/1:m)
all_dist_five_prime_partners = data.table(num_five_prime_partners=61:2)
all_dist_five_prime_partners = merge(all_dist_five_prime_partners,dist_five_prime_partners,by.x="num_five_prime_partners",by.y="num_five_prime_partners",all.x=TRUE,all.y=FALSE)
all_dist_five_prime_partners[is.na(N),N:=0]
all_dist_five_prime_partners = setorder(all_dist_five_prime_partners,-num_five_prime_partners)
dist_five_prime_partners = all_dist_five_prime_partners

cum_sum_three_prime_recurrent_genes = cumsum(dist_five_prime_partners$N)
p_value_three_prime_genes = c()
for (counter in 1:nrow(dist_five_prime_partners)){
  k = k + 1
  k_n = num_fusions
  c = dist_five_prime_partners$num_five_prime_partners[counter]
  t = k_n/(n^(1-1/c))
  p_value_three_prime_genes[counter] = ppois(cum_sum_three_prime_recurrent_genes[counter],t^c/factorial(c),lower.tail = FALSE)
}


###########################################################################
######################### compute the statsitics for the recurrent fusions 
n = 20000 *(20000-1)
k_n = nrow(TCGA_fusions_unique)
c = dist_fusions$occurrence[1]
t = k_n/(n^(1-1/c))
prob_recurrent_fusions = log10(ppois(dist_fusions$N[1]-1, t^c/factorial(c), log = FALSE,lower.tail = FALSE)) #we compute the logarithm of the p-value for all recurrent partner genes

m = nrow(dist_fusions)   #the number of hypotheses needed for BHY method
c_m = sum(1/1:m)
k = 1

upper_CI_recurrent_fusions = c()  #upper CI based on qpois command 
expected_freq_recurrent_fusions = c()

k = 0
for (counter in 1:nrow(dist_fusions)){
  n = 20000*(20000-1)
  k_n = num_fusions
  k = k + 1
  c = dist_fusions$occurrence[counter]
  t = k_n/(n^(1-1/c))
  upper_CI_recurrent_fusions[counter] = qpois(alpha*k/m/c_m,t^c/factorial(c), lower.tail = FALSE)
  expected_freq_recurrent_fusions [counter] = t^c/factorial(c)
}


for (counter in 2:nrow(dist_fusions)){
  n = 20000*(20000-1) - dist_fusions$N[counter-1]
  k_n = num_fusions
  c = dist_fusions$occurrence[counter]
  t = k_n/(n^(1-1/c))
  prob_recurrent_fusions =  prob_recurrent_fusions + log10(ppois(dist_fusions$N[counter]-1, t^c/factorial(c), log = FALSE,lower.tail = FALSE))
}
cum_sum_recurrent_fusions = cumsum(dist_fusions$N)

############################################################
############################################################
#here we want to find the percentage of significant samples (samples that have at least one fusion with significant recurrent gene) in GTEx and TCGA as a function of k (number of partners)
num_partners = 7:51

#first tcga sample
num_sig_genes = c()
num_sig_tcga_samples = c()
significant_five_prime_recurrent_genes = list()
significant_three_prime_recurrent_genes = list()
tcga_fusions_with_significant_recurrent_genes = list()
gtex_fusions_with_significant_recurrent_genes = list()
significant_recurrent_genes = list()
tcga_fusions_extreme_fraction = c()
gtex_fusions_extreme_fraction = c()
tcga_samples_extreme_fraction = c()
gtex_samples_extreme_fraction = c()
tcga_samples_nonextreme_fraction = c()
gtex_samples_nonextreme_fraction = c()
for (k in 1:45){
  significant_five_prime_recurrent_genes[[k]] = five_prime_genes[num_three_prime_partners >= num_partners[k]]$all_gene1
  significant_three_prime_recurrent_genes[[k]] = three_prime_genes[num_five_prime_partners >= num_partners[k]]$all_gene2
  tcga_fusions_with_significant_recurrent_genes[[k]] = unique(DEEPEST_fusions_with_ID[(gene1_unified%in%significant_five_prime_recurrent_genes[[k]]) |(gene2_unified%in%significant_three_prime_recurrent_genes[[k]])])
  tcga_fusions_extreme_fraction[k]=nrow(unique(tcga_fusions_with_significant_recurrent_genes[[k]][type=="machete",list(sample_name,fusion)])) / nrow(unique(tcga_fusions_with_significant_recurrent_genes[[k]][,list(sample_name,fusion)]))
  tcga_samples_extreme_fraction[k]=nrow(unique(tcga_fusions_with_significant_recurrent_genes[[k]][type=="machete",list(sample_name)])) / 9946*100
  tcga_samples_nonextreme_fraction[k]=nrow(unique(tcga_fusions_with_significant_recurrent_genes[[k]][type=="knife",list(sample_name)])) / 9946*100
  significant_recurrent_genes[[k]] = unique(c(significant_five_prime_recurrent_genes[[k]],significant_three_prime_recurrent_genes[[k]]))
  num_sig_tcga_samples[k] = length(unique(tcga_fusions_with_significant_recurrent_genes$sample_name))/9946*100  # normalize to the total number of tcga samples we have run
  num_sig_genes[k] = length(significant_recurrent_genes)
}

#now gtex samples
num_sig_gtex_samples = c()   # a vector with the fraction of significant samplesas a function of k (number of partners)
for (k in 1:45){
  significant_five_prime_recurrent_genes[[k]] = five_prime_genes[num_three_prime_partners >= num_partners[k]]$all_gene1
  significant_three_prime_recurrent_genes[[k]] = three_prime_genes[num_five_prime_partners >= num_partners[k]]$all_gene2
  gtex_fusions_with_significant_recurrent_genes[[k]] = unique(DEEPEST_SBT_GTEx_hits[(gene1%in%significant_five_prime_recurrent_genes[[k]]) |(gene2%in%significant_three_prime_recurrent_genes[[k]])])
  gtex_fusions_extreme_fraction[k]=nrow(unique(gtex_fusions_with_significant_recurrent_genes[[k]][type=="machete",list(sample_name,fusion)])) / nrow(unique(gtex_fusions_with_significant_recurrent_genes[[k]][,list(sample_name,fusion)]))
  gtex_samples_extreme_fraction[k]=nrow(unique(gtex_fusions_with_significant_recurrent_genes[[k]][type=="machete",list(sample_name)])) / 287*100
  gtex_samples_nonextreme_fraction[k]=nrow(unique(gtex_fusions_with_significant_recurrent_genes[[k]][type=="knife",list(sample_name)])) / 287*100
  significant_recurrent_genes[[k]] = unique(c(significant_five_prime_recurrent_genes[[k]],significant_three_prime_recurrent_genes[[k]]))
  num_sig_gtex_samples[k] = length(unique(gtex_fusions_with_significant_recurrent_genes$sample_name)) / 287 *100  # we normalize to the total number of gtex samples we have run
}

marks = c(7,10,20,30,40,50)
plot(num_partners,tcga_fusions_extreme_fraction*100,xat = "n",xlab = "Significantly fused gene",ylab = "% of extreme fusions",pch = 19,ylim = c(0,100),xlim=c(8,49), cex = 1.1,cex.lab=1.3,col="darkcyan")
points(num_partners,gtex_fusions_extreme_fraction*100,pch = 19,cex.lab = 1.3,col = "darkorange")
legend(x = 34, y = 66, legend = c("TCGA tumor","GTEx"), pch = c(19,19), col = c("darkcyan", "darkorange"),cex = 1)
axis(1,at = marks,labels = marks)

plot(num_partners,tcga_samples_extreme_fraction,xlab = "Significantly fused gene",ylab= "% of samples",pch=19,ylim = c(0,50),cex = 1.1,cex.lab=1.3,col="darkcyan")
points(num_partners,gtex_samples_extreme_fraction,pch=19,cex.lab=1.3,col="darkorange")
legend(x = 34, y = 66, legend = c("TCGA tumor","GTEx"), pch = c(19,19), col = c("darkcyan", "darkorange"),cex=1)

plot(num_partners,num_sig_tcga_samples,xlab = "Significantly fused gene",ylab= "Fraction of samples",pch=19,ylim = c(0.01,50),cex = 1.1,cex.lab=1.3,col="darkcyan")
points(num_partners,num_sig_gtex_samples,pch=19,cex.lab=1.3,col="darkorange")
legend(x = 34, y = 46, legend = c("TCGA tumor","GTEx"), pch = c(19,19), col = c("darkcyan", "darkorange"),cex=1)


###################################################################################
####################################################################################
##################Plotting the results############################################

#plot the graph for recurrent genes both 5' and 3'recurrent genes
marks = c(0.01,0.1,1,10,100,1000,10000)
marks_labels = c(0,0.1,1,10,100,1000,10000)
marks_x=c(2,10,20,30,40,50,60,65)
plot(dist_three_prime_partners$num_three_prime_partners,cum_sum_five_prime_recurrent_genes,xlab = "Number of partners",ylab= "Number of genes", yaxt="n",pch=19,log = "y",ylim = c(0.01,30000),xlim=c(4,65),cex = 1.15,cex.lab=1.4,xaxt = "n")
points(dist_five_prime_partners$num_five_prime_partners,cum_sum_three_prime_recurrent_genes,pch = 18,cex = 1.1,cex.lab = 1.4,col = "orange")
points(65:2,expected_freq[65:2],col = "blue",pch=15,cex = 1.2)
points(65:2,upper_CI_num_partners[65:2],col = "red",pch=17,cex = 1.3)
legend(x = 43, y = 15000, legend = c("5' Recurrent genes","3' Recurrent genes","Upper 99% CI","Expected"), pch = c(19,18,17,15), col = c("black", "orange","red","blue"),cex=1)
axis(2,las=1,at=marks,labels=marks_labels)
axis(1,at=marks_x,labels=marks_x)


#plot the graph for the recurrent fusions graph
marks = c(0.01,0.1,1,10,100,1000,10000)
marks_labels = c(0,0.1,1,10,100,1000,10000)
marks_x=c(0,50,100,150,182)
expected_freq_recurrent_fusions[expected_freq_recurrent_fusions<0.013] = 0.013
upper_CI_recurrent_fusions[upper_CI_recurrent_fusions < 0.013] = 0.013
plot(dist_fusions$occurrence,cum_sum_recurrent_fusions,xlab = "Number of occurrences",ylab= "Number of fusions", yaxt="n",pch=19,log = "y",ylim = c(0.01,10000),cex = 1.1,cex.lab=1.4)
points(dist_fusions$occurrence,expected_freq_recurrent_fusions,col = "blue",pch=15,cex = 1.2)
points(dist_fusions$occurrence,upper_CI_recurrent_fusions,col = "red",pch=18,cex = 0.9)
legend(x = 120, y = 3500, legend = c("Recurrent fusions","Upper 99% CI","Expected"), pch = c(19,18,15), col = c("black", "red","blue"),cex=1)
axis(2,las=1,at=marks,labels=marks)
axis(1,at=marks_x,labels=marks_x)
