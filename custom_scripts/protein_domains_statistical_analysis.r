library(data.table)
library(binom)


#All input files are available in DEEPEST-Fusion Github repository, custom scripts folder at: https://github.com/salzmanlab/DEEPEST-Fusion/tree/master/custom_scripts/files 


######### Input files ###############
fusion_prot_domain_file = "DEEPEST-Fusion/custom_scripts/files/fusion_protein_domain.txt"  #the file containing protein domains for each detected fusion obtained by AGFuison
gene_num_domain_file = "DEEPEST-Fusion/custom_scripts/files/num_domains_per_gene.txt"  #the file containg the number of domains in each gene of the genome
fusions_domain_pair_file = "DEEPEST-Fusion/custom_scripts/files/fusion_domain_cooccurrence.txt" #the file containing the domain interactions observed in fusions
transcriptome_domain_pair_file = "GDEEPEST-Fusion/custom_scripts/files/transcriptome_domain_pair.txt"  #the file containg domain interaction for linear transcriptome
transcriptome_single_domain_file = "DEEPEST-Fusion/custom_scripts/files/transcriptome_single_domain.txt"  #the file containing single domain information across linear transcriptome
fusions_single_domain_file = "DEEPEST-Fusion/custom_scripts/files/fusion_single_domain.txt"  #the file containing single domain information across linear transcriptome
###########################################


############## read in input files ###############
fusion_prot_domain = fread(fusion_prot_domain_file,header = TRUE,sep = "\t")
gene_num_domain = fread(gene_num_domain_file,header = TRUE,sep = "\t")
fusions_domain_pair = fread(fusions_domain_pair_file,header = TRUE,sep = "\t")
transcriptome_domain_pair = fread(transcriptome_domain_pair_file,header = TRUE,sep = "\t")   
transcriptome_single_domain = fread(transcriptome_single_domain_file,header = TRUE,sep = "\t")  
fusions_single_domain = fread(fusions_single_domain_file,header = TRUE,sep = "\t")
################################################



########## domain analysis: 
fusion_prot_domain = fusion_prot_domain[!(fusion_effect=="Outside transcript boundry")]  #we want to focus on only those fusions for which we have domain information
fusion_prot_domain[,three_prime_protein_domain:=gsub(":None","",three_prime_protein_domain),by=1:nrow(fusion_prot_domain)]
fusion_prot_domain[,five_prime_protein_domain:=gsub(":None","",five_prime_protein_domain),by=1:nrow(fusion_prot_domain)]
fusion_prot_domain[,three_prime_protein_domain:=gsub("None:","",three_prime_protein_domain),by=1:nrow(fusion_prot_domain)]
fusion_prot_domain[,five_prime_protein_domain:=gsub("None:","",five_prime_protein_domain),by=1:nrow(fusion_prot_domain)]
fusion_prot_domain[,three_prime_protein_domain:=gsub("None","",three_prime_protein_domain),by=1:nrow(fusion_prot_domain)]
fusion_prot_domain[,five_prime_protein_domain:=gsub("None","",five_prime_protein_domain),by=1:nrow(fusion_prot_domain)]
fusion_prot_domain[,gene1:=strsplit(fusion,split="--")[[1]][1],by=1:nrow(fusion_prot_domain)]
fusion_prot_domain[,gene2:=strsplit(fusion,split="--")[[1]][2],by=1:nrow(fusion_prot_domain)]

## to count number of domains that are unique, merge files by renaming
fusion_prot_domain = merge(fusion_prot_domain,gene_num_domain[,list(`Gene stable ID`,number_of_domains)],by.x="gene_id_1",by.y="Gene stable ID",all.x=TRUE,all.y = FALSE)
setnames(fusion_prot_domain,old = "number_of_domains",new = "gene1_num_domains")
fusion_prot_domain = merge(fusion_prot_domain,gene_num_domain[,list(`Gene stable ID`,number_of_domains)],by.x="gene_id_2",by.y="Gene stable ID",all.x=TRUE,all.y = FALSE)
setnames(fusion_prot_domain,old = "number_of_domains",new = "gene2_num_domains")

## analysis for one domain parents
num_fusions_with_onedomain_parents = dim(unique(fusion_prot_domain[gene1_num_domains==1&gene2_num_domains==1,list(sample_name,gene1,gene2)]))[1] ## all fusions with 1 domain parents
num_fusions_with_twodomains_and_onedomain_parents = dim(unique(fusion_prot_domain[gene1_num_domains==1&gene2_num_domains==1&!(three_prime_protein_domain== "")&!(five_prime_protein_domain== ""),list(gene1,gene2,sample_name)]))[1] ## all fusions with two domains and one-domain parental genes
#-- should be 1/3*1/2*1/2 of the above under null
#################################### DONE

## protein domain composition
## for each domain pair, compute prob. of each domain marginally at 5' and 3' ends
fusion_prot_domain$five_prime_protein_domain[fusion_prot_domain$five_prime_protein_domain == ""] = "empty"
fusion_prot_domain$three_prime_protein_domain[fusion_prot_domain$three_prime_protein_domain == ""] = "empty"

fusion_prot_domain[,num_distinct_five_prime_domains:=length(unique(strsplit(five_prime_protein_domain,split = ":")[[1]])),by=1:nrow(fusion_prot_domain)]
fusion_prot_domain[,num_distinct_three_prime_domains:=length(unique(strsplit(three_prime_protein_domain,split = ":")[[1]])),by=1:nrow(fusion_prot_domain)]

#first find the frequency for five prime domains
fusion_junction_repeat = rep(fusion_prot_domain$junction,fusion_prot_domain$num_distinct_five_prime_domains)
domain_seq_vector = list()
for (counter in 1:nrow(fusion_prot_domain)){
  domain_seq_vector[[counter]] = unique(strsplit(fusion_prot_domain[counter,`five_prime_protein_domain`],split = ":")[[1]])  #need the number of distinct domains for each protein
}
domain_seq_vector_unlist = unlist(domain_seq_vector)
fusion_prot_domain_single_dt = data.table(fusion_junction_repeat,domain_seq_vector_unlist)
fusion_five_prime_domains_marginal_counts = fusion_prot_domain_single_dt[,.N,by=domain_seq_vector_unlist]
fusion_five_prime_domains_marginal_counts[,N1 := N/length(unique(fusion_prot_domain$junction)),by = 1:nrow(fusion_five_prime_domains_marginal_counts)]
setnames(fusion_five_prime_domains_marginal_counts,old=c("domain_seq_vector_unlist","N1"),new=c("domain","freq"))

#now find the marginal frequency for three prime domains
fusion_junction_repeat = rep(fusion_prot_domain$junction,fusion_prot_domain$num_distinct_three_prime_domains)
domain_seq_vector = list()
for (counter in 1:nrow(fusion_prot_domain)){
  domain_seq_vector[[counter]] = unique(strsplit(fusion_prot_domain[counter,`three_prime_protein_domain`],split = ":")[[1]])  #need the number of distinct domains for each protein
}
domain_seq_vector_unlist = unlist(domain_seq_vector)
fusion_prot_domain_single_dt = data.table(fusion_junction_repeat,domain_seq_vector_unlist)
fusion_three_prime_domains_marginal_counts = fusion_prot_domain_single_dt[,.N,by=domain_seq_vector_unlist]
fusion_three_prime_domains_marginal_counts[,N1 := N/length(unique(fusion_prot_domain$junction)),by = 1:nrow(fusion_three_prime_domains_marginal_counts)]
setnames(fusion_three_prime_domains_marginal_counts,old=c("domain_seq_vector_unlist","N1"),new=c("domain","freq"))


fusions_domain_pair = merge(fusions_domain_pair,fusion_five_prime_domains_marginal_counts,by.x="domain_1",by.y="domain",all.x=TRUE,all.y=FALSE)
setnames(fusions_domain_pair,old=c("freq","N"),new=c("marginal_freq_domain1_on_5prime_side","marginal_count_domain1_on_5prime_side"))
fusions_domain_pair = merge(fusions_domain_pair,fusion_three_prime_domains_marginal_counts,by.x="domain_2",by.y="domain",all.x=TRUE,all.y=FALSE)
setnames(fusions_domain_pair,old=c("freq","N"),new=c("marginal_freq_domain2_on_3prime_side","marginal_count_domain2_on_3prime_side"))
fusions_domain_pair = merge(fusions_domain_pair,fusion_five_prime_domains_marginal_counts,by.x="domain_2",by.y="domain",all.x=TRUE,all.y=FALSE)
setnames(fusions_domain_pair,old=c("freq","N"),new=c("marginal_freq_domain2_on_5prime_side","marginal_count_domain2_on_5prime_side"))
fusions_domain_pair = merge(fusions_domain_pair,fusion_three_prime_domains_marginal_counts,by.x="domain_1",by.y="domain",all.x=TRUE,all.y=FALSE)
setnames(fusions_domain_pair,old=c("freq","N"),new=c("marginal_freq_domain1_on_3prime_side","marginal_count_domain1_on_3prime_side"))
fusions_domain_pair[,same_domains:=0]
fusions_domain_pair[domain_1==domain_2,same_domains:=1]


## set up for normal domain distribution: MERGE files with null rates
fusions_domain_pair = merge(fusions_domain_pair,transcriptome_single_domain[,list(Pfam_name,freq1)],by.x = "domain_1",by.y = "Pfam_name", all.x = TRUE,all.y = FALSE)
setnames(fusions_domain_pair,"freq1","normal_freq_domain_1")
fusions_domain_pair = merge(fusions_domain_pair,transcriptome_single_domain[,list(Pfam_name,freq1)],by.x = "domain_2",by.y= "Pfam_name",all.x = TRUE,all.y = FALSE)
setnames(fusions_domain_pair,"freq1","normal_freq_domain_2")

## Not sure why NAs, but remove them:
fusions_domain_pair[is.na(normal_freq_domain_1),normal_freq_domain_1:=0]
fusions_domain_pair[is.na(normal_freq_domain_2),normal_freq_domain_2:=0]
fusions_domain_pair[is.na(marginal_count_domain1_on_5prime_side),marginal_count_domain1_on_5prime_side:=0]
fusions_domain_pair[is.na(marginal_count_domain1_on_3prime_side),marginal_count_domain1_on_3prime_side:=0]
fusions_domain_pair[is.na(marginal_count_domain2_on_5prime_side),marginal_count_domain2_on_5prime_side:=0]
fusions_domain_pair[is.na(marginal_count_domain2_on_3prime_side),marginal_count_domain2_on_3prime_side:=0]
fusions_domain_pair[is.na(marginal_freq_domain1_on_5prime_side),marginal_freq_domain1_on_5prime_side:=0]
fusions_domain_pair[is.na(marginal_freq_domain1_on_3prime_side),marginal_freq_domain1_on_3prime_side:=0]
fusions_domain_pair[is.na(marginal_freq_domain2_on_5prime_side),marginal_freq_domain2_on_5prime_side:=0]
fusions_domain_pair[is.na(marginal_freq_domain2_on_3prime_side),marginal_freq_domain2_on_3prime_side:=0]

# set up for statistical test of whether the counts
## Roozbeh, please check this:
ep = 1/10^5
fusions_single_domain = merge(fusions_single_domain,transcriptome_single_domain[,list(Pfam_name,Pfam_description,freq1)],by.x = "domain",by.y = "Pfam_name",all.x = TRUE,all.y = FALSE)
setnames(fusions_single_domain,"freq1","normal_freq")
fusions_single_domain[,binom.lower:=binom.confint(absolute_number,dim(fusion_prot_domain)[1], method = "exact", conf.level = 1-ep)$lower,by=1:nrow(fusions_single_domain)]
fusions_single_domain[,fold_change:=freq/normal_freq,by = 1:nrow(fusions_single_domain)]
fusions_single_domain[,p_value:=pbinom(absolute_number, round(absolute_number/freq), normal_freq, lower.tail = FALSE),by=1:nrow(fusions_single_domain)]

fusions_domain_pair[,pairbinom.lower:=binom.confint(count_new,dim(fusion_prot_domain)[1], method = "exact",conf.level = 1-10^-6)$lower,by=1:nrow(fusions_domain_pair)] ## pair analysis
fusions_domain_pair[same_domains==1,null_prob:=normal_freq_domain_1*normal_freq_domain_2]
fusions_domain_pair[same_domains==0,null_prob:=normal_freq_domain_1*normal_freq_domain_2+normal_freq_domain_1*normal_freq_domain_2]
fusions_domain_pair[,foldchange := frequency/null_prob]
fusions_domain_pair[,p_value:=pbinom(count_new, round(count_new/frequency), null_prob, lower.tail = FALSE),by=1:nrow(fusions_domain_pair)]

## use a hockberg correction--- right now bonferonni
## this gives which domains are enirched in fusions and their point estimate of fold enrichment in second column
enriched_single_domains = fusions_single_domain[binom.lower > normal_freq] # domains that are enriched compared to 1-10-^5
enriched_domain_pairs = fusions_domain_pair[(count_new>1)& !(domain_domain_interaction %like%"empty") & (pairbinom.lower > null_prob)] ## enriched domain pairs in fusions

write.table(enriched_single_domains,"enriched_fusion_single_domain.txt",row.names = FALSE,quote = FALSE,sep="\t")
write.table(enriched_domain_pairs,"enriched_fusion_domain_pair.txt",row.names = FALSE,quote = FALSE,sep="\t")

## many domain-domain pair enriched domains are chromatin remodeling
fusions_domain_pair[(count_new>4) &!(domain_domain_interaction %like%"empty")&(pairbinom.lower> normal_freq_domain_2*normal_freq_domain_1)]