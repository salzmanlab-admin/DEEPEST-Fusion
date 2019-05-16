require(data.table)
require(stringr)
require(compare)
library(plyr)
library(dplyr)
library(dcGOR)

#All input files are available in DEEPEST-Fusion Github repository, custom scripts folder at: https://github.com/salzmanlab/DEEPEST-Fusion/tree/master/custom_scripts/files 

####### Inputs files ############
enriched_single_domains_file = "DEEPEST-Fusion/custom_scripts/files/enriched_fusion_single_domain.txt" 
enriched_domain_pairs_file =   "DEEPEST-Fusion/custom_scripts/files/enriched_fusion_domain_pair.txt"
Pfam_A_hmm_file = "DEEPEST-Fusion/custom_scripts/files/Pfam-A.hmm.dat"
######################################

###### read in input files ####################
enriched_single_domains = fread(enriched_single_domains_file,header = TRUE,sep="\t")
enriched_domain_pairs_file = fread(enriched_domain_pairs_file,header = TRUE,sep="\t")
Pfam_A_hmm = fread(Pfam_A_hmm_file,sep="\t",header = FALSE)    # we extract the pfam_ids and pfa_names based on this files
################################################



########################################################
############## building a data frame that pairs each pfam_id with its corresponding pfam_name
PF_ids = Pfam_A_hmm[V1%like%"#=GF AC"]    
PF_names = Pfam_A_hmm[V1%like%"#=GF ID"]
PF_description = Pfam_A_hmm[V1%like%"#=GF DE"]
PF_ids[,V1:=strsplit(V1,split = "#=GF AC")[[1]][2],by=1:nrow(PF_ids)]
PF_names[,V1:=strsplit(V1,split = "#=GF ID")[[1]][2],by=1:nrow(PF_names)]
PF_description[,V1:=strsplit(V1,split = "#=GF DE")[[1]][2],by=1:nrow(PF_description)]
PF_ids[,V1:=gsub(" ","",V1),by = 1:nrow(PF_ids)]
PF_names[,V1:=gsub(" ","",V1),by = 1:nrow(PF_names)]
PF_description[,V1:=gsub(" ","",V1),by = 1:nrow(PF_description)]
PF_ids[,V1:=strsplit(V1,split = ".",fixed = TRUE)[[1]][1],by=1:nrow(PF_ids)]
PF_id_PF_name = data.table(PF_ids,PF_names,PF_description)  #the data frame that shows the pfam name for each pfam_id
names(PF_id_PF_name) = c("Pfam_id","Pfam_name","Pfam_description")
#######################################################


##GO enrichment analysis based on dcGOR package
Pfam = dcRDataLoader('Pfam')
eoutput_BY_binom = dcEnrichment(enriched_single_domains$Pfam_id, domain="Pfam", ontology="GOBP",p.adjust.method="BY",test = "BinomialTest")
eoutput_BH_binom = dcEnrichment(enriched_single_domains$Pfam_id, domain="Pfam", ontology="GOBP",p.adjust.method="BH",test = "BinomialTest")
eoutput_BH_hyper = dcEnrichment(enriched_single_domains$Pfam_id, domain="Pfam", ontology="GOBP",p.adjust.method="BH",test = "HypergeoTest")
eoutput_BY_hyper = dcEnrichment(enriched_single_domains$Pfam_id, domain="Pfam", ontology="GOBP",p.adjust.method="BY",test = "HypergeoTest")

Biological_processes_GO_analysis_fusion_domains=data.table(view(eoutput_BY_binom, top_num=1000, sortBy="pvalue", details=TRUE))
Biological_processes_GO_analysis_fusion_domains[,fold_enrichement:=nOverlap/nGroup/(nAnno/3241),by=1:nrow(Biological_processes_GO_analysis_fusion_domains)]
write.table(Biological_processes_GO_analysis_fusion_domains,"GO_analysis_for_enriched_fusion_domains.txt",quote = FALSE,row.names = FALSE,sep="\t")
