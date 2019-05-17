# Custom scripts
The provided custom scripts in this folder have been used to perform downstream systems-level analysis on the fusions detected by DEEPEST on the TCGA corpus. The scripts also generate the figures that have been used to depict the analysis. All scripts are in R and all needed files for running them are provided in "files" folder. The following packages needed to be installed in R before running the scripts:  

- data.table
- stringr
- compare
- dplyr
- ggplot2
- binom
- plyr
- dcGOR
- DescTools

# Description of scripts:

- tp53_analysis.R: this script identifies the correlation between the abudnance of fusions and TP53 mutation rate as two orthogonal measures of genomic instability. In addition to DEEPEST, the script calculates the correlation coefficients for TumorFusions and (Gao et al., 2018) studies as well (correspondin to the results in Figure 2A). 
- kinase_cosmic_analysis.R: this script performs a statistical analysis to identify TCGA tumor types with high enrichment for fusions involving kinase or COSMIC genes (corresponding to the results presented in Figure 2B and 5A in the paper)

- recurrent_partners_statistical_analysis.R: this script performs the statistical analysis based on generalized birthday problem to find significantly-fused genes and also show the recurrent fusions are highly enriched in fuions. Then it compares the distinct profiles of significantly fused genes in TCGA and GTEx samples (corresponding to results in Figure 4).   

- protein_domains_statistical_analysis.r: This script finds enriched single protein domains and protein domains across fusion proteins when their frequencies are compared to their null frequencies in the reference transcriptome (corresponding to Dataset S3).  

- protein_domains_GO_enrichment_analysis.R: This script performs GO enrichment analysis to find enriched biological functions in the  enriched protein domains(corresponding to Figure 5B).  
