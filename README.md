# DEEPEST-Fusion

DEEPEST-Fusion (Data-Enriched Efficient PrEcise STatistical Fusion detection) is a statistical fusion detection algorithm developed in the Salzman Lab (http://salzmanlab.stanford.edu/).  DEEPEST is a completele suit of splice detection algorithms that can detect linear junctions, backsplice junctions (for circRNA detection), and fusion junctions. This repository contains the instructions on how to use the online tool or run the tool on a local cluster.

# Online tool with web interface:
DEEPEST-Fusion online tool with a web interface is now publicly available on [Cancer Genomics Cloud (CGC)](http://www.cancergenomicscloud.org/) at: https://cgc.sbgenomics.com/public/apps#jordanski.milos/deepest-fusion/deepest-fusion/

The tool is easy-to-use with only few clicks and can be run either on the dataset uploaded by the user or on public RNA-Seq datasets already available on CGC (such as TCGA, TARGET, CCLE, ...). 

To run the online tool on the cloud, you need to login to CGG [here](https://cgc-accounts.sbgenomics.com/auth/login?next=https%3A%2F%2Fcgc-accounts.sbgenomics.com%2Foauth2%2Fauthorization%3Fresponse_type%3Dcode%26client_id%3D08bbb98f354e4554bd7fd315de64d955%26redirect_uri%3Dhttps%253A%252F%252Fcgc.sbgenomics.com%252Foauth2%252Fredirect%26scope%3Dopenid%26state%3DDlQ4PIZFvqpWnrod5lOzyVG6M9qcLf%26nonce%3D2AKOsefdeicsyDctFCyug2LBl6KyL8). Each new account would automatically have $300 pilot funds. Our current estimate for running DEEPEST on typical RNA-Seq datasets (such as TCGA or GTEx) is $4-$5 per sample.    


# Software requirements for running on a local cluster

- Anaconda
- SciPy 1.1.0
- Bowtie2 2.2.9
- Python 3.4 with Biopython/1.70 installed
- java
- R 3.5.1
- Trim Galore! 0.4.4

# DEEPEST-Fusion main script

All computational steps including any alignment step required for running DEEPEST-Fusion has been packaged in a single JSON file "DEEPEST-Fusion_pipeline.json". This script can be run on any local cluster using an input JSON file that provides the paths for reference index files and input RNA-seq files. To run the script, Rabix toolkit should be pre-installed on the local cluster.  

# Input file

All input parameters required for running DEEPEST-Fusion JSON script should be provided via an input JSON file "DEEPEST-Fusion_input.json". The following parameters should be set in the input JSON file:

- Bowtie2 index files for the reference genome
- Bowtie2 index files for the regular junctions
- Bowtie2 index files for the scrambled junctions
- Bowtie2 index files for ribosome
- Bowtie2 index files for the transcriptome
- Bowtie2 index files for indel junctions (for junctions with up to 5 symmetric indels at the splice site)
- Fastq files for the input RNA-Seq data
- Pickle files for genes/exons annotation
- Pickle file for the known fusions list (a list of known fusions constructed based on ChimerDB 3.0 curated list of known cancer fusions)
- Bowtie2 index files for known fusions

Note: Since index files are too large to be uploaded to githiub, we provide original reference fasta files (genome, transcriptome, regular junctions, scrambled junctions, ribosome) along with needed scripts/instructions for generating index files here: 
https://github.com/salzmanlab/DEEPEST-Fusion/tree/master/reference_files 

# Toolkit for executing DEEPEST-Fusion JSON script

For running DEEPEST-Fusion JSON script on a local cluster, Rabix should be installed first. Rabix is an open-source tool developed by Seven Bridges, that can be used to run a computational workflow written in Common Workflow Language (CWL) on a locul cluster. More information on how to install Rabix can be found in this GitHub repositiory: https://github.com/rabix/bunny  

# Example batch script for running DEEPEST-Fusion on a local cluster

An example batch script "DEEPEST-Fusion_submit_job.sbatch", based on job scheduler Slurm has been provided. In the batch script file the path to where Rabix has been installed, DEEPEST-Fusion pipeline JSON file (DEEPEST-Fusion_pipeline.json), and DEEPEST-Fusion input JSON file (DEEPEST-Fusion_input.json) should be provided. 

# Output files

Three primary report files containing reported fusion junction with their corresponding statsitical scores and number of various types of aligned reads can be found as follows:

- modified-MACHETE report file (based on FarJunctions database): Knife_and_MACHETE_Known_fusions_parallel_rev_\*/root/MACHETE_AppendNaiveReptParallel/\[Dataset name\]\_naive_report_AppendedMACHETE\_Parallel.txt
- modified-MACHETE report file (based on known fusion database): Knife_and_MACHETE_Known_fusions_parallel_rev_\*/root/MACHETE_AppendNaiveReptParallel_Known/\[Dataset name\]\_naive_report_Appended_MACHETE_Parallel\_Known.txt
- KNIFE report file: Knife_and_MACHETE_Known_fusions_parallel_rev_\*/root/KNIFE_GLM_model/\[Dataset name\]\_1\_circJuncProbs.txt_cdf  

In each report file, high quality junctions are called by imposing thresholds on the statistical scores (as described in the paper). The list of reported fusion junctions which is the output of the junction nomination component (first computational component in DEEPEST-Fusion) is the union of fusion junctions called from these 3 report files. This list of reported fusions should then undergo the statistical refinement step (the second componenet of DEEPEST-Fusion) which is based on Sequence Bloom Trees.

# CWL scripts for implementing Sequence Bloom Trees

All scripts needed for implementing Sequence Bloom Tree (SBT) filters for an RNA-Seq database can be found in this github repository: https://github.com/elehnertSBG/SBT-Apps. Mor information on the SBT algorithm and the order in which the SBT CWL scripts should be run can be found in the original SBT manual: https://www.cs.cmu.edu/~ckingsf/software/bloomtree/sbt-manual.pdf 

# Citation

Dehghannasiri, R., Freeman, D.E., Jordanski, M., Hsieh, G.L., Damljanovic, A., Lehnert, E. and Salzman, J., 2019. `Improved detection of gene fusions by applying statistical methods reveals oncogenic RNA cancer drivers`, Proceedings of the National Academy of Sciences, Jul 2019, 201900391; DOI: https://doi.org/10.1073/pnas.1900391116

# Contact

In case of any questions, please contact Roozbeh Dehghannasiri (rdehghan@stanford.edu) or Milos Jordanski (milos.jordanski@stanford.edu).
