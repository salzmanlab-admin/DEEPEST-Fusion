# Building index files needed for DEEPEST-Fusion

DEEPEST-Fusion needs the following index files:

- Bowtie2 index files for the reference genome
- Bowtie2 index files for the regular junctions
- Bowtie2 index files for the scrambled junctions
- Bowtie2 index files for ribosome
- Bowtie2 index files for the transcriptome
- Bowtie2 index files for indel junctions (up to 5 symmetric indels at the splice site)
- Pickle files for genes/exons annotation
- Pickle file for the known fusions list (a list of known fusions constructed based on ChimerDB 3.0 curated list of known cancer fusions)
- Bowtie2 index files for known fusions

We provide two options for index files: 

   1- Download the tar folder containing all ready-to-use bowtie2 index files based on human genome build hg38 directly from (file size: 152 GB): https://drive.google.com/open?id=1ftT5UjKafcsDzUFlV49szOD8SiH3-Qws

   2- Download the reference sequence and annotation files for human genome build hg38 from: https://drive.google.com/a/stanford.edu/file/d/1PhZA63Jf4mHZr02pOR5wwadSPc2UKuNj/view?usp=sharing and follow the instructions below to build the index files:   
# Software Requirements

- Bowtie2 2.2.9
- Python 3.4 with Biopython/1.70 installed

# Creating known_fusions fasta and pickle files

The known_fusions.fa file should be used to build known fusions index. In order to create known_fusions.pickle, please run the following command:

    python create_pickle_file.py --fasta known_fusions.fa

# Creating regular/scrambled junction fasta files

In order to create regular/scrambled junction fasta files, please run the following commands:
    
    1. mkdir index && python makeExonDB.py -f /path/to/fasta_file -a path to /path/to/gtf_file -o /path/to/makeExonsDB_output_dir 
    2. sh createJunctionIndex.sh path/to/output_dir  /path/to/makeExonsDB_output_dir prefix_string
    
# Creating pickle_tar file for annotation

In order to create pickle_tar file, please run the following command:

    tar -czf reference_name.tar.gz reference_name
    
"reference_name" folder has to contain the following subfolders: exons, genes and records. These subfolders should contain the output files of the previous step


# Creating transcriptome file

In order to create transcriptome file, please use Tophat2 gtf_to_fasta tool:

    gtf_to_fasta /path/to/gtf_file /path/to/fasta_file gtf.fasta && awk '{if ($1 ~ /^>/) print ">"$2; else print $0}' gtf.fasta > transcripts.fa

# Creating indel junctions

In order to create indel junctions files, please run the following command:
    
    python subsample_reg_indels.py --input_file /path/to/reg_junctions.fa --probability p --reference_string reference_name
    
Parameter p represents the probability that regular junctions will be kept in indel junctions files. Since regular junctions fasta file for human genome is large, we have used p = 0.1

# Build bowtie2 indices

In order to build Bowtie2 index files for the reference genome, regular junctions, scrambled junctions, ribosome, transcriptome, indel junctions, known fusions, run the following command:

    bowtie2-build -f --threads numTreads ref.fa ./ref
 

