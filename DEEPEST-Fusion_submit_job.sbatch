#!/bin/bash
#SBATCH --job-name=DEEPEST-Fusion
#SBATCH --output=DEEPEST-Fusion.%j.out
#SBATCH --error=DEEPEST-Fusion.%j.err
#SBATCH --time=24:00:00
#SBATCH --account=horence
#SBATCH --cpus-per-task=32
#SBATCH --mem=256GB
#SBATCH --export=RABIX.CORES=32,RABIX.MEMORY=262144,ALL
#SBATCH --mail-type=BEGIN,END,FAIL

module purge
module load anaconda
module load scipy
module load bowtie/2.2.9
module load biopython/1.70
module load java
module load R
module load python/3.4
module load trim_galore/0.4.4
Rscript $HOME/packages.R
pip install --user pandas
pip install --user docopt
pip3 install --user docopt
pip3 install --user typing
java -jar /srv/gsfs0/projects/salzman/jorda/bunny-1.0.5-1/rabix-cli/trim_galore_new/rabix-cli-1.0.5.jar --no-container --verbose /scratch/DEEPEST-Fusion_pipeline.json /scratch/DEEPEST-Fusion_input.json
