#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:20:00
#SBATCH -J 02_FastQC_wgs
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load FastQC/0.11.9
# Your commands
fastqc /home/juliaa/genomanalys/Data/Preprocessed_data/trimmed_data/WGS_samples_trimmed/SRR6058604_scaffold_10.2P.fastq.gz -o /home/juliaa/genomanalys/Data/Preprocessed_data/quality_checked_data 
