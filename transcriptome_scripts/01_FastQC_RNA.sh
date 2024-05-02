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

input=/home/juliaa/genomanalys/Data/Raw_data/RNA_samples/untrimmed/SRR6040095_scaffold_10.1.fastq.gz
output=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/rna_qc

fastqc $input -o $output
