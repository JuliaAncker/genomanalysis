#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:20:00
#SBATCH -J 01_trim
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load trimmomatic/0.39

input_forward=/home/juliaa/genomanalys/Data/Raw_data/RNA_samples/untrimmed/SRR6040095_scaffold_10.1.fastq.gz
input_reverse=/home/juliaa/genomanalys/Data/Raw_data/RNA_samples/untrimmed/SRR6040095_scaffold_10.2.fastq.gz
output=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/rna_qc/trimmed
path_adapters=/home/juliaa/genomanalys/Data/Raw_data/adapters.fa

# Your commands

trimmomatic PE \
    $input_forward $input_reverse\
    ${output}/02_forward_paired_output.fastq.gz ${output}/02_forward_unpaired_output.fastq.gz \
    ${output}/02_reverse_paired_output.fastq.gz ${output}/02_reverse_unpaired_output.fastq.gz \
    ILLUMINACLIP:adapters/TruSeq2-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 MINLEN:36

