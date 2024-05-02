#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:20:00
#SBATCH -J checkbam
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load samtools/0.1.19

# Your commands

samtools flagstat /home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/star_rp/mapping_rp/SRR6040092_scaffold_10.1.fastq.gzAligned.out.bam
