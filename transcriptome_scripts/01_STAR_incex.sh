#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -J 01_star_index
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load star/2.7.11a

output_index=/home/juliaa/genomanalys/Data/Annotation/STAR/index
path_pilon=/home/juliaa/genomanalys/Data/Assembly/genome_assembly/pilon_assembly/pilon_assembly.fasta


star --runThreadN 8 --runMode genomeGenerate --genomeDir $output_index --genomeFastaFiles $path_pilon
