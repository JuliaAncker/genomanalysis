#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -J 02_star_index
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load star/2.7.11a

output_index=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/star_rp/star_index_rp
path_rp=/home/juliaa/genomanalys/Data/Assembly/genome_assembly/repeat_masker_assembly/pilon_assembly.fasta.masked


star --runThreadN 8 --runMode genomeGenerate --genomeDir $output_index --genomeFastaFiles $path_rp
