#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J 05_pilon_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load Pilon/1.24


path_genome="/home/juliaa/genomanalys/Data/Assembly/genome_assembly/pac_bio_assembly/assembly.fasta"
path_bwa_bam="/home/juliaa/genomanalys/Data/Assembly/genome_assembly/illumina_mapped/02_illumina_mapped_sorted.bam"
path_bwa_bai="/home/juliaa/genomanalys/Data/Assembly/genome_assembly/illumina_mapped/02_illumina_mapped_sorted.bam.bai"

java -Xmx16G -jar "$PILON_HOME/pilon.jar" --genome "$path_genome" --frags "$path_bwa_bam" --output /home/juliaa/genomanalys/Data/Assembly/genome_assembly/pilon_assembly --changes




 
