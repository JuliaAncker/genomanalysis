#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:20:00
#SBATCH -J 01_sam2bam
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load samtools/0.1.19

# Your commands

input_dir=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/STAR/mapping_02

for sam_file in "$input_dir"/*.sam; do
	samtools view -bS "$sam_file" > "${sam_file%.sam}.bam"
	rm $sam_file
done
