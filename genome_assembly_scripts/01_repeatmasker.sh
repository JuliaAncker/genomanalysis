#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 02:00:00
#SBATCH -J 01_repeatmasker
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out

# Load modules
module load bioinfo-tools
module load RepeatMasker/4.1.5

pilon_assembly="/home/juliaa/genomanalys/Data/Assembly/genome_assembly/pilon_assembly/pilon_assembly.fasta"
output_dir="/home/juliaa/genomanalys/Data/Assembly/genome_assembly/repeat_masker_assembly"

RepeatMasker -pa 4 -dir "$output_dir" "$pilon_assembly"

