#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J 03_01_quast
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load python/3.7.2
module load quast/5.0.2

pilon_assembly="/home/juliaa/genomanalys/Data/Assembly/genome_assembly/pilon_assembly/pilon_assembly.fasta"
output_dir="/home/juliaa/genomanalys/Data/Assembly/assembly_QC/QC_3"

quast.py $pilon_assembly -o $output_dir


