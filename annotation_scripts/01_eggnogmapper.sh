#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:20:00
#SBATCH -J eggnogmapper
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load eggNOG-mapper/2.1.9

input=/home/juliaa/genomanalys/Data/Assembly/genome_assembly/repeat_masker_assembly/pilon_assembly.fasta.masked
output_dir=/home/juliaa/genomanalys/Data/Annotation/eggmapper


emapper.py -m hmmer -i $input --itype CDS -o ${output_dir}/emapper1
