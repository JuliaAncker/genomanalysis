#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -J eggnogmapper
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load eggNOG-mapper/2.1.9
# module load gffread/0.12.7

ref_genome=/home/juliaa/genomanalys/Data/Assembly/genome_assembly/repeat_masker_assembly/pilon_assembly.fasta.masked
genemark=/home/juliaa/genomanalys/Data/Annotation/braker_02/GeneMark-ET/genemark.gtf
output_dir=/home/juliaa/genomanalys/Data/Annotation/eggmapper


# gffread $genemark -g $ref_genome -y ${output_dir}/predicted_proteins.faa
emapper.py -i ${output_dir}/updated_predicted_proteins.faa -o ${output_dir}/emapper1 --override
