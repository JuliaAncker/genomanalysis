#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 06:00:00
#SBATCH -J 01_genome_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load Flye

#Your commands

flye --pacbio-raw /home/juliaa/genomanalys/Data/Raw_data/WGS_pacbio_samples/SRR6037732_scaffold_10.fq.gz --out-dir /home/juliaa/genomanalys/Data/Assembly/genome_assembly/pac_bio_assembly --threads 4


