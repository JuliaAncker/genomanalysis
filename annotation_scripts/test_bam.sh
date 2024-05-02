#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 00:20:00
#SBATCH -J 01_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools

bam_files=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/star_rp/mapping_rp


bam_list=""

for bamfile in "$bam_files"/*.bam; do
        bam_list+=",$bamfile"
done

bam_list=${bam_list:1}


echo "bam" $bam_list

