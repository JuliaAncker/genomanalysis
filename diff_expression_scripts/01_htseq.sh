#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J 01_htseq
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load htseq/2.0.2

# Your commands

bam_files=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/02_star_rp/mapping_rp
genemark_file=/home/juliaa/genomanalys/Data/Annotation/braker_02/GeneMark-ET/genemark.gtf 
output_dir=/home/juliaa/genomanalys/Data/Differential_expression/htseq


for bamfile in "$bam_files"/*.bam; do

    base_name=$(basename "$bamfile" .bam)

    output_file="$output_dir/${base_name}_counts.txt"    

    htseq-count -f bam -s no "$bamfile" "$genemark_file" > "$output_file"
done
