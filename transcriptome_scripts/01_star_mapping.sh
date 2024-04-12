#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -J 01_star_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load star/2.7.11a

path_index=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/STAR/index
path_to_run1=/home/juliaa/genomanalys/Data/Raw_data/RNA_samples/trimmed/forward 
path_to_run2=/home/juliaa/genomanalys/Data/Raw_data/RNA_samples/trimmed/reverse 
output_dir=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/STAR/mapping

forward=""
reverse=""

for file in "$path_to_run1"/*; do
    forward+=",$file"
done

for file in "$path_to_run2"/*; do
    reverse+=",$file"
done

forward="${forward:1}"
reverse="${reverse:1}"

star --runThreadN 8 --genomeDir $path_index --readFilesCommand zcat --readFilesIn $forward $reverse --outFileNamePrefix $output_dir/01
