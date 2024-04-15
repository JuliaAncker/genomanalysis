#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -J 02_star_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load star/2.7.11a
module load samtools/0.1.19

path_index=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/STAR/index
path_to_run1=/home/juliaa/genomanalys/Data/Raw_data/RNA_samples/trimmed/forward 
path_to_run2=/home/juliaa/genomanalys/Data/Raw_data/RNA_samples/trimmed/reverse
output_dir=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/STAR/mapping_02


forward=""
reverse=""

for file in "$path_to_run1"/*; do
    forward_files+=("$file")
done

for file in "$path_to_run2"/*; do
    reverse_files+=("$file")
done


for ((x=0; x<${#forward_files[@]}; x++)); do
	forward_file="${forward_files[x]}"
	reverse_file="${reverse_files[x]}"
	file_name=$(basename "$forward_file")
	star --runThreadN 8 --genomeDir $path_index --readFilesCommand zcat --readFilesIn "$forward_file" "$reverse_file" --outFileNamePrefix $output_dir/$file_name
	samtools view -bS -o "$output_dir/${forward_file_name%.fastq.gz}Aligned.out.sam.bam" "$output_dir/${forward_file_name%.fastq.gz}Aligned.out.sam"
	rm "$output_dir/${forward_file_name%.fastq.gz}Aligned.out.sam"
done 

