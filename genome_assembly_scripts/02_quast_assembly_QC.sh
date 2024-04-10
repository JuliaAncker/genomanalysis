#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J 01_02_quast
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load python/3.7.2
module load quast/5.0.2

pilon_assembly="/home/juliaa/genomanalys/Data/Assembly/genome_assembly/pilon_assembly/pilon_assembly.fasta"
flye_assembly_ref="/home/juliaa/genomanalys/Data/Assembly/genome_assembly/pac_bio_assembly/assembly.fasta"
pacbio_lr="/home/juliaa/genomanalys/Data/Raw_data/WGS_pacbio_samples/SRR6037732_scaffold_10.fq.gz"
illumina_1="/home/juliaa/genomanalys/Data/Raw_data/WGS_illumina_samples/SRR6058604_scaffold_10.1P.fastq.gz"
illumina_2="/home/juliaa/genomanalys/Data/Raw_data/WGS_illumina_samples/SRR6058604_scaffold_10.2P.fastq.gz"
illumina_map_bam="/home/juliaa/genomanalys/Data/Assembly/genome_assembly/illumina_mapped/02_illumina_mapped_sorted.bam"
output_dir="/home/juliaa/genomanalys/Data/Assembly/assembly_QC"

quast.py "${pilon_assembly}" -r "${flye_assembly_ref}" --pacbio "${pacbio_lr}" --pe1 "${illumina_1}" --pe2 "${illumina_2}" --bam "${illumina_map_bam}" -o "${output_dir}"


# Your commands
