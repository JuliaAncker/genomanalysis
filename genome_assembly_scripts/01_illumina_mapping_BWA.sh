#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J 06_BWA_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules

module load bioinfo-tools
module load htslib/1.19
module load bwa/0.7.17
module load samtools

bwa index /home/juliaa/genomanalys/Data/Assembly/genome_assembly/pac_bio_assembly/assembly.fasta

bwa mem -t 2  /home/juliaa/genomanalys/Data/Assembly/genome_assembly/pac_bio_assembly/assembly.fasta /home/juliaa/genomanalys/Data/Raw_data/WGS_illumina_samples/SRR6058604_scaffold_10.1P.fastq.gz /home/juliaa/genomanalys/Data/Raw_data/WGS_illumina_samples/SRR6058604_scaffold_10.2P.fastq.gz  > /home/juliaa/genomanalys/Data/Assembly/genome_assembly/illumina_mapped/02_illumina_mapped.sam

samtools view -bS /home/juliaa/genomanalys/Data/Assembly/genome_assembly/illumina_mapped/02_illumina_mapped.sam \
    > /home/juliaa/genomanalys/Data/Assembly/genome_assembly/illumina_mapped/02_illumina_mapped.bam

# Sort BAM file
samtools sort /home/juliaa/genomanalys/Data/Assembly/genome_assembly/illumina_mapped/02_illumina_mapped.bam \
    -o /home/juliaa/genomanalys/Data/Assembly/genome_assembly/illumina_mapped/02_illumina_mapped_sorted.bam

# Index the sorted BAM file
samtools index /home/juliaa/genomanalys/Data/Assembly/genome_assembly/illumina_mapped/02_illumina_mapped_sorted.bam
