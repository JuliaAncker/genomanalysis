#!/bin/bash -l

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -J 04_braker
#SBATCH --mail-type=ALL
#SBATCH --mail-user julia.ancker@gmail.com
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools
module load braker/2.1.1_Perl5.24.1
module load augustus/3.2.3_Perl5.24.1
module load bamtools/2.5.1
module load blast/2.9.0+
module load GenomeThreader/1.7.0
module load samtools/1.8
module load GeneMark/4.33-es_Perl5.24.1

source $AUGUSTUS_CONFIG_COPY
export AUGUSTUS_CONFIG_PATH=/home/juliaa/genomanalys/git_repository/annotation_scripts/augustus_config

chmod a+w -R /home/juliaa/genomanalys/git_repository/annotation_scripts/augustus_config/species/
cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key

export AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin
export AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts
export GENEMARK_PATH=/sw/bioinfo/GeneMark/4.33-es/snowy


path_genome=/home/juliaa/genomanalys/Data/Assembly/genome_assembly/repeat_masker_assembly/pilon_assembly.fasta.masked
bam_files=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/02_star_rp/mapping_rp
output_dir=/home/juliaa/genomanalys/Data/Annotation/braker
bam1=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/02_star_rp/mapping_rp/SRR6040092_scaffold_10.1.fastq.gzAligned.out.bam


bam_list=""

for bamfile in "$bam_files"/*.bam; do
	bam_list+=",$bamfile"
done

bam_list=${bam_list:1}	

braker.pl --species=durian  --useexisting --genome=$path_genome --bam=$bam_list --softmasking
