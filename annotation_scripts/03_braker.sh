#!/bin/bash -l
#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J BRAKER_annotation
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julia.ancker@gmail.com
#SBATCH --output=BRAKER_annotation.%j.out

# Load required modules
module load bioinfo-tools
module load braker/2.1.1_Perl5.24.1
module load augustus/3.2.3_Perl5.24.1
module load bamtools/2.5.1
module load blast/2.9.0+
module load GenomeThreader/1.7.0
module load samtools/1.8
module load GeneMark/4.33-es_Perl5.24.1

OUT_DIR=/home/juliaa/genomanalys/Data/Annotation/braker_02

# Set environment variables
export AUGUSTUS_CONFIG_PATH=/home/juliaa/genomanalys/Data/Annotation/braker/augustus_config
export AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin
export AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts
export GENEMARK_PATH=/sw/bioinfo/GeneMark/4.33-es/snowy

# Update AUGUSTUS configuration permissions
mkdir -p ${OUT_DIR}/augustus_config/species
chmod a+w -R ${OUT_DIR}/augustus_config/species

# Copy GeneMark key
cp -vf /sw/bioinfo/GeneMark/4.33-es/snowy/gm_key $HOME/.gm_key

# Copy AUGUSTUS configuration file
source $AUGUSTUS_CONFIG_COPY
echo "AUGUSTUS_CONFIG_COPY is set to: $AUGUSTUS_CONFIG_COPY"


# Replace the following paths with the correct paths to your input files
GENOME=/home/juliaa/genomanalys/Data/Assembly/genome_assembly/repeat_masker_assembly/pilon_assembly.fasta.masked

bam_files=/home/juliaa/genomanalys/Data/Assembly/transcriptome_assembly/star_rp/mapping_rp

bam_list=""

for bamfile in "$bam_files"/*.bam; do
	bam_list+=",$bamfile"
done

bam_list=${bam_list:1}	

# Running BRAKER
braker.pl \
   --genome=${GENOME} \
   --bam=${bam_list} \
   --species=your_species \
   --softmasking \
   --useexisting \
   --cores=8	\
   --workingdir=${OUT_DIR}
