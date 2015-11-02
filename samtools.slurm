#!/bin/bash
#
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-24:00 # Runtime in D-HH:MM
#SBATCH -p serial_requeue # Partition to submit to
#SBATCH --mem=4000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=seandavi@gmail.com # Email to which notifications will be sent

cd /n/regal/ac290r/$USER 
source new-modules.sh
module load samtools
module load bcftools
samtools mpileup --reg ${1} -ugf ../genome.fa *sorted.bam \
  | bcftools call -vmO z -o nci60.${1}.vcf.gz