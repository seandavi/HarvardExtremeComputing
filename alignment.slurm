#!/bin/bash
#
#SBATCH -n 8 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-8:00 # Runtime in D-HH:MM
#SBATCH -p general # Partition to submit to
#SBATCH --mem=8000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=seandavi@gmail.com # Email to which notifications will be sent

cd /n/regal/ac290r/$USER 
echo $1

# load up modules
echo "loading modules"
source new-modules.sh
module load bwa
module load samtools

# do the actual alignment using 8 cores (threads)
echo "aligning $1"
bwa mem -t 8 /n/regal/informatics_public/ref/igenome/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \
  /n/regal/ac290r/tmp/${1}_R1.fq.gz \
  /n/regal/ac290r/tmp/${1}_R2.fq.gz | samtools view -bS - | samtools sort - ${1}.sorted

# and index the bam file
echo "indexing $1"
samtools index $1.sorted.bam

