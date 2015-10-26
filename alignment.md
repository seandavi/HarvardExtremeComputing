# Alignment on odyssey

## setup

You should set up your own directory for capturing results.

```bash
cd /n/regal/ac290r/
mkdir $USER
cd $USER
```

## Interactive

```bash
# Start an interactive job with 8 cores (parallelization) and 8GB of RAM for 6 hours
srun -p interact --pty -c 8 --mem 8000 -t 0-6:00 bash

# switch to your own directory
cd /n/regal/ac290r/$USER

# load up modules
source new-modules.sh
module load bwa
module load samtools

# do the actual alignment using 8 cores (threads)
bwa mem -t 8 /n/regal/informatics_public/ref/igenome/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \
  /n/regal/ac290r/tmp/SHORT_TEST_R1.fq.gz \
  /n/regal/ac290r/fastq/SHORT_TEST_R2.fq.gz > TEST.sam

# convert from text-based SAM format to binary (and compressed) BAM format
samtools view -bS TEST.sam > TEST.bam

#sort into chromosome order
samtools sort TEST.bam TEST.sorted

# and index the bam file
samtools index TEST.sorted.bam

#cleanup 
rm NCI-H460.sam
rm NCI-H460.bam
```
