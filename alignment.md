# Alignment on odyssey

## Interactive

```bash
# Start an interactive job with 8 cores (parallelization) and 8GB of RAM for 6 hours
srun -p interact --pty -c 8 --mem 8000 -t 0-6:00 bash

# load up modules
source new-modules.sh
module load bwa
module load samtools

# do the actual alignment
bwa mem -t 4 /n/regal/informatics_public/ref/igenome/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome \
  /n/regal/ac290r/fastq/NCI-H460_R1.fq.gz \
  /n/regal/ac290r/fastq/NCI-H460_R2.fq.gz > NCI-H460.sam

# convert from text-based SAM format to binary (and compressed) BAM format
samtools view -bS NCI-H460.sam > NCI-H460.bam

#sort into chromosome order
samtools sort NCI-H460.bam NCI-H460.sorted

# and index the bam file
samtools index NCI-H460.sorted.bam

#cleanup 
rm NCI-H460.sam
rm NCI-H460.bam
```
