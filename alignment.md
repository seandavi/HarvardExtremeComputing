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
  /n/regal/ac290r/tmp/SHORT_TEST_R2.fq.gz > TEST.sam

# convert from text-based SAM format to binary (and compressed) BAM format
samtools view -bS TEST.sam > TEST.bam

#sort into chromosome order
samtools sort TEST.bam TEST.sorted

# and index the bam file
samtools index TEST.sorted.bam

#cleanup 
rm TEST.sam
rm TEST.bam
```

## Batch

Now, create a slurm batch script for aligning all the "pairs" of files in the `tmp` directory. The goal is to produce
sorted `.bam` files for each of the 61 samples.

## Variant calling

After you have generated all the sorted, indexed BAM files, it is time to call variants. Find the help for the `samtools` package (online) and perform basic variant calling across ALL the samples. 

## Bonus

A lot of computational problems involve stepwise processing of data where tasks rely on previous tasks; they form a directed acyclic graph. Consider converting these couple of steps into a workflow that can be run reproducibly and in parallel while. See here for some ideas of frameworks that help with this problem.

https://www.biostars.org/p/115745/

