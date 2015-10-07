# HarvardExtremeComputing

Materials for the Fall 2015 Harvard Extreme Computing course.

## Lecture 1

This lecture will serve as a basic review of the 
core concepts of normal biology and also to introduce
students to technologies that allow researchers to probe
the genome to better understand normal and human biology.

[Slides](http://watson.nci.nih.gov/~sdavis/courses/HarvardExtremeComputing/BasicBio.pdf)

## Lab 1

In this lab, the idea is to get familiar with the 
[R statistical programming environment](http://www.r-project.org). To get
started:

1. Install the [R software](https://cran.r-project.org/).
2. Consider installing [RStudio](https://www.rstudio.com/products/rstudio/download/) as a convenient environment for accessing R

I have created some introductor materials here. The "lecture slides" are meant to be used "interactively" and
additional material is useful to get a sense of some of the basic capabilities of R.

- [Introduction to R](http://watson.nci.nih.gov/~sdavis/tutorials/IntroToR/)


## Lecture 2

We will continue our discussion of biological principles and begin
to delve more deeply into technologies that allow us to 
examine normal and disease biology at a molecular level. We will
also begin to discuss data analysis approaches to high-throughput
biological data analysis.

## Lab 2

In this lab, we will work with [the Bioconductor project](http://bioconductor.org) as
an environment for biological data analysis. We will introduce some useful data containers
and look at our first cancer biology dataset involving gene copy number.

## Lecture 3

In this lecture, we will focus on a particularly important biological
application, RNA gene expression analysis.  We will work interactively
to analyze an RNA-seq dataset and you can continue to work on this
on your own.

## Lab 3

In this lab, we will introduce our "final project". We will be using data from 
the TCGA glioblastoma multiforme (GBM) project and some basic machine learning
to describe relationships between gene expression, gene copy number, and 
DNA methylation. The goal of this lab is to access and prepare our data prior
to applying large-scale, embarrasingly parallel computation.

- [Regulome Explorer](http://explorer.cancerregulome.org/all_pairs/)
- [Youtube video describing project](https://youtu.be/tPtJd6AzU8c?t=40m24s)
- [Behind the Compute Engine demo at Google I/O 2012 Keynote](https://cloud.google.com/compute/io)

## Lecture 4

In this lecture, we will discuss [the RF-ACE algorithm](https://code.google.com/p/rf-ace/) and 
show example data output from applying this approach to our TCGA GBM data.

## Lab 4

In this lab, we will be setting up and running the large-scale RF-ACE analysis. At this point,
students 

### Background materials

- Public resources
    - [The Cancer Genome Atlas (TCGA) project](http://cancergenome.nih.gov/)
    - [The National Center for Biomedical Informatics](http://www.ncbi.nlm.nih.gov/)
        - [PUBMED](http://www.ncbi.nlm.nih.gov/pubmed/)
        - [Entrez Gene](http://www.ncbi.nlm.nih.gov/gene/)
    - [The UCSC Genome Browser](http://genome.ucsc.edu)
    - [BioStars](https://biostars.org)
    - [Bioconductor](http://bioconductor.org)
    
- Genomic assays
    - [DNA microarrays]()
    - [Illumina sequencing overview (pdf)](https://www.illumina.com/documents/products/techspotlights/techspotlight_sequencing.pdf)
- Classic papers
    - [Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring](http://rileylab.bio.umb.edu/sites/g/files/g1314676/f/201502/Golub1999Molecular.pdf)
    - [Classification and diagnostic prediction of cancers using gene expression profiling and artificial neural networks](http://www.nature.com/nm/journal/v7/n6/pdf/nm0601_673.pdf)
    - [Comprehensive genomic characterization defines human glioblastoma genes and core pathways](http://www.nature.com/nature/journal/v455/n7216/full/nature07385.html)
    - [Gene expression profiling predicts clinical outcome of breast cancer](http://www.nature.com/nature/journal/v415/n6871/full/415530a.html)
    - [Genome-scale hypothesis scanning](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0000015)