# HarvardExtremeComputing

Materials for the Fall 2015 Harvard Extreme Computing course.

- [Sean Davis](http://watson.nci.nih.gov/~sdavis/), [@seandavis12](https://twitter.com/seandavis12)
- Eric Stahlberg

Course email list (approval, but just ask): https://groups.google.com/forum/#!forum/harvardextremecomputing2015


## Lecture 1

This lecture will serve as a basic review of the 
core concepts of normal biology with a focus on the
current state of the understanding of the human
genome.

- [Slides](http://1.usa.gov/1RHVnEf)
- [Webinar link](http://webmeeting.nih.gov/bioifx/)
- [UCSC Exercises](http://www.hsl.unc.edu/collections/Bioinformatics/ClassMaterials/UCSC_exercises_V14a_unc.pdf)

## Lab 1

This lab will be run by the Research Computing group at Harvard and will cover:

1. Introduction to the R programming language
2. Accessing and using the Odyssey High Performance Computing (HPC) cluster

### Resources of interest

- [R statistical programming environment](http://www.r-project.org).
- [Introduction to Harvard's Odyssey cluster](https://rc.fas.harvard.edu/training/intro-to-odyssey/)
- [Odyssey module list](https://rc.fas.harvard.edu/resources/module-list/)

I have created some introductory materials here that you may use (not required) 
to augment what you do during the lab. The "lecture slides" are meant to be used "interactively" and
additional material is useful to get a sense of some of the basic capabilities of R.

- [Introduction to R](http://watson.nci.nih.gov/~sdavis/tutorials/IntroToR/)

### Getting started with R

1. Install the [R software](https://cran.r-project.org/).
2. Consider installing [RStudio](https://www.rstudio.com/products/rstudio/download/) as a convenient environment for accessing R




## Lecture 2

We will continue our discussion of biological principles and begin
to delve more deeply into technologies that allow us to 
examine normal and disease biology at a molecular level. We will
also begin to discuss data analysis approaches to high-throughput
biological data analysis.

- [Slides](http://1.usa.gov/1GRMzoY)
- [Class exercises](http://watson.nci.nih.gov/~sdavis/tutorials/IntroToR/)

## Lab 2

### Section 1

- Gene expression microarray data exercise
- [Class exercises](http://watson.nci.nih.gov/~sdavis/tutorials/IntroToR/)

### Section 2

- Begin working on DNA alignment project, based on [this paper](http://cancerres.aacrjournals.org/content/73/14/4372.full).
- [Starter materials](alignment.md)

## Lecture 3

This will be a combined lecture/lab.  Topics covered include:

- Basics of parallel computing (1 hour)
- Work on parallelizing and running full alignment dataset
- Develop template for bam --> vcf (variant calling)

## Lab 3

- Perform/review variant calling
- Visualize results in IGV
   - http://watson.nci.nih.gov/projects/nci60/wes/VCF/SK-MEL-2.vcf
   - http://watson.nci.nih.gov/projects/nci60/wes/BAMS/SK-MEL-2_reord_mdups_ralgn_fmate_recal.bam
- [Drug sensitivity prediction](DrugResponsePrediction.Rmd)

## Lecture 4

- Wrapup of final project, including unanswered and open questions
- Additional topics in Extreme Biocomputing
- Parting questions, comments

## Background materials

### Biology

- Genetics and gene regulation
    - [GEUVADIS paper](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3918453/)

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
- Manuscripts of interest and classic papers
    - [The cancer genome](http://www.nature.com/nature/journal/v458/n7239/full/nature07943.html)
    - [Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression Monitoring](http://rileylab.bio.umb.edu/sites/g/files/g1314676/f/201502/Golub1999Molecular.pdf)
    - [Classification and diagnostic prediction of cancers using gene expression profiling and artificial neural networks](http://www.nature.com/nm/journal/v7/n6/pdf/nm0601_673.pdf)
    - [Comprehensive genomic characterization defines human glioblastoma genes and core pathways](http://www.nature.com/nature/journal/v455/n7216/full/nature07385.html)
    - [Gene expression profiling predicts clinical outcome of breast cancer](http://www.nature.com/nature/journal/v415/n6871/full/415530a.html)
    - [Genome-scale hypothesis scanning](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0000015)

### Computation and machine learning

- [bwa aligner manuscript](http://www.ncbi.nlm.nih.gov/pubmed/19451168)
- [Bioconductor](https://bioconductor.org/)
