# setup for dada2 microbiome lectures
# these are the homework 6 tasks

# make sure latest version of bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()

# install dada2
# see info at http://benjjneb.github.io/dada2/index.html
#source("https://bioconductor.org/biocLite.R")
biocLite("dada2")

library(ggplot2)

# me side note - all packages updated
# I had to reinstall Rcpp and digest - these are now OK
sessionInfo()

# see dada2 walkthrough
# http://benjjneb.github.io/dada2/tutorial.html 

library(dada2)
packageVersion("dada2")
## Loading required package: Rcpp
## Warning messages:
##  1: package ‘dada2’ was built under R version 3.3.3 
##  2: package ‘Rcpp’ was built under R version 3.3.3 
## [1] ‘1.2.2’

library(ShortRead)
packageVersion("ShortRead")
## [1] '1.32.1'

library(ggplot2)
packageVersion("ggplot2")
## [1] ‘2.2.1’

# downloaded the Mothur MiSeq SOP ZIP file
# "MiSeqSOPData.zip" 36MB file
# unzipped the files
# see C:\MyGithub\N741hmwk6\MiSeqSOPData\MiSeq_SOP 
# 45 files 163 MB - lots of *.fastq files

# also downloaded rdp_train_set_14.fa.gz
# and rdp_species_assignment_14.fa.gz
# and gg_13_8_train_set_97.fa.gz
# from https://zenodo.org/record/158955#.WNAsyvnyut9

#path <- "~/MiSeq_SOP" 
# note my current directory is
getwd()

# set path on Windows 10 system
path <- "./MiSeqSOPData/MiSeq_SOP"
# CHANGE ME to the directory containing the fastq files after unzipping.

fns <- list.files(path)
fns



