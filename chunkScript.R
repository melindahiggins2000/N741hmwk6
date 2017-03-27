
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")

# Set the path to the data files
path <- "./MiSeqSOPData/MiSeq_SOP"
fileNames <- list.files(path)
fileNames


# Read in the names of the .fastq files

fastqFiles <- fileNames[grepl(".fastq$", fileNames)]

# Sort so that forward / reverse reads are in the same order

fastqFiles <- sort(fastqFiles)

# Get just the forward read files

fileNameForwards <- fastqFiles[grepl("_R1", fastqFiles)]

# Get just the reverse read files

fileNameReverses <- fastqFiles[grepl("_R2", fastqFiles)]

# Get the sample names, assuming file naming convention is: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(fileNameForwards, "_"), `[`, 1)

# Specify the full file path to the forward and reverse read.s

fileNameForwards <- file.path(path, fileNameForwards)
fileNameReverses <- file.path(path, fileNameReverses)

# Visualize the quality profile of the first two files containing forward reads

plotQualityProfile(fileNameForwards[[1]])
plotQualityProfile(fileNameForwards[[2]])

# Visualize the quality profile of the first two files containing reverse reads

plotQualityProfile(fileNameReverses[[1]])
plotQualityProfile(fileNameReverses[[2]])

# Make a directory and filenames for the filtered fastqs

filt.path <- file.path(path, "filtered")
if(!file_test("-d", filt.path)) dir.create(filt.path)
filtForwards <- file.path(filt.path, paste0(sample.names, "_F_filt.fastq.gz"))
filtReverses <- file.path(filt.path, paste0(sample.names, "_R_file.fastq.gz"))

# Now filter

for (i in seq_along(fileNameForwards)) {
  fastqPairedFilter(c(fileNameForwards[i], fileNameReverses[i]), 
                    c(filtForwards[i], filtReverses[i]),
                    truncLen=c(240,160),
                    maxN = 0, maxEE = c(2,2),
                    truncQ = 2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)
}

# Dereplicate

derepForwards <- derepFastq(filtForwards, verbose=TRUE)
derepReverses <- derepFastq(filtReverses, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepForwards) <- sample.names
names(derepReverses) <- sample.names

# note - this takes a while to run....
dadaForwards.learn <- dada(derepForwards, err=NULL, selfConsist = TRUE, multithread=TRUE)
dadaReverses.learn <- dada(derepReverses, err=NULL, selfConsist = TRUE, multithread=TRUE)


# Store initial error estimates for Forward reads

errForwards <- dadaForwards.learn[[1]]$err_out

# Now for the Reverse reads

errReverses <- dadaReverses.learn[[1]]$err_out

# Plot the estimated error rates for the Forward reads

plotErrors(dadaForwards.learn[[1]], nominalQ=TRUE)

# And for the Reverse reads

plotErrors(dadaReverses.learn[[1]], nominalQ = TRUE)

# First with the Forward reads

dadaForwards <- dada(derepForwards, err = errForwards, multithread = TRUE)

# Then with the Reverse reads

dadaReverses <- dada(derepReverses, err = errReverses, multithread = TRUE)

# Inspect the dada-class objects returned by the dada function

dadaForwards[[1]]
dadaReverses[[1]]


# Merge the denoised forward and reverse reads

mergedPairs <- mergePairs(dadaForwards, derepForwards, dadaReverses, derepReverses, verbose = TRUE )

# Inspect the merged data.frame from the first sample

head(mergedPairs[[1]])


# Construct sequence table

seqtab <- makeSequenceTable(mergedPairs[names(mergedPairs) != "Mock"])

# Consider the table

dim(seqtab)
class(seqtab)

# Inspect the distribution of sequence lengths

table(nchar(getSequences(seqtab)))


# Remove chimeric sequences

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab)



# Assign taxonomy

# First initialize random number generator for reproducibility

set.seed(100)

path
list.files(path)
reference_fasta <- "rdp_train_set_14.fa.gz"
taxa <- assignTaxonomy(seqtab.nochim, refFasta = reference_fasta)
unname(head(taxa))

# this takes a while to run...
# Assign species
genus.species <- assignSpecies(seqtab.nochim, "rdp_species_assignment_14.fa.gz")

# install DECIPHER
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DECIPHER")


library(DECIPHER)
seqs <- getSequences(seqtab.nochim)

# This next command will allow propagation of sequence names to the tip labels of the tree
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

# install phangorn from CRAN

library(phangorn)

# Construct the tree
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Tip order will not equal sequence order
fit <- pml(treeNJ, data=phang.align)

## negative edges length changed to 0.

fitGTR <- update(fit, k=4, inv=0.2)

# this next step takes a while to run...
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    rearrangement = "stochastic", control=pml.control(trace=0))
detach("package:phangorn", unload=TRUE)

