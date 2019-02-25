library(magrittr)
library(tidyverse)
library(CAGEfightR)
library(GenomicFeatures)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm9)

# Script wide settings
register(MulticoreParam(3))

# Rename for easier access
bsg <- BSgenome.Mmusculus.UCSC.mm9
odb <- org.Mm.eg.db

#### Reading the data ####

# Some example data
#design <- read.table("~/Desktop/mm9_nanotubes/design_matrix.txt", header=TRUE, stringsAsFactors=FALSE)
design <- read.table("~/BINF/CAGE_data/mm9_nanotubes/design_matrix.txt", header=TRUE, stringsAsFactors=FALSE)

# Better formated design
design$Samples <- gsub(x=design$Samples, pattern=".fastq", replacement="")
design$BigWigPlus <- paste0("mm9.", design$Samples, ".plus.bw")
design$BigWigMinus <- paste0("mm9.", design$Samples, ".minus.bw")
design$Class <- factor(design$Class)

# Remove technical replicates and inflamed
design <- subset(design, Name != "C548_2", select=-c(Samples))
#design <- subset(design, Class == "control")

# Sort and name
design <- design[order(design$Class, design$Name),]
rownames(design) <- design$Name
design <- DataFrame(design)

# Set up BigWig files
bw_plus <- BigWigFileList(file.path("~/BINF/CAGE_data/mm9_nanotubes", design$BigWigPlus))
bw_minus <- BigWigFileList(file.path("~/BINF/CAGE_data/mm9_nanotubes", design$BigWigMinus))
names(bw_plus) <- rownames(design)
names(bw_minus) <- rownames(design)

#### Quantifying CAGE TSSs (CTSSs) and finding clusters ####

# Load CTSSs
CTSSs <- quantifyCTSSs(bw_plus, bw_minus, genome = seqinfo(bsg), design = design); invisible(gc())

# Calculate pooled CTSS
CTSSs <- calcTPM(CTSSs, outputColumn = "CTSSTags") %>%
    calcPooled()

# Default clustering
TCs <- quickTSSs(CTSSs)
BCs <- quickEnhancers(CTSSs)

#### Save ####

usethis::use_data(TCs, internal = TRUE)


