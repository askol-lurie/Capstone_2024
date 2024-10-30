## SET UP 

exonfile <- 'C:/Users/sramachandran/Documents/GitHub/Capstone_2024/Resources/CodingLocs_hg38.rds'
readRDS(exonfile)
exonpositions <- readRDS(exonfile)

# Install BiocManager to handle Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install VariantAnnotation for reading VCF files
BiocManager::install("VariantAnnotation")

# Install GenomicRanges for handling genomic regions
BiocManager::install("GenomicRanges")

# Install dplyr for data manipulation
install.packages("dplyr")

# Load the VariantAnnotation package for VCF file handling
library(VariantAnnotation)

# Load GenomicRanges for defining and handling regions of interest
library(GenomicRanges)

# Load dplyr for data manipulation
library(dplyr)

# Define region of interest
  # Specific for CHR 22
rng <- GRanges(seqnames="chr22", ranges=IRanges(
  start=c(50301422, 50989541), 
  end=c(50312106, 51001328),
  names=c("gene_79087", "gene_644186")))

# Define region of CHEK2: 
rng_chek2 <- GRanges(seqnames = "chr22", ranges = IRanges(
  start = 28636957,
  end = 28689186,
  names = "CHEK2"
))

# Path to gnomAD VCF for chromosome 22
gnomad_vcf_chr22 <- "F:/Capstone/Resources/gnomAD/gnomad.joint.v4.1.sites.chr22.vcf.bgz"

# Open the VCF with TabixFile for subsetting
tab_gnomad <- TabixFile(gnomad_vcf_chr22)
vcf_rng_gnomad <- readVcf(tab_gnomad, "hg38", param=rng)

# Try for just CHEK2 
vcf_rng_gnomad_chek2 <- readVcf(tab_gnomad, "hg38", param=rng_chek2)

# Extract fixed fields for CHROM, POS, REF, ALT
gnomad_fixed <- as.data.frame(rowRanges(vcf_rng_gnomad_chek2))

# Rename columns for clarity
gnomad_fixed <- gnomad_fixed %>%
  rename(
    CHROM = seqnames,
    POS = start
  )

# Check contents
head(gnomad_fixed)

### Having trouble addin this fafmax

##
info <- info(vcf_rng_gnomad_chek2)
faf <- unlist(info$fafmax_faf95_max_joint)

# Extract the FAFmax_faf95_max_joint allele frequency from INFO fields
gnomad_fixed$FAFmax_faf95_max_joint <- faf

# Check resulting gnomAD data frame
head(gnomad_fixed)

## Choose one of the large pop 
# Should NA be considered rare or is it number of samples
# how many sampels from each pop and what is the count 
# maybe get rid of small indels 







#######

# ClinVar Extraction 
clinvar_vcf <- "F:/Capstone/Resources/ClinVar/clinvar.vcf.gz"
vcf_clinvar <- readVcf(clinvar_vcf, "hg38")

# Open the VCF with TabixFile for subsetting
tab_clinvar <- TabixFile(clinvar_vcf)
vcf_rng_clinvar <- readVcf(tab_clinvar, "hg38", param=rng)

# Extract fixed fields for CHROM, POS, REF, ALT
clinvar_fixed <- as.data.frame(rowRanges(vcf_rng_clinvar))[, c("seqnames", "start", "REF", "ALT")]
