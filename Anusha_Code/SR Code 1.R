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

# Define Clinvar path to VCF 
clinvar_vcf <- "F:/Capstone/Resources/ClinVar/clinvar.vcf.gz"

# Load the VCF file
#change to 38?
vcf_clinvar <- readVcf(clinvar_vcf, "hg19") 

# Check the contents of the VCF file
vcf_clinvar
summary(vcf_clinvar)

# Check the header of the VCF file
header(vcf_clinvar)

# View the fixed fields (chromosome, position, reference, alternate alleles)
fixed(vcf_clinvar)

# Check the INFO fields (annotations like clinical significance, allele frequency, etc.)
info(vcf_clinvar)

# Just try on CHR 22 for workflow 
## Path to the gnomAD VCF for chromosome 22
gnomad_vcf_chr22 <- "F:/Capstone/Resources/gnomAD/gnomad.joint.v4.1.sites.chr22.vcf.bgz"

rng <- GRanges(seqnames="chr22", ranges=IRanges(
  start=c(50301422, 50989541), 
  end=c(50312106, 51001328),
  names=c("gene_79087", "gene_644186")))

tab <- TabixFile(gnomad_vcf_chr22)
vcf_rng <- readVcf(tab, "hg38", param=rng)

  
# Load the gnomAD VCF file for chromosome 22
vcf_gnomad_chr22 <- readVcf(gnomad_vcf_chr22, "hg38")


# Check the contents
head(vcf_gnomad_chr22)

#### Where things get messed up 

# Define areas of interest 
library(GenomicRanges)

# Define the region of interest on chromosome 22
regions <- GRanges(seqnames = c("chr22"), ranges = IRanges(start = c(100000), end = c(200000)))

# Filter the gnomAD VCF for this region
filtered_gnomad_chr22 <- subsetByOverlaps(vcf_gnomad_chr22, regions)

# Filter the ClinVar VCF for the same region
filtered_clinvar <- subsetByOverlaps(vcf_clinvar, regions)

# Convert filtered gnomAD and ClinVar VCFs to data frames
gnomad_df <- as.data.frame(filtered_gnomad_chr22)
clinvar_df <- as.data.frame(filtered_clinvar)

### Merge filtered data
# unable to open 

