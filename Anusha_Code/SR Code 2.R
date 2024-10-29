exonfile <- 'C:/Users/sramachandran/Documents/GitHub/Capstone_2024/Resources/CodingLocs_hg38.rds'
exonfile <- '/Users/anushaakhtar/Documents/GitHub/Capstone_2024/Resources/CodingLocs_hg38.rds'
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
clinvar_vcf <- "/Volumes/Seagate/Capstone/clinvar.vcf.gz"

# Load the VCF file
#change to 38?
vcf_clinvar <- readVcf(clinvar_vcf, "hg38") 

# Check the contents of the VCF file
vcf_clinvar
summary(vcf_clinvar)

# Check the header of the VCF file
header(vcf_clinvar)

# View the fixed fields (chromosome, position, reference, alternate alleles)
fixed(vcf_clinvar)

# Check the INFO fields (annotations like clinical significance, allele frequency, etc.)
info(vcf_clinvar)

# Subset file 
vcf_clinvar_subset <- vcf_clinvar[1:1000, ]

# Check the INFO fields (annotations like clinical significance, allele frequency, etc.) and convert to data frame
info_data <- info(vcf_clinvar_subset)
info_data
clinvar_info <- as.data.frame(info_data)


# Extract the genomic ranges and convert to data frame
genomic_ranges <- rowRanges(vcf_clinvar_subset)
clinvar_gr <- as.data.frame(genomic_ranges)

# Combine both data frames 
clinvar_combined_df <- cbind(clinvar_gr, clinvar_info)
head(clinvar_combined_df)


# Selecting specific columns 
# Should be able to merge with ref, alt, and start
clinvar_filtered <- clinvar_combined_df %>%
  select(start, REF, ALT, AF_ESP, ALLELEID, CLNHGVS, CLNSIG, MC, GENEINFO, RS)

# : seperated list with clinsigs
# include any variables you think are related to pathogenicity 
# seq names might be different between gnomad and clinvar
# keep p. if you see it 




# Just try on CHR 22 for workflow 
## Path to the gnomAD VCF for chromosome 22
gnomad_vcf_chr22 <- "F:/Capstone/Resources/gnomAD/gnomad.joint.v4.1.sites.chr22.vcf.bgz"

rng <- GRanges(seqnames="chr22", ranges=IRanges(
  start=c(50301422, 50989541), 
  end=c(50312106, 51001328),
  names=c("gene_79087", "gene_644186")))

tab <- TabixFile(gnomad_vcf_chr22)
vcf_rng <- readVcf(tab, "hg38", param=rng)


# We would want to exlude everything in the filter column that doesn't say PASS
info(vcf_rng)

header(vcf_rng)

# row ranges to find the position and other fixed things we will need to keep and the filter info 
head(rowRanges(vcf_rng), 3)

# 
A <- info(header(vcf_rng))
write.table(file="header.txt",A, quote= FALSE, sep="\t")

# Load the gnomAD VCF file for chromosome 22
vcf_gnomad_chr22 <- readVcf(gnomad_vcf_chr22, "hg38")

# Check the contents
head(vcf_gnomad_chr22)

# Look for the what the fixed fields are 
fixed_fields <- fixed(vcf_gnomad_chr22)

# Filter the gnomAD VCF for this region
filtered_gnomad_chr22 <- subsetByOverlaps(vcf_gnomad_chr22, regions)

# Filter the ClinVar VCF for the same region
filtered_clinvar <- subsetByOverlaps(vcf_clinvar, regions)

# Convert filtered gnomAD and ClinVar VCFs to data frames
gnomad_df <- as.data.frame(filtered_gnomad_chr22)
clinvar_df <- as.data.frame(filtered_clinvar)

### Oct 25 
# Look for fixed fields within GnomAD
fixed_fields <- fixed(vcf_rng)
head(fixed_fields)

# Check the info fields 
info_fields <- info(header(vcf_rng))
info_fields

# Since we couldnt not fine CHR or POS in fixed file - use row ranges 
head(rowRanges(vcf_rng))

# Convert rowRanges to a data frame to extract CHROM and POS
chrom_pos_df <- as.data.frame(rowRanges(vcf_rng))[, c("seqnames", "start", "REF", "ALT")]

# Rename the columns so POS and CHROM are clear 
chrom_pos_df <- chrom_pos_df %>%
  rename(
    CHROM = seqnames,
    POS = start
  )

# Review DF 
head(chrom_pos_df)

# Extract the FAFmax_faf95_max_joint column from the INFO fields
fafmax <- info(vcf_rng)$FAFmax_faf95_max_joint

# Add the INFO field to chrom_pos_df
chrom_pos_df$FAFmax_faf95_max_joint <- fafmax

#View updated DF
head(chrom_pos_df)

# Just try on CHR 22 for workflow 
## Path to the gnomAD VCF for chromosome 22
gnomad_vcf_chr22 <- "F:/Capstone/Resources/gnomAD/gnomad.joint.v4.1.sites.chr22.vcf.bgz"
gnomad_vcf_chr22 <- "/Volumes/Seagate/Capstone/gnomad.joint.v4.1.sites.chr22.vcf.bgz"

rng <- GRanges(seqnames="chr22", ranges=IRanges(
  start=c(50301422, 50989541), 
  end=c(50312106, 51001328),
  names=c("gene_79087", "gene_644186")))

tab <- TabixFile(gnomad_vcf_chr22)
vcf_rng <- readVcf(tab, "hg38", param=rng)


# We would want to exlude everything in the filter column that doesn't say PASS
info(vcf_rng)

header(vcf_rng)

# row ranges to find the position and other fixed things we will need to keep and the filter info 
head(rowRanges(vcf_rng), 3)

# 
A <- info(header(vcf_rng))
write.table(file="header.txt",A, quote= FALSE, sep="\t")
  
# Load the gnomAD VCF file for chromosome 22
vcf_gnomad_chr22 <- readVcf(gnomad_vcf_chr22, "hg38")

# Check the contents
head(vcf_gnomad_chr22)

# Look for the what the fixed fields are 
fixed_fields <- fixed(vcf_gnomad_chr22)

# Filter the gnomAD VCF for this region
filtered_gnomad_chr22 <- subsetByOverlaps(vcf_gnomad_chr22, regions)

# Filter the ClinVar VCF for the same region
filtered_clinvar <- subsetByOverlaps(vcf_clinvar, regions)

# Convert filtered gnomAD and ClinVar VCFs to data frames
gnomad_df <- as.data.frame(filtered_gnomad_chr22)
clinvar_df <- as.data.frame(filtered_clinvar)

### Oct 25 
# Look for fixed fields within GnomAD
fixed_fields <- fixed(vcf_rng)
head(fixed_fields)

# Check the info fields 
info_fields <- info(header(vcf_rng))
info_fields

# Since we couldnt not fine CHR or POS in fixed file - use row ranges 
head(rowRanges(vcf_rng))

# Convert rowRanges to a data frame to extract CHROM and POS
chrom_pos_df <- as.data.frame(rowRanges(vcf_rng))[, c("seqnames", "start", "REF", "ALT")]

# Rename the columns so POS and CHROM are clear 
chrom_pos_df <- chrom_pos_df %>%
  rename(
    CHROM = seqnames,
    POS = start
  )

# Review DF 
head(chrom_pos_df)

# Extract the FAFmax_faf95_max_joint column from the INFO fields
fafmax <- info(vcf_rng)$FAFmax_faf95_max_joint

# Add the INFO field to chrom_pos_df
chrom_pos_df$FAFmax_faf95_max_joint <- fafmax

#View updated DF
head(chrom_pos_df)

## unable the add the FAFMAX

## Do something similar for ClinVar 

#Specify range
rng <- GRanges(seqnames="chr22", ranges=IRanges(
  start=c(50301422, 50989541), 
  end=c(50312106, 51001328),
  names=c("gene_79087", "gene_644186")))

# Use TabixFile to read only the specified region
clinvar_tab <- TabixFile(clinvar_vcf)
clinvar_rng <- readVcf(clinvar_tab, "hg38", param=rn
                       
# Need to check for REF, ALT, POS, CHROM and make that into a df for clinvar 


