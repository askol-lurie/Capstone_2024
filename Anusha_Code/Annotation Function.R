library(readr)
library(dplyr)

# Replace 'path/to/your/file.tsv.gz' with the actual path
annotations <- read_tsv("AlphaMissense_hg38.tsv", skip=3)

# Check the column names
colnames(annotations)

# Look at which chromosomes are represented 
unique_chrs <- unique(annotations$'#CHROM')
print(unique_chrs) # Only chromosome 1 is represented 
print(unique(annotations$POS)) # Includes genes in positions 69094-69506

# Create a function for the desired annotation
extract_variants <- function(data, chr, start, end) {
# Filter the data based on the specified region
  filtered_data <- data %>%
    filter(`#CHROM` == chr & POS >= start & POS <= end)
  return(filtered_data)
}

# Example usage 
region_variants <- extract_variants(annotations, "chr1", 69094, 69415)
print(region_variants)

# github change
