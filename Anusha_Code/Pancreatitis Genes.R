library(dplyr) # will be used for data manipulation 
library(readr) # used to read in csv files
library(VariantAnnotation) # for variant annotation 
library(stringr) # split string to extract ref and alt allele


##ClinVar

# Combining all of ClinVar gene data into a single data frame
file_paths <- list.files(path = "ClinVar Files", pattern = "\\.txt$", full.names = TRUE)

# Reading each file into a data frame and storing them as a list 
data_frames_list <- lapply(file_paths, function(file_path) {
  # Read the data frame
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Convert the GRCh37Chromosome & Location column to character (change other columns as needed)
  df$GRCh37Chromosome <- as.character(df$GRCh37Chromosome)
  df$GRCh37Location <- as.character(df$GRCh37Location)
  df$GRCh38Chromosome <- as.character(df$GRCh38Chromosome)
  df$GRCh38Location <- as.character(df$GRCh38Location)
  df$AlleleID.s. <- as.character(df$AlleleID.s.)
  df$VariationID <- as.character(df$VariationID)
  
  # Return the modified data frame
  return(df)
})

# Using dplyr's bind_rows to combine each file into a single data frame
clinvar_combined_data_frame <- bind_rows(data_frames_list)

# Create a subset data frame keeping only relevant columns and reordering them
clinvar_subset <- clinvar_combined_data_frame %>% select(c(
  Gene.s.,
  Name,
  GRCh38Chromosome,
  GRCh38Location,
  VariationID,
  AlleleID.s.,
  dbSNP.ID,
  Germline.classification
))
View(clinvar_subset)  

# Display the chromosomes that are represented 
print(unique(clinvar_subset$GRCh38Chromosome)) # 1, 4, 5, 6, 7, 8, 11, 15





## GnomAD

# Combining all of GnomAD gene data into a single data frame using the same steps as above
# Unlike ClinVar, these files are a .csv file
file_paths <- list.files(path = "GnomAD Files", pattern = "\\.csv$", full.names = TRUE)
data_frames_list <- lapply(file_paths, function(file_path) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  return(df)
})

# Combine into a single data frame
gnomad_combined_data_frame <- bind_rows(data_frames_list)
print(gnomad_combined_data_frame)






# String split to extract ref and alt allele in ClinVar data set 
clinvar_subset$Transcript.Consequence <- sapply(strsplit(clinvar_subset$Name, ":"), function(x) {
  if (length(x) >= 2) {
    gsub("\\(.*\\)", "", x[2])
  } else {
    NA
  }
})
clinvar_subset$Prot_Tran.Consequence <- sapply(strsplit(clinvar_subset$Name, ":"), `[`, 2)

# Adding p. and c. information into one column in the GnomAD dataset 
gnomad_combined_data_frame <- gnomad_combined_data_frame %>%
  mutate(Prot_Tran.Consequence = ifelse(nchar(trimws(Protein.Consequence)) > 0,
                                        paste0(Transcript.Consequence, " (", trimws(Protein.Consequence), ")"),
                                        Transcript.Consequence))

head(gnomad_combined_data_frame$Prot_Tran.Consequence)
View(gnomad_combined_data_frame$Prot_Tran.Consequence)

# Checking structure of Prot_Tran.Consequence columns in both data sets to make sure they match
str(clinvar_subset$Prot_Tran.Consequence)
str(gnomad_combined_data_frame$Prot_Tran.Consequence) #they do match

# Getting rid of any trailing whitespace and capitalization issues 
clinvar_subset$Prot_Tron.Consequence <- trimws(clinvar_subset$Prot_Tron.Consequence)
gnomad_combined_data_frame$Prot_Tron.Consequence <- trimws(gnomad_combined_data_frame$Prot_Tron.Consequence)




# Comvine ClinVar and GnomAD data based on ref and alt allele and position
combined_df <- full_join(clinvar_subset, gnomad_combined_data_frame, by = "Prot_Tran.Consequence")


any(duplicated(gnomad_combined_data_frame$Prot_Tran.Consequence))
print(sum(duplicated(clinvar_subset$Prot_Tron.Consequence)))
print(sum(duplicated(gnomad_combined_data_frame$Prot_Tron.Consequence)))
## Annotate variants so that we can merge ClinVar and GnomAD data sets 
## Can combine them based on c.


# Need to first extract the ref and alt allele in the clinvar_subset dataset 
# Extract reference and alternate alleles for all variants in clinvar_subset
#clinvar_subset$ref <- str_extract(clinvar_subset$Name, "([A-Z])>([A-Z])")[, 1]
#clinvar_subset$alt <- str_extract(clinvar_subset$Name, "([A-Z])>([A-Z])")[, 2]

#extracted_values <- str_extract(clinvar_subset$Name, "([A-Z])>([A-Z])")
#print(extracted_values)

# Create a GRanges object from your data frame
#granges <- GRanges(
#  seqnames = your_data_frame$chromosome,
#  ranges = IRanges(start = your_data_frame$start, end = your_data_frame$end),
# ref = your_data_frame$ref,
#  alt = your_data_frame$alt
#)

# Annotate the variants
#annotated_granges <- annotateVariants(granges, TxDb("Ensembl", "hsapiens", version = "99"))







## Andrew's Example Code 
data <- tidyr::full_join(clinvar, gnomad, by=c("chr" = "chromosome", "position"="position", "ref"="ref", "alt"="alt")
                         
# Split ref and alt allele on clinvar bc it'll be needed for VEP 

# Next step get both files merged and pathogenicity defined based on clinvar and allele frequency info 









                         