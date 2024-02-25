setwd("THE_LIST_of_DATAFRAMES_to_merge/")

library(dplyr)
library(tidyverse)

#####################################################################################
#####################################################################################

file_names=c("simple.the.samples.number.of.SNV.INDEL.strelka.tnscope.txt",
"simple.the.samples.number.of.SV.on.bedpe.files.SV.txt",
"simple.the.samples.results.CHORD.hr.status.txt",
"simple.the.samples.results.CHORD.hrd.type.txt",
"simple.the.samples.results.NUMBER.KATAEGIS.LOCI.default.settings.txt",
"simple.the.samples.results.p.hrd.txt",
"simple.the.samples.results.p_BRCA1.5.txt",
"simple.the.samples.results.p_BRCA1.50.txt",
"simple.the.samples.results.p_BRCA1.95.txt",
"simple.the.samples.results.p_BRCA1.txt",
"simple.the.samples.results.p_BRCA2.5.txt",
"simple.the.samples.results.p_BRCA2.50.txt",
"simple.the.samples.results.p_BRCA2.95.txt",
"simple.the.samples.results.p_BRCA2.txt",
"simple.the.samples.results.p_hrd.5.txt",
"simple.the.samples.results.p_hrd.50.txt",
"simple.the.samples.results.p_hrd.95.txt",
"simple.the.samples.results.remarks.hr.status.txt",
"simple.the.samples.tnscope.strelka2.hg38.snv.signatures.2015.signature.first.txt",
"simple.the.samples.tnscope.strelka2.hg38.snv.signatures.2015.signature.second.txt",
"simple.the.samples.tnscope.strelka2.hg38.snv.signatures.2015.signature.third.txt",
"simple.the_list_files_with_results_of_scarHRD_SCORES_on_SCLUST_RESULTS.txt",
"simple.the.samples.number.of.SNV.INDEL.strelka.tnscope.for.MSI.MS.CHORD.analysis.txt")

# if it is necessary, to shorten the name of the SAMPLES
# Function to shorten strings in the first column of a data frame

shorten_strings <- function(df) {
  df[[1]] <- substr(df[[1]], 1, 18)  # Shorten to the first 18 characters (adjust as needed)
  return(df)
}

# Function to apply the above function to each data frame in a list

shorten_strings_in_list <- function(df_list) {
  df_list_shortened <- lapply(df_list, shorten_strings)
  return(df_list_shortened)
}

file_list <- list.files(path = ".", pattern = "\\.txt", full.names = TRUE)
# for (file in file_list) { print(file) }

for (file in file_list) {
    
  # Read the text file into a data frame
  df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

  # Shorten strings in the first column
  df_shortened <- shorten_strings(df)

  # Save the modified data frame to a new text file
  output_file <- paste0("shortened_", basename(file))
  write.table(df_shortened, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

}

#####################################################################################
#####################################################################################
# to combine these files
list.files()
#####################################################################################
#####################################################################################

# Specify the prefix and suffix for file filtering
prefix <- "shortened"
suffix <- "txt"

# Create a regular expression pattern for filtering files
pattern <- sprintf("^%s.*%s$", prefix, suffix)

# List files that match the pattern
matching_files <- list.files(path = ".", pattern = pattern, full.names = TRUE)

# Print the matching files
print(matching_files)

#####################################################################################
file.list = matching_files
#####################################################################################

# Initialize an empty list to store data frames
df_list <- list()

# Loop over each text file and read into a data frame
for (file in file.list) {
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  df_list[[basename(file)]] <- df
}

head(df_list[[1]],2)
head(df_list[[22]],2)

list.df = df_list 

merged.df <- Reduce(function(x, y) merge(x, y, by = "SAMPLE", all = FALSE), df_list)
head(merged.df,2)
dim(merged.df)

write.table(merged.df, 
            file = "integrated.files.with.all.the.outputs.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


#####################################################################################
library("tidyverse")
#####################################################################################

merged_df <- reduce(df_list, full_join, by = "SAMPLE")

dim(merged_df)
head(merged_df,2)

write.table(merged_df, 
            file = "integrated.files.with.all.the.outputs.full.join.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#####################################################################################
#####################################################################################