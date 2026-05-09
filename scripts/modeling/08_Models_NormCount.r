###############
###############
#Prep for community models
# Load required libraries
library(dplyr)
library(readr)

# 1. Read OTU table
otu_file <- "C:/Users/jonat/OneDrive - University of Glasgow/PFAS_Microbiome/Manuscript/Inputs/otu_table_clean.tsv"
otu_df <- read.delim(otu_file, check.names = FALSE, stringsAsFactors = FALSE)

# 2. Trim whitespace from species names
otu_df$Taxa <- trimws(otu_df$Taxa)

# 3. Set species names as row names and remove Taxa column
rownames(otu_df) <- otu_df$Taxa
otu_df$Taxa <- NULL

# 4. Convert all counts to numeric (in case read as character)
otu_df[] <- lapply(otu_df, function(x) as.numeric(as.character(x)))

# 5. Normalize each sample (column) to relative abundance
otu_rel <- sweep(otu_df, 2, colSums(otu_df), FUN = "/")

# 6. Write to CSV
out_path <- "C:/Users/jonat/OneDrive - University of Glasgow/PFAS_Microbiome/Manuscript/Inputs/relative_abundance.csv"
write.csv(otu_rel, file = out_path, row.names = TRUE)

cat("✅ Relative abundance table saved to:", out_path, "\n")

#I've now manually added GEM names, now need to clean
library(dplyr)
library(readr)

# 1. Read the normalized coverage table
norm_file <- "C:/Users/jonat/OneDrive - University of Glasgow/PFAS_Microbiome/Manuscript/Inputs/NormCoverage.csv"
df <- read.csv(norm_file, check.names = FALSE, stringsAsFactors = FALSE)

# 2. Ensure the first column is named "Species" (adjust if needed)
colnames(df)[1] <- "Species"

# 3. Sum duplicates by Species
df_summed <- df %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()

# 4. Re-normalize relative abundances for each sample (column) so that they sum to 1
sample_cols <- colnames(df_summed)[-1]  # exclude Species column
df_norm <- df_summed
df_norm[sample_cols] <- sweep(df_summed[sample_cols], 2, colSums(df_summed[sample_cols]), FUN = "/")

# 5. Write output
out_path <- "C:/Users/jonat/OneDrive - University of Glasgow/PFAS_Microbiome/Manuscript/Inputs/NormCoverage_Cleaned.csv"
write.csv(df_norm, out_path, row.names = FALSE)

cat("✅ Cleaned and normalized coverage table saved to:", out_path, "\n")

