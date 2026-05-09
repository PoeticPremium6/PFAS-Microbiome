############################################################
# 02_suppfig_s2_family_and_elimination_boxplots.R
#
# Supplementary Figure S2: family-level abundance comparisons and PFAS elimination-rate boxplot.
############################################################

library(ggplot2)
library(dplyr)

# Define the output directory for saving plots
output_dir <- file.path("results", "supplementary_figures", "family_comparison")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Subset data for Case and Control groups
physeq_case <- subset_samples(physeq_rel_abund, Condition == "PFAS_Exposed")
physeq_control <- subset_samples(physeq_rel_abund, Condition == "Reference")

# Extract OTU and Taxonomy tables
otu_case <- as.data.frame(otu_table(physeq_case))
otu_control <- as.data.frame(otu_table(physeq_control))
taxonomy_case <- as.data.frame(tax_table(physeq_case))
taxonomy_control <- as.data.frame(tax_table(physeq_control))

# Identify unique phyla in both datasets
unique_phyla <- unique(c(taxonomy_case$Family, taxonomy_control$Family))
unique_phyla <- unique_phyla[!is.na(unique_phyla)]  # Remove any NAs

# Loop through each Family and generate a plot
for (Family in unique_phyla) {
  # Get indices of OTUs belonging to the current Family
  case_idx <- which(taxonomy_case$Family == Family)
  control_idx <- which(taxonomy_control$Family == Family)
  
  # Extract OTU abundance for the Family from both groups
  Family_case <- otu_case[case_idx, , drop = FALSE]
  Family_control <- otu_control[control_idx, , drop = FALSE]
  
  # Sum the abundance for the Family in each sample
  abundance_case <- rowSums(Family_case)
  abundance_control <- rowSums(Family_control)
  
  # Create a data frame for plotting
  Family_data <- data.frame(
    Abundance = c(abundance_case, abundance_control),
    Group = rep(c("PFAS-Exposed", "Reference"), c(length(abundance_case), length(abundance_control)))
  )
  
  # Skip if no data for this Family
  if (nrow(Family_data) == 0) next
  
  # Perform statistical test (Wilcoxon test)
  wilcox_test <- wilcox.test(Abundance ~ Group, data = Family_data)
  
  # Generate the plot
  p <- ggplot(Family_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste("Relative Abundance of", Family), 
         x = "Group", y = "Abundance") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10, face = "bold")
    ) +
    annotate("text", x = 1.5, y = max(Family_data$Abundance, na.rm = TRUE), 
             label = paste("p-value =", round(wilcox_test$p.value, 4)), 
             size = 5, hjust = 0.5)
  
  # Save the plot
  filename <- paste0(output_dir, "/", gsub(" ", "_", Family), "_abundance_comparison.png")
  ggsave(filename = filename, plot = p, width = 12, height = 8, dpi = 300)
  
  # Print message for tracking progress
  message("Saved plot for: ", Family)
}
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

sample_metadata_path <- file.path("data", "metadata", "Sample_Metadata_New.tsv")

# Load and Preprocess Data
sample_metadata <- read.delim(sample_metadata_path, check.names = FALSE, sep = "\t", row.names = 1)
names(sample_metadata) <- make.names(names(sample_metadata), unique = TRUE)

# Filter PFAS_Exposed group
pfas_exposed_data <- sample_metadata %>% 
  filter(Condition == "PFAS_Exposed")

# Select degradation rate columns
pfas_exposed_data_k <- pfas_exposed_data %>%
  as_tibble() %>%
  dplyr::select(starts_with("K_"))

# Reshape to long format
pfas_exposed_data_long <- pfas_exposed_data_k %>%
  pivot_longer(cols = everything(), 
               names_to = "PFAS_Type", 
               values_to = "Elimination_Rate_Value")

# Define dark purple-based palette
dark_palette <- c(
  "#3E1F47", "#5D3A6B", "#7A4D89", "#9A67A1",
  "#BC8CC3", "#6E2C7B", "#421C52", "#8E5BA3"
)

# Plot
plot <- ggplot(pfas_exposed_data_long, aes(x = PFAS_Type, y = Elimination_Rate_Value, fill = PFAS_Type)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black", width = 0.6) +
  scale_fill_manual(values = dark_palette) +
  theme_minimal(base_size = 14) +
  labs(
    x = "PFAS Type", 
    y = "Elimination Rate Value", 
    title = ""
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "#2D0033"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12, color = "#2D0033"),
    axis.text.y = element_text(face = "bold", size = 12, color = "#2D0033"),
    axis.title = element_text(face = "bold", size = 14, color = "#2D0033"),
    legend.position = "none", # remove redundant legend
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10) # tighter margins, less whitespace
  )

# Save
ggsave(
  file.path("results", "supplementary_figures", "Supp2_pfas_elimination_boxplot.png"), 
  plot = plot, 
  width = 9, height = 5, dpi = 300, bg = 'white'
)


#########################################
#########################################
