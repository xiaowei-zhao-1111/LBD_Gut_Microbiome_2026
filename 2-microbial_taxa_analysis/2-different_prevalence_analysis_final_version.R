library(pheatmap)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

####################################################################
setwd("../2-microbial_taxa_analysis")
####################################################################

# 1. Read the preprocessed species relative abundance data
species_all <- read.csv(file = "species_preprocessed.csv", row.names = 1, header = T)

# 2. Read the taxonomy information
taxonomy <- read.csv(file = "../0-raw_data/taxonomy_info.csv")

# 3. Read the metadata
metadata <- read.csv(file = "../0-raw_data/imputed_BMI_metadata.csv", sep = ",", header = T, row.names = 1)

# 4. Change the relative abundance data to prevalence data
species_only <- species_all[, c(1:(ncol(species_all)-33))]
dim(species_only)
species_prevalence <- as.data.frame(ifelse(species_only > 0, 1, 0))
dim(species_prevalence)

# 7. Add metadata to the prevalence data
prevalence_all <- merge(species_prevalence, metadata, by = "row.names", all = T)
rownames(prevalence_all) <- prevalence_all$Row.names
prevalence_all <- prevalence_all[, -1]
View(prevalence_all)

###############################################################################

# Find the bacterial species that are passing the prevalence cut-off (10%)

prevalence_cutoff <- function(prev_data, metadata) {

  # Filter species columns based on prevalence cut-off (10%)
  filtered_data <- prev_data[, c(1:(ncol(prev_data)-33))][, colSums(prev_data[, c(1:(ncol(prev_data)-33))] > 0) >= 0.1*nrow(prev_data)]
  
  # Add metadata back to the filtered data
  result_df <- merge(filtered_data, metadata, by = "row.names", all = F)
  
  # Set row names and remove the first column (Row.names)
  rownames(result_df) <- result_df$Row.names
  result_df <- result_df[, -1]
  
  return(result_df)
  
}

fisher_exact_test <- function(filtered_data, file_name, name1, name2, condition_levels = c("lbd", "lbd_control")) {
  
  # Extract species columns and condition column
  # Assumes last 33 columns are metadata; keep species + condition
  fisher_data <- filtered_data[, c(1:(ncol(filtered_data) - 33), ncol(filtered_data) - 30)]
  
  # Ensure condition column is a factor with specified levels
  fisher_data$condition <- factor(fisher_data$condition, levels = condition_levels)
  
  # Initialize output vectors
  p_values <- character()
  species_names <- character()
  group1_perc <- numeric()
  group2_perc <- numeric()
  
  # Iterate over each species column
  for (i in 1:(ncol(fisher_data) - 1)) {
    
    species_names[i] <- colnames(fisher_data)[i]
    
    # Build contingency table: presence/absence vs condition
    fisher_table <- as.matrix(table(fisher_data[, i], fisher_data$condition))
    
    # Normalize table to get column-wise percentages
    perc <- sweep(fisher_table, 2, colSums(fisher_table), FUN = "/")
    
    # Handle case when only one unique value (e.g., all 0s or 1s)
    if (nrow(fisher_table) < 2) {
      p_values[i] <- NA
      group1_perc[i] <- perc[1, 1]
      group2_perc[i] <- perc[1, 2]
    } else {
      p_values[i] <- fisher.test(fisher_table)[["p.value"]]
      group1_perc[i] <- perc[2, 1]
      group2_perc[i] <- perc[2, 2]
    }
  }
  
  # Combine results into a dataframe
  results <- data.frame(
    species = species_names,
    p_value = p_values,
    group1_perc,
    group2_perc,
    stringsAsFactors = FALSE
  )
  
  results$qval <- p.adjust(results$p_value, method = "BH")
  # Rename group columns using user-defined names
  colnames(results)[3:4] <- c(name1, name2)
  
  # Save results to CSV
  write.csv(results, file = file_name, row.names = FALSE)
  
  # Return results
  return(results)
}

###############################################################################

# LBD vs. Control

###############################################################################

# 1. Create the folder if not exist
folder_name <- "3-lbd_vs_control_diff_prevalence_results"
folder <- file.path(".", folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

# 2. Extract data for LBD and controls
lbd_vs_control_prev <- prevalence_all[prevalence_all$con %in% c("lbd", "lbd_control"), ]
dim(lbd_vs_control_prev)

# 3. Select microbial taxa that are passing the prevalence cut-off (10%)
lbd_vs_control_prev_filtered <- prevalence_cutoff(lbd_vs_control_prev, metadata)
dim(lbd_vs_control_prev_filtered)

# 4. Fisher's exact test
lbd_vs_control_prev_results <- fisher_exact_test(lbd_vs_control_prev_filtered, "./3-lbd_vs_control_diff_prevalence_results/lbd_vs_control_prevalence_results.csv", "LBD", "LBD_control", condition_levels = c("lbd", "lbd_control"))

# 5. Select significant prevalent species
sig_lbd_vs_control_prev <- lbd_vs_control_prev_results[lbd_vs_control_prev_results$p_value < 0.05, ]

# 5.1 Order the significant prevalent species by p-value and remove rows with all NA values
sig_lbd_vs_control_prev_ordered <- sig_lbd_vs_control_prev[order(sig_lbd_vs_control_prev$p_value), ]
sig_lbd_vs_control_prev_ordered_done <- sig_lbd_vs_control_prev_ordered[rowSums(is.na(sig_lbd_vs_control_prev_ordered)) != ncol(sig_lbd_vs_control_prev_ordered), ]

# 5.2 Save the significant prevalent species to a CSV file
write.csv(sig_lbd_vs_control_prev_ordered_done, file = "./3-lbd_vs_control_diff_prevalence_results/sig_lbd_vs_control_prev.csv", row.names = F)

###############################################################################

# iRBD vs. Control

###############################################################################

# 1. Create the folder if not exist
folder_name <- "4-irbd_vs_control_diff_prevalence_results"
folder <- file.path(".", folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

# 2. Extract data for iRBD and their controls
irbd_vs_control_prev <- prevalence_all[prevalence_all$condition %in% c("irbd", "irbd_control"), ]
dim(irbd_vs_control_prev)

# 3. Select microbial taxa that are passing the prevalence cut-off (10%)
irbd_vs_control_prev_filtered <- prevalence_cutoff(irbd_vs_control_prev, metadata)
dim(irbd_vs_control_prev_filtered)

# 4. Fisher's exact test
irbd_vs_control_prev_results <- fisher_exact_test(irbd_vs_control_prev_filtered, "./4-irbd_vs_control_diff_prevalence_results/irbd_vs_control_prevalence_results.csv", "iRBD", "iRBD_control", condition_levels = c("irbd", "irbd_control"))

# 5. Select significant prevalent species
sig_irbd_vs_control_prev <- irbd_vs_control_prev_results[irbd_vs_control_prev_results$p_value < 0.05, ]

# 5.1 Order the significant prevalent species by p-value and remove rows with all NA values
sig_irbd_vs_control_prev_ordered <- sig_irbd_vs_control_prev[order(sig_irbd_vs_control_prev$p_value), ]
sig_irbd_vs_control_prev_ordered_done <- sig_irbd_vs_control_prev_ordered[rowSums(is.na(sig_irbd_vs_control_prev_ordered)) != ncol(sig_irbd_vs_control_prev_ordered), ]

# 5.2 Save the significant prevalent species to a CSV file
write.csv(sig_irbd_vs_control_prev_ordered_done, file = "./4-irbd_vs_control_diff_prevalence_results/sig_irbd_vs_control_prev.csv", row.names = F)

# 6. Species are more prevalent in both LBD and iRBD
intersect(sig_lbd_vs_control_prev_ordered_done$species, sig_irbd_vs_control_prev_ordered_done$species)

###############################################################################

# Combine LBD and iRBD results together

###############################################################################

# 1. Combine the significant species from both LBD and iRBD analyses
sig_species_all <- c(sig_lbd_vs_control_prev_ordered_done$species, sig_irbd_vs_control_prev_ordered_done$species)

# 2. Filter the LBD vs Control prevalence data to include only the significant species and the condition column
lbd_df <- lbd_vs_control_prev_filtered[, colnames(lbd_vs_control_prev_filtered) %in% c(sig_species_all, "condition")]

# 3. Calculate prevalence for each species in LBD and Control groups
lbd_results <- lbd_df %>%
  pivot_longer(-condition, names_to = "species", values_to = "presence") %>%
  group_by(condition, species) %>%
  summarise(prevalence = mean(presence) * 100, .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = prevalence)

# 4. Filter the iRBD vs Control prevalence data to include only the significant species and the condition column
irbd_df <- irbd_vs_control_prev_filtered[, colnames(irbd_vs_control_prev_filtered) %in% c(sig_species_all, "condition")]

# 5. Calculate prevalence for each species in iRBD and Control groups
irbd_results <- irbd_df %>%
  pivot_longer(-condition, names_to = "species", values_to = "presence") %>%
  group_by(condition, species) %>%
  summarise(prevalence = mean(presence) * 100, .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = prevalence)

# 6. Merge the LBD and iRBD results into a single dataframe
lbd_irbd_df <- merge(lbd_results, irbd_results, by = "species", all = T)

# 7. Rename the columns for clarity
lbd_irbd_long_df <- melt(lbd_irbd_df, id.vars = "species")

# 8. Rename the levels of the 'variable' column for clarity
lbd_irbd_long_done <- lbd_irbd_long_df %>%
  mutate(variable = case_when(
    variable == "lbd" ~ "LBD",
    variable == "irbd" ~ "iRBD",
    variable == "lbd_control" ~ "LBD_Control",
    variable == "irbd_control" ~ "iRBD_Control",
    TRUE ~ variable  # leave other values unchanged
  ))

# 9. Set the order of the 'variable' factor for consistent plotting
lbd_irbd_long_done$variable <- factor(lbd_irbd_long_done$variable, levels = c("LBD", "LBD_Control", "iRBD", "iRBD_Control"))

###############################################################################

# Species gradually increase or decrease from iRBD to LBD

###############################################################################

# 1. Define species increase gradually from iRBD to LBD
group1 <- c("Oscillospiraceae_bacterium_CLA_AA_H250", "Akkermansia_muciniphila", "Actinomyces_oris", "GGB3730_SGB5060", "GGB9627_SGB15081")

# 2. Define species decrease gradually from iRBD to LBD
group2 <- c("Clostridia_bacterium_UC5_1_1D1", "Longicatena_caecimuris", "GGB9719_SGB15272")

# 3. Combine the two groups into a single vector
group <- c(group1, group2)
lbd_df <- lbd_vs_control_prev_results[lbd_vs_control_prev_results$species %in% group, ]
colnames(lbd_df)[2] <- "lbd_p_value"
irbd_df <- irbd_vs_control_prev_results[irbd_vs_control_prev_results$species %in% group, ]
colnames(irbd_df)[2] <- "irbd_p_value"

# 4. Merge the LBD and iRBD dataframes based on species
lbd_irbd_df <- merge(lbd_df, irbd_df, by = "species", all = T)
lbd_irbd_df_done <- lbd_irbd_df[, c(1, 3, 4, 7, 8)]

# 5. Set the order of the 'species' factor for consistent plotting
lbd_irbd_df_done$species <- factor(lbd_irbd_df_done$species, levels = c(group))
lbd_irbd_df_done_ordered <- lbd_irbd_df_done %>% arrange(species)

# 6. Save the combined LBD and iRBD prevalence data to a CSV file
write.csv(lbd_irbd_df_done_ordered, file = "./5-others/lbd_and_irbd_prev_species.csv", row.names = F)

# 7. Create a heatmap of the combined LBD and iRBD prevalence data
pdf(file = "./5-others/heatmap_all_lbd_and_irbd.pdf", width = 5.5, height = 3)
pheatmap(data.matrix(lbd_irbd_df_done_ordered[2:5]), 
         color = colorRampPalette(c("#f7f3e8", "#b1182d"))(100), 
         border_color = "white",
         cluster_rows = F, 
         cluster_cols = F, 
         cellwidth = 30, cellheight = 15,
         scale = "none", 
         main = "",
         angle_col = "45",
         labels_row = lbd_irbd_df_done_ordered$species,
         breaks = seq(0, 1, length.out = 101),
         legend_breaks = seq(0, 1, 0.2),
         legend_labels = seq(0, 100, 20))
dev.off()