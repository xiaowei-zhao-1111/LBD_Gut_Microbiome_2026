library(pheatmap)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

###############################################################################
folder_path <- "~/Desktop/LBD/analysis/differential_prevalence_analysis/"
setwd(folder_path)
###############################################################################

# 1. Read the species relative abundance data
species_data <- read.table(file = "~/Desktop/LBD/analysis/metaphlan_results_new.tsv", sep = "\t", header = T, row.names = 1)
species <- species_data[grepl("s__", rownames(species_data)) & !grepl("t__", rownames(species_data)), ]

taxonomy <- read.csv(file = "~/Desktop/LBD/analysis/taxonomy_info.csv")

# 2. Modify taxonomy names
new_rownames <- gsub(".*s__", "", rownames(species))
rownames(species) <- new_rownames

# 3. Modify sample names
col_names <- names(species)
new_colnames <- gsub("metaphlan_|_S.*", "", col_names)
colnames(species) <- new_colnames

# 4. Read the metadata
metadata <- read.csv(file = "/Users/M306307/Desktop/LBD/analysis/imputed_BMI_metadata.csv", sep = ",", header = T, row.names = 1)
species_clean <- species[, c(rownames(metadata))]

# 5. Change the taxonomic data from percentage into proportion
colSums(species_clean)
species_prop <- data.frame(apply(species_clean, 2, function(x) x/sum(x)))
colSums(species_prop)

# 6. Use prevalence cut-off (10^-4.5) to remove bacterial species that are present in very low abundance
## Number of low abundance species in each sample
cut_off <- 10^-4.7
species_cutoff <- species_prop
species_cutoff[species_cutoff <= cut_off] <- 0

# 7. Change the relative abundance data to prevalence data
prevalence_data <- as.data.frame(ifelse(species_cutoff > 0, 1, 0))

# 7. Add metadata to the prevalence data
prevalence_t <- as.data.frame(t(prevalence_data))
prevalence_all <- merge(prevalence_t, metadata, by = "row.names", all = F)
rownames(prevalence_all) <- prevalence_all$Row.names
prevalence_all <- prevalence_all[, -1]

###############################################################################

# Find the bacterial species that are passing the prevalence cut-off (10%)

prevalence_cutoff <- function(prev_data, metadata) {
  
  filtered_data <- prev_data[, c(1:(ncol(prev_data)-33))][, colSums(prev_data[, c(1:(ncol(prev_data)-33))] > 0) >= 0.1*nrow(prev_data)]
  result_df <- merge(filtered_data, metadata, by = "row.names", all = F)
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

heatmap_z <- function(sig, 
                      label1, 
                      label2, 
                      sig_file_name, 
                      heatmap_file_name, 
                      taxonomy_table = taxonomy) {
  
  # Merge with taxonomy by species name
  sig_done <- merge(sig, taxonomy_table, by = "species", all = FALSE)
  
  # Store selected species (optional)
  selected_species <- sig_done$species
  
  # Create 'condition' column based on comparison of the 3rd and 4th columns
  sig_done$condition <- ifelse(sig_done[, 3] > sig_done[, 4], label1, label2)
  
  # Split and sort by descending abundance
  sig_1 <- sig_done[sig_done$condition == label1, ]
  sig_1 <- sig_1[order(-sig_1[, 3]), ]
  
  sig_2 <- sig_done[sig_done$condition == label2, ]
  sig_2 <- sig_2[order(-sig_2[, 4]), ]
  
  # Combine and set row names
  sig_done_1 <- rbind(sig_1, sig_2)
  row.names(sig_done_1) <- sig_done_1$species
  
  # Save full table
  write.csv(sig_done_1, file = sig_file_name, row.names = FALSE)
  
  # Plot heatmap for 3rd and 4th columns
  pdf(file = heatmap_file_name)
  pheatmap(data.matrix(sig_done_1[3:4]), 
           color = colorRampPalette(c("#f7f3e8", "#b1182d"))(100), 
           border_color = "white",
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           cellwidth = 25, 
           cellheight = 10,
           scale = "none", 
           main = "",
           angle_col = 45)
  dev.off()
}


#. Create barplot for each significant species

barplot_z <- function(sig, level1, level2, sig_names, color1, color2, name) {
  
  barplot_data <- melt(sig[, c(1, 3:4)], id.vars = "species")
  barplot_data$variable <- factor(barplot_data$variable, levels = c(level1, level2))
  
  for (i in 1:nrow(sig)) {
    plot <- ggplot(barplot_data[barplot_data$species == sig_names[i], ], 
                   aes(x = variable, y = value*100, fill = variable)) + 
      geom_bar(stat = "identity", color = "black") +
      scale_fill_manual(values = c(color1, color2)) +
      theme_bw() + 
      theme(
        panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
      labs(x = "", 
           y = "Prevalence (%)",
           title = sig_names[i])
    
    pdf(file = paste0("~/Desktop/LBD/analysis/differential_prevalence_analysis/", name, "_", lbd_vs_control_sig_names[i], ".pdf"), height = 6, width = 4.5)
    print(plot)
    dev.off()
  }
  
}

barplot_z <- function(sig, level1, level2, sig_names, color1, color2, name) {
  
  # Reshape input to long format for ggplot
  barplot_data <- melt(sig[, c(1, 3:4)], id.vars = "species")
  
  # Set factor levels to control bar order
  barplot_data$variable <- factor(barplot_data$variable, levels = c(level1, level2))
  
  # Loop over significant species to create individual barplots
  for (i in 1:nrow(sig)) {
    
    # Subset data for the current species
    single_species <- barplot_data[barplot_data$species == sig_names[i], ]
    
    # Create barplot
    plot <- ggplot(single_species, 
                   aes(x = variable, y = value * 100, fill = variable)) + 
      geom_bar(stat = "identity", color = "black") +
      scale_fill_manual(values = c(color1, color2)) +
      theme_bw() + 
      theme(
        panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
      labs(
        x = "", 
        y = "Prevalence (%)",
        title = sig_names[i]
      )
    
    # Save plot as PDF
    pdf(file = paste0("~/Desktop/LBD/analysis/differential_prevalence_analysis/", name, "_", sig_names[i], ".pdf"), 
        height = 6, width = 4.5)
    print(plot)
    dev.off()
  }
}


###############################################################################

# LBD vs. Control#

# 1. Extract data for LBD and their controls
lbd_vs_control_prev <- prevalence_all[prevalence_all$condition %in% c("lbd", "lbd_control"), ]

# 2. Find the bacterial prevalence that are passing the cut-off
lbd_vs_control_prev_filtered <- prevalence_cutoff(lbd_vs_control_prev, metadata)

# 3. Fisher's exact test
lbd_vs_control_prev_results <- fisher_exact_test(lbd_vs_control_prev_filtered, "lbd_vs_control_prevalence_results.csv", "LBD", "LBD_control", condition_levels = c("lbd", "lbd_control"))

write.csv(lbd_vs_control_prev_results, file = "lbd_vs_control_prev_results.csv")


# 4. Significant prevalent species
sig_lbd_vs_control_prev <- lbd_vs_control_prev_results[lbd_vs_control_prev_results$p_value < 0.05, ]
sig_lbd_vs_control_prev_ordered <- sig_lbd_vs_control_prev[order(sig_lbd_vs_control_prev$p_value), ]
sig_lbd_vs_control_prev_ordered_done <- sig_lbd_vs_control_prev_ordered[rowSums(is.na(sig_lbd_vs_control_prev_ordered)) != ncol(sig_lbd_vs_control_prev_ordered), ]
write.csv(sig_lbd_vs_control_prev_ordered_done, file = "sig_lbd_vs_control_prev_new.csv", row.names = F)
lbd_vs_control_sig_names <- sig_lbd_vs_control_prev_ordered$species

# 5. Heatmap
heatmap_z(sig_lbd_vs_control_prev_ordered, "LBD>Control", "LBD<Control", "sig_results_lbd_vs_control.csv", "lbd_vs_control_diff_prev_V2.pdf")

# 6. Create barplot for each significant species
barplot_z(sig_lbd_vs_control_prev_ordered, "LBD", "LBD_control", lbd_vs_control_sig_names, "#ff9274", "#55b7e6", "lbd_vs_control_boxplot")

 ###############################################################################

# iRBD vs. Control#

# 1. Extract data for iRBD and their controls
irbd_vs_control_prev <- prevalence_all[prevalence_all$condition %in% c("irbd", "irbd_control"), ]

# 2. Find the bacterial prevalence that are passing the cut-off
irbd_vs_control_prev_filtered <- prevalence_cutoff(irbd_vs_control_prev, metadata)

# 3. Fisher's exact test
irbd_vs_control_prev_results <- fisher_exact_test(irbd_vs_control_prev_filtered, "irbd_vs_control_prevalence_results.csv", "iRBD", "iRBD_control", condition_levels = c("irbd", "irbd_control"))

write.csv(irbd_vs_control_prev_results, file = "irbd_vs_control_prev_results.csv")

# 4. Significant prevalent species
sig_irbd_vs_control_prev <- irbd_vs_control_prev_results[irbd_vs_control_prev_results$p_value < 0.05, ]
sig_irbd_vs_control_prev_ordered <- sig_irbd_vs_control_prev[order(sig_irbd_vs_control_prev$p_value), ]
# write.csv(sig_irbd_vs_control_prev_ordered, file = "sig_irbd_vs_control_prev.csv", row.names = F)
irbd_vs_control_sig_names <- sig_irbd_vs_control_prev$species

# 5. Heatmap
heatmap_z(sig_irbd_vs_control_prev_ordered, "iRBD>Control", "iRBD<Control", "sig_results_irbd_vs_control.csv", "irbd_vs_control_diff_prev.pdf")

# 6. Create barplot for each significant species
barplot_z(sig_irbd_vs_control_prev_ordered, "iRBD", "iRBD_control", irbd_vs_control_sig_names, "#2DA248", "#fdc848", "irbd_vs_control_boxplot")

###############################################################################

# LBD vs. iRBD #

# 1. Extract data for LBD and iRBD
lbd_vs_irbd_prev <- prevalence_all[prevalence_all$condition %in% c("lbd", "irbd"), ]

# 2. Find the bacterial prevalence that are passing the cut-off
lbd_vs_irbd_prev_filtered <- prevalence_cutoff(lbd_vs_irbd_prev, metadata)

# 3. Fisher's exact test
lbd_vs_irbd_prev_results <- fisher_exact_test(lbd_vs_irbd_prev_filtered, "lbd_vs_irbd_prevalence_results.csv", "LBD", "iRBD")

# 4. Significant prevalent species
sig_lbd_vs_irbd_prev <- lbd_vs_irbd_prev_results[lbd_vs_irbd_prev_results$p_value < 0.05 & !is.na(lbd_vs_irbd_prev_results$p_value), ]
sig_lbd_vs_irbd_prev_ordered <- sig_lbd_vs_irbd_prev[order(sig_lbd_vs_irbd_prev$p_value), ]
write.csv(sig_lbd_vs_irbd_prev_ordered, file = "sig_lbd_vs_irbd_prev.csv", row.names = F)

# 5. Heatmap
heatmap_z(sig_lbd_vs_irbd_prev_ordered, "LBD>iRBD", "iRBD<LBD", "sig_results_lbd_vs_irbd.csv", "lbd_vs_irbd_diff_prev.pdf")

# 6. Create barplot for each significant species
barplot_z(sig_irbd_vs_control_prev_ordered, "iRBD", "iRBD_control", irbd_vs_control_sig_names, "#2DA248", "#fdc848", "irbd_vs_control_boxplot")

###############################################################################

## species are more prevalent in both LBD and iRBD
intersect(sig_lbd_vs_control_prev$species, sig_irbd_vs_control_prev$species)

###############################################################################

sig_species_all <- c(lbd_vs_control_sig_names, irbd_vs_control_sig_names)

lbd_df <- lbd_vs_control_prev_filtered[, colnames(lbd_vs_control_prev_filtered) %in% c(sig_species_all, "condition")]

lbd_results <- lbd_df %>%
  pivot_longer(-condition, names_to = "species", values_to = "presence") %>%
  group_by(condition, species) %>%
  summarise(prevalence = mean(presence) * 100, .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = prevalence)

irbd_df <- irbd_vs_control_prev_filtered[, colnames(irbd_vs_control_prev_filtered) %in% c(sig_species_all, "condition")]
irbd_results <- irbd_df %>%
  pivot_longer(-condition, names_to = "species", values_to = "presence") %>%
  group_by(condition, species) %>%
  summarise(prevalence = mean(presence) * 100, .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = prevalence)

lbd_irbd_df <- merge(lbd_results, irbd_results, by = "species", all = T)

lbd_irbd_long_df <- melt(lbd_irbd_df, id.vars = "species")

lbd_irbd_long_done <- lbd_irbd_long_df %>%
  mutate(variable = case_when(
    variable == "lbd" ~ "LBD",
    variable == "irbd" ~ "iRBD",
    variable == "lbd_control" ~ "LBD_Control",
    variable == "irbd_control" ~ "iRBD_Control",
    TRUE ~ variable  # leave other values unchanged
  ))

lbd_irbd_long_done$variable <- factor(lbd_irbd_long_done$variable, levels = c("LBD", "LBD_Control", "iRBD", "iRBD_Control"))

for (i in 1:nrow(lbd_irbd_df)) {

  name <- lbd_irbd_df$species[i]
  pdf(file = paste0("barplot_", name, ".pdf"), width = 4, height = 5.5)
  
  plot <- ggplot(filter(lbd_irbd_long_done, species == name), 
         aes(x = variable, y = value, fill = variable)) +
    geom_bar(stat = "identity", color = "black") +
    labs(x = "", y = "Prevalence (%)", title = name) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.grid.major.y = element_line(linetype = "dashed", color = "gray70", linewidth = 0.1),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = c("LBD" = "#ff9274", 
                                 "LBD_Control" = "#55b7e6",
                                 "iRBD" = "#2DA248",
                                 "iRBD_Control" = "#fdc848")) + 
    rremove("legend")
  print(plot)
  dev.off()
  
}

# Define species of interest
group1 <- c("Oscillospiraceae_bacterium_CLA_AA_H250", "Akkermansia_muciniphila", "Actinomyces_oris", "GGB3730_SGB5060", "GGB9627_SGB15081")

# Filter data for selected species
lbd_irbd_long_done_g1 <- lbd_irbd_long_done[lbd_irbd_long_done$species %in% group1, ]

# Get species order based on decreasing value in LBD group
order1 <- lbd_irbd_long_done_g1 %>%
  filter(variable == "LBD") %>%
  arrange(desc(value)) %>%
  pull(species)

# Apply factor levels to reorder species
lbd_irbd_long_done_g1$species <- factor(lbd_irbd_long_done_g1$species, levels = order1)

pdf(file = "barplot_group1.pdf", width = 10, height = 6)
plot <- ggplot(lbd_irbd_long_done_g1, 
               aes(x = species, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  labs(x = "", y = "Prevalence (%)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.y = element_line(linetype = "dashed", color = "gray70", linewidth = 0.1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("LBD" = "#ff9274", 
                               "LBD_Control" = "#55b7e6",
                               "iRBD" = "#2DA248",
                               "iRBD_Control" = "#fdc848"))
print(plot)
dev.off()

group2 <- c("Clostridia_bacterium_UC5_1_1D1", "Longicatena_caecimuris", "GGB9719_SGB15272")

# Filter data for selected species
lbd_irbd_long_done_g2 <- lbd_irbd_long_done[lbd_irbd_long_done$species %in% group2, ]

# Get species order based on decreasing value in LBD group
order1 <- lbd_irbd_long_done_g2 %>%
  filter(variable == "LBD") %>%
  arrange(desc(value)) %>%
  pull(species)

# Apply factor levels to reorder species
lbd_irbd_long_done_g2$species <- factor(lbd_irbd_long_done_g2$species, levels = order1)

pdf(file = "barplot_group2.pdf", width = 8, height = 6)
plot <- ggplot(lbd_irbd_long_done_g2, 
               aes(x = species, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  labs(x = "", y = "Prevalence (%)") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.y = element_line(linetype = "dashed", color = "gray70", linewidth = 0.1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("LBD" = "#ff9274", 
                               "LBD_Control" = "#55b7e6",
                               "iRBD" = "#2DA248",
                               "iRBD_Control" = "#fdc848"))
print(plot)
dev.off()


group <- c(group1, group2)
lbd_df <- lbd_vs_control_prev_results[lbd_vs_control_prev_results$species %in% group, ]
colnames(lbd_df)[2] <- "lbd_p_value"
irbd_df <- irbd_vs_control_prev_results[irbd_vs_control_prev_results$species %in% group, ]
colnames(irbd_df)[2] <- "irbd_p_value"

lbd_irbd_df <- merge(lbd_df, irbd_df, by = "species", all = T)
lbd_irbd_df_done <- lbd_irbd_df[, c(1, 3, 4, 6, 7)]

lbd_irbd_df_done$species <- factor(lbd_irbd_df_done$species, levels = c(group))
lbd_irbd_df_done_ordered <- lbd_irbd_df_done %>% arrange(species)
write.csv(lbd_irbd_df_done_ordered, file = "lbd_and_irbd_prev_species.csv", row.names = F)

pdf(file = "heatmap_all_lbd_and_irbd.pdf", width = 5.5, height = 3)
pheatmap(data.matrix(lbd_irbd_df_done_ordered[2:5]), 
         color = colorRampPalette(c("#f7f3e8", "#b1182d"))(100), 
         border_color = "white",
         cluster_rows = F, 
         cluster_cols = F, 
         cellwidth = 30, cellheight = 25,
         scale = "none", 
         main = "",
         angle_col = "45",
         labels_row = lbd_irbd_df_done_ordered$species,
         breaks = seq(0, 1, length.out = 101),
         legend_breaks = seq(0, 1, 0.2),                   # ticks at 0,0.2,…,1
         legend_labels = seq(0, 100, 20))
dev.off()


library(ComplexHeatmap)
library(circlize)

mat <- data.matrix(lbd_irbd_df_done_ordered[2:5])

Heatmap(mat,
        name = "Abundance",
        col = colorRamp2(c(0, 1), c("#f7f3e8", "#b1182d")),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_labels = lbd_irbd_df_done_ordered$species,
        heatmap_legend_param = list(
          at = seq(0, 1, 0.2),          # tick mark positions
          labels = seq(0, 1, 0.2),      # label text
          title = "Relative abundance"
        ))

