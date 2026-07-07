library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)
library(lmerTest)
library(tidyr)
library(dplyr)
library(effectsize)

####################################################################
setwd("../2-microbial_taxa_analysis")
####################################################################

# 1. Read the species relative abundance data
species_data <- read.table(file = "../0-raw_data/metaphlan_results_new.tsv", sep = "\t", header = T, row.names = 1)
species <- species_data[grepl("s__", rownames(species_data)) & !grepl("t__", rownames(species_data)), ]

# 2. Read the taxonomy information
taxonomy <- read.csv(file = "../0-raw_data/taxonomy_info.csv")

# 3. Modify taxonomy names
new_rownames <- gsub(".*s__", "", rownames(species))
rownames(species) <- new_rownames

# 4. Modify sample names
col_names <- names(species)
new_colnames <- gsub("metaphlan_|_S.*", "", col_names)
colnames(species) <- new_colnames

# 5. Read the metadata
metadata <- read.csv(file = "../0-raw_data/imputed_BMI_metadata.csv", sep = ",", header = T, row.names = 1)

# 6. Clean the species data to only include samples present in the metadata
species_clean <- species[, c(rownames(metadata))]

# 7. Change the taxonomic data from percentage into proportion
colSums(species_clean)
species_prop <- data.frame(apply(species_clean, 2, function(x) x/sum(x)))
colSums(species_prop)

# 8. Find the appropriate presence cutoff for bacterial species data
relab_number <- unlist(species_prop, use.name = FALSE)
relab_number_ordered <- relab_number[order(relab_number, decreasing = T)]
## Create a rank plot to visualize the ordered relative abundance of microbial species and determine a suitable cutoff for low abundance species.
pdf("./rank_plot_with_threshold_microbial_species.pdf", height = 6, width = 8)
plot(x=c(1:length(relab_number_ordered)), y=log10(relab_number_ordered), pch = 20, xlab = "Rank", ylab = "Relative abundance of microbial species", main = "Rank plot of ordered relative abundance of microbial species", xlim = c(0, 14000))
abline(h=-4.7, col = "orange", lwd = 3)
dev.off()

## It seems 10^-4.7 is the good cutoff

# 9. Apply the cutoff to the species proportion data, setting values below the threshold to zero. This step helps to filter out low-abundance species that may not be biologically relevant or may introduce noise into downstream analyses.
cut_off <- 10^-4.7
species_cutoff <- species_prop
species_cutoff[species_cutoff < cut_off] <- 0

# 10. Add metadata to the proportion data
species_t <- as.data.frame(t(species_cutoff))
species_all <- merge(species_t, metadata, by = "row.names", all = F)
rownames(species_all) <- species_all$Row.names
species_all <- species_all[, -1]
write.csv(species_all, "species_preprocessed.csv")

####################################################################

## The following function performs a transformation on a dataset and merges it with metadata and group information to produce a final structured data frame.
trs <- function(dataset) {
  transformed <- asin(sqrt(dataset[, 1:(ncol(dataset)-33)]))
  return(transformed)
} 

# Find the bacterial species that are passing the prevalence cut-off (10%)
prevalence_cutoff <- function(transformed_data, metadata) {

  # Filter the transformed data to include only species that are present in more than 10% of the samples. 
  filtered_data <- transformed_data[, colSums(transformed_data > 0) >= 0.1*nrow(transformed_data)]

  # Merge the filtered data with the metadata to create a final structured data frame that includes both the transformed abundance values and relevant sample information.
  result_df <- merge(filtered_data, metadata, by = "row.names", all = F)

  # Set the row names of the resulting data frame to the sample identifiers and remove the redundant 'Row.names' column to clean up the final output.
  rownames(result_df) <- result_df$Row.names

  # Remove the 'Row.names' column from the data frame to avoid duplication and maintain a clean structure for further analysis.
  result_df <- result_df[, -1]

  return(result_df) 

}

# The following function performs a mixed-effects model analysis to compare the abundance of species between two time points. It calculates statistical significance (p-value) and trend direction (difference) for each species across the given time points.
mixed_effect <- function(filtered_data, file_name) {
  # Use a mixed-effects linear model to analyze the relationship between species abundance and the condition variable, while accounting for potential confounding factors such as age and BMI. The model includes a random effect for household_id to account for within-household correlations.
  p_diff_relab <- c()
  for (i in 1:(ncol(filtered_data)-33)) {
    p_diff_relab[i] <- tryCatch({summary(lmer(filtered_data[, i] ~ condition + age + BMI + (1|household_id), data = filtered_data, REML = F))[["coefficients"]][2,5]}, error = function(e) NA) # If error occurs, assign NA and continue
  }

  # Create a data frame to store the results
  results <- data.frame("species" = names(filtered_data)[1:(ncol(filtered_data)-33)], "p_value" = p_diff_relab)

  # Adjust the p-values for multiple testing using the Benjamini-Hochberg (BH) method to control the false discovery rate (FDR).
  results$q_value <- p.adjust(results$p_value, method = "BH")

  # Merge the results with the taxonomy information to provide a comprehensive view of the species analyzed, including their taxonomic classification.
  results_df <- merge(results, taxonomy, by = "species", all = F)

  # Write the results to a CSV file for further analysis or reporting.
  write.csv(results_df, file = file_name, row.names = F)

  return(results_df)

}

# The following function calculates the log2 fold-change (log2fc) of species abundance between two groups, integrates p-values from a mixed-effects linear model, and saves the results to a CSV file. It also returns a list containing the cleaned log2fc results and the merged data frame with p-values.
log2_fc_results <- function(filtered_data, group1, group2, p_val_results, file_name) {

  # Re-name the column names so it contains the group labels
  fc_data <- filtered_data[, c(1:(ncol(filtered_data)-33), ncol(filtered_data)-30)]

  # Create a new identifier for each row by combining the condition and the original row names. 

  fc_data$new_id <- paste(fc_data$condition, rownames(fc_data), sep = "_")

  # Set the row names of the data frame to the new identifiers and remove the redundant 'new_id' column to clean up the final output.
  rownames(fc_data) <- fc_data$new_id

  # Remove the 'condition' column from the data frame to avoid duplication and maintain a clean structure for further analysis.
  fc_data <- fc_data[, 1:(ncol(filtered_data)-33)]

  # Transpose the data frame to have species as rows and samples as columns.
  fc_data_t <- as.data.frame(t(fc_data))

  # Compute medians for each group
  fc_data_t$median_g1 <- apply(fc_data_t[, grepl(group1, colnames(fc_data_t))], 1, function(x) median(x))
  fc_data_t$median_g2 <- apply(fc_data_t[, grepl(group2, colnames(fc_data_t))], 1, function(x) median(x))

  # Compute log2 median fold change
  pseudocount <- cut_off
  fc_data_t$log2fc <- log2((fc_data_t$median_g1 + pseudocount)/(fc_data_t$median_g2 + pseudocount))

  # Ensure species names from p_val and log2fc results are the same
  if (!all(rownames(fc_data_t) %in% p_val_results$species)) {
    stop("Error: Some species in `fc_data_t` are missing in `p_results`.")
  }

  # Ensure there is a column named "p_value" in the p_val_results
  if (!"p_value" %in% colnames(p_val_results)) {
    stop("Error: `p_results` must contain a `p_value` column.")
  }

  # Merge p-values with log2fc results
  fc_data_t_p <- merge(fc_data_t, p_val_results, by.x = "row.names", by.y = "species", all = F)
  fc_data_t_p_clean <- fc_data_t_p[, c(1, (ncol(fc_data_t_p)-14):ncol(fc_data_t_p))]
  colnames(fc_data_t_p_clean)[1] <- "species"

  # Write the cleaned log2fc results with p-values to a CSV file for further analysis or reporting
  write.csv(fc_data_t_p_clean, file = file_name, row.names = F)

  return(list(fc_data_t_p_clean, fc_data_t_p))

}

## The following function processes the filtered dataset containing species abundance values for two conditions, computes the log2 fold-change, integrates p-values from mixed-effects linear model, and categorizes significant differences, and generates a volcano plot. 
volcano_plot <- function(fc_data_t_p, label1, label2, xlim1, xlim2) {

  # Categorize species based on log2 fold-change and p-value thresholds. Species with a log2 fold-change greater than 0 and a p-value less than 0.05 are categorized as "g1", while those with a log2 fold-change less than 0 and a p-value less than 0.05 are categorized as "g2". All other species are categorized as "All".
  g1 <- paste0(label1, ">", label2)
  g2 <- paste0(label1, "<", label2)
  fc_data_t_p$category <- "All"
  fc_data_t_p$category[fc_data_t_p$log2fc > 0 & fc_data_t_p$p_value < 0.05] <- g1
  fc_data_t_p$category[fc_data_t_p$log2fc < 0 & fc_data_t_p$p_value < 0.05] <- g2
  fc_data_t_p$category <- factor(fc_data_t_p$category, levels = c("All", g1, g2))
  
  # Extract significant data and make them in an increasing order based on p_value
  positive_fold <- fc_data_t_p %>% filter(category == g1) %>% arrange(p_value)
  negative_fold <- fc_data_t_p %>% filter(category == g2) %>% arrange(p_value)
  
  # Write the significant species with positive and negative fold changes to separate CSV files for further analysis or reporting
  write.csv(positive_fold, file = paste0("species_higher_in_", label1, "_than_", label2, ".csv"), row.names = F)
  write.csv(negative_fold, file = paste0("species_higher_in_", label2, "_than_", label1, ".csv"), row.names = F)
  
  # Create a volcano plot to visualize the relationship between log2 fold-change and p-value for the species analyzed
  plot <- ggplot(data = fc_data_t_p, aes(x = log2fc, y = -log10(p_value), color = category)) + 
    geom_point(size = 1, color = "lightgrey") + 
    xlim(as.numeric(xlim1), as.numeric(xlim2)) + 
    geom_hline(yintercept= -log10(0.05), color = "black", linetype = "dashed") +
    geom_vline(xintercept = log2(1), color = "black",linetype = "dashed") +
    
    ## Add points that have significant p-values in positive fold change
    geom_point(data = positive_fold, aes(x = log2fc, y = -log10(p_value)), size = 1.5) +
    geom_text_repel(data = positive_fold[1:(min(nrow(positive_fold), 5)), ], 
                    aes(x = log2fc, y = -log10(p_value), label= Row.names), 
                    max.overlaps = 15) +
    
    ## Add points that have significant p-values in negative fold change
    geom_point(data = negative_fold, aes(x = log2fc, y = -log10(p_value)), size = 1.5) +
    geom_text_repel(data = negative_fold[1:(min(nrow(negative_fold), 5)), ], 
                    aes(x = log2fc, y = -log10(p_value), label = Row.names), 
                    max.overlaps = 15) +
    
    labs(x = paste0("log2fold-change of relative abundance (", label1, "/", label2, ")"), 
         y = "-log10(p-value)") + ggtitle("") + 
    
    # Define color mapping
    scale_color_manual(values = setNames(c("brown2", "#45B1E9"), c(g1, g2)), 
                       labels = c(g1, g2)) +
    theme_bw()
  
  pdf(file = paste0("volcano_plot_", label1, "_vs_", label2, ".pdf"), height = 7, width = 8)
  print(plot)
  dev.off()
  
  return(list(positive_fold, negative_fold, fc_data_t_p))
  
}

####################################################################

# LBD vs. Control

####################################################################

# 1. Create the folder if not exist
folder_name <- "1-lbd_vs_control_diff_abundance_results"
folder <- file.path(".", folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

# 2. Extract data for LBD and controls
lbd_vs_control_relab <- species_all[species_all$condition %in% c("lbd", "lbd_control"), ]

# 3. Use arcsine square root transformation for the data
lbd_vs_control_relab_transformed <- trs(lbd_vs_control_relab)

# 4. Find the bacterial species that are passing the prevalence cut-off
lbd_vs_control_relab_transformed_filtered <- prevalence_cutoff(lbd_vs_control_relab_transformed, metadata)

# 5. Use the mixed-effect linear model to get p values condition variable. household_id will be random effect, and age and BMI are confounding variables. 
lbd_vs_control_relab_results <- mixed_effect(lbd_vs_control_relab_transformed_filtered, "./1-lbd_vs_control_diff_abundance_results/lbd_vs_control_differential_abundance_p_values.csv")

# 6. Use median of each group to calculate the log2fc
lbd_vs_control_diff_abun_median <- log2_fc_results(lbd_vs_control_relab_transformed_filtered, "lbd_BIOME", "lbd_control", lbd_vs_control_relab_results, "./1-lbd_vs_control_diff_abundance_results/lbd_vs_control_diff_abun_p_fc_done.csv")

# 7. Only focus on the differentially abundant species that have log2fc of median relative abundance not 0, and p < 0.05
sig_lbd_vs_control_diff_abun_species <- lbd_vs_control_diff_abun_median[[1]][lbd_vs_control_diff_abun_median[[1]]$p_value < 0.05 & lbd_vs_control_diff_abun_median[[1]]$log2fc != 0, ]
nrow(sig_lbd_vs_control_diff_abun_species)

# 8. Create a new column in the significant species data frame to indicate whether the species is more abundant in LBD or Control based on the log2 fold-change. If the log2 fold-change is greater than 0, it indicates that the species is more abundant in LBD; otherwise, it is more abundant in Control.
sig_lbd_vs_control_diff_abun_species$group <- ifelse(sig_lbd_vs_control_diff_abun_species$log2fc>0, "LBD>Control", "LBD<Control")

## 8.1 Count the number of species that are more abundant in LBD
nrow(sig_lbd_vs_control_diff_abun_species[sig_lbd_vs_control_diff_abun_species$group == "LBD>Control", ])

## 8.2 Count the number of species that are more abundant in Control
nrow(sig_lbd_vs_control_diff_abun_species[sig_lbd_vs_control_diff_abun_species$group == "LBD<Control", ])

# 9. Create the volcano plot
lbd_vs_control_fc_results <- volcano_plot(lbd_vs_control_diff_abun_median[[2]], "LBD", "Control", -12, 12) 

# 10. Create the overall boxplot
# 10.1 Select species that are used for boxplots
positive_fold <- lbd_vs_control_fc_results[[1]][1:4, ]
negative_fold <- lbd_vs_control_fc_results[[2]][1:5, ]

# 10.2 Order species based on median relative abundance in each group and combine the ordered species names for plotting
order_of_species_pos <- positive_fold[order(positive_fold$median_g1, decreasing = T), ]$Row.names
order_of_species_neg <- negative_fold[order(negative_fold$median_g2, decreasing = T), ]$Row.names
order_of_species <- c(order_of_species_pos, order_of_species_neg)

# 10.3 Combine the results from positive and negative fold changes for boxplot visualization and filter the data to include only the selected species for plotting.
boxplot_data <- rbind(lbd_vs_control_fc_results[[1]], lbd_vs_control_fc_results[[2]])
boxplot_data_done <- boxplot_data[boxplot_data$Row.names %in% order_of_species, ]

# 10.4 Melt the data frame to long format for ggplot2 visualization, create a new column to indicate the condition (LBD or Control) based on the variable names, and set the factor levels for species and condition to ensure proper ordering in the boxplot.
boxplot_long_df <- melt(boxplot_data_done[, 1:51], id.vars = "Row.names")
boxplot_long_df$condition <- ifelse(grepl("control", boxplot_long_df$variable), "Control", "LBD")
colnames(boxplot_long_df)[1] <- "species"
boxplot_long_df$species <- factor(boxplot_long_df$species, levels = order_of_species)
boxplot_long_df$condition <- factor(boxplot_long_df$condition, levels = c("LBD", "Control"))

# 10.5 Create the boxplot for the selected species
pdf(file = "./1-lbd_vs_control_diff_abundance_results/lbd_vs_control_boxplot_species.pdf", height = 8, width = 10)
plot <- ggplot(boxplot_long_df, 
                     aes(x = species, y = value, fill = condition)) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white", 
               outlier.size = 2, 
               outlier.shape = 21) +
  scale_fill_manual(values = c("#ff9274", "#55b7e6")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Relative abundance (arcsine square root transformed)")
print(plot)
dev.off()

# 10. Calculate Cohen's d
# 10.1 Create a new data frame that includes only the relevant columns for calculating Cohen's d, including the 'condition' column and any other non-microbe columns.
cohen_df <- lbd_vs_control_relab_transformed_filtered[, c(1:493,495)]

# 10.2 Create a list of microbe column names to be used for calculating Cohen's d. This list includes all columns except for the 'household_id' and 'condition' columns, which are not relevant for the effect size calculation.
p_list <- colnames(lbd_vs_control_relab_transformed_filtered)[1:492]

# 10.3 Create a vector of microbe column names by excluding the 'household_id' and 'condition' columns from the list of all column names in the 'cohen_df' data frame.
microbe_cols <- setdiff(names(cohen_df), c("household_id", "condition"))

# 10.4 Calculate Cohen's d for each microbe
cohen_results <- lapply(microbe_cols, function(microbe) {
  
  df_wide <- cohen_df %>%
    select(household_id, condition, all_of(microbe)) %>%
    pivot_wider(names_from = condition, values_from = all_of(microbe)) %>%
    mutate(diff = lbd - lbd_control)
  
  d <- cohens_d(df_wide$diff, mu = 0, paired = FALSE)  # one-sample on diff scores
  
  data.frame(
    microbe = microbe,
    cohens_d = d$Cohens_d,
    ci_low   = d$CI_low,
    ci_high  = d$CI_high
  )}) %>%
  
  bind_rows()

# 10.5 Merge the Cohen's d results with the log2 fold-change results.
cohens_results_all <- merge(lbd_vs_control_fc_results[[3]], cohen_results, by.x = "Row.names", by.y = "microbe", all = TRUE)

# 10.6 Create a new data frame that includes the species, p-value, log2 fold-change, and Cohen's d results for reporting.
cohens_results_all <- cohens_results_all %>%
  mutate(p_effect_size = (-log10(cohens_results_all$p_value))*abs(cohens_results_all$cohens_d)) %>%
  arrange(desc(p_effect_size))
write.csv(cohens_results_all, file = "./1-lbd_vs_control_diff_abundance_results/cohensD_lbd_vs_control_species.csv", row.names = F)

# 10.7 Select species that are used for boxplots and order them based on the previously defined order of species.
sig_cohen_results <- cohen_results[cohen_results$microbe %in% order_of_species, ]
sig_cohen_results$microbe <- factor(sig_cohen_results$microbe, levels = order_of_species)

# 10.8 Create a bar plot of Cohen's d for each microbe
pdf(file = "./1-lbd_vs_control_diff_abundance_results/cohensD_LBD_species.pdf", height = 3.5, width = 8)
ggplot(sig_cohen_results, aes(x = microbe, y = cohens_d)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  ylim(-0.8, 0.8) + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Cohen's d")
dev.off()

####################################################################

# iRBD vs. Control

####################################################################

# 1. Create the folder if not exist
folder_name <- "2-irbd_vs_control_diff_abundance_results"
folder <- file.path(".", folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

# 2. Extract data for iRBD and their controls
irbd_vs_control_relab <- species_all[species_all$condition %in% c("irbd", "irbd_control"), ]

# 3. Use arcsine square root transformation for the data
irbd_vs_control_relab_transformed <- trs(irbd_vs_control_relab)

# 4. Find the bacterial species that are passing the prevalence cut-off
irbd_vs_control_relab_transformed_filtered <- prevalence_cutoff(irbd_vs_control_relab_transformed, metadata)

# 5. Use the mixed-effect linear model to get p values condition variable. household_id will be random effect, and age and BMI are confounding variables. 
irbd_vs_control_relab_results <- mixed_effect(irbd_vs_control_relab_transformed_filtered, "./2-irbd_vs_control_diff_abundance_results/irbd_vs_control_differential_abundance_p_values.csv")

# 6. Use median of each group to calculate the log2fc
irbd_vs_control_diff_abun_median <- log2_fc_results(irbd_vs_control_relab_transformed_filtered, "irbd_BIOME", "irbd_control", irbd_vs_control_relab_results, "./2-irbd_vs_control_diff_abundance_results/irbd_vs_control_diff_abun_p_fc_done.csv")

# 7. Only focus on the differentially abundant species that have log2fc of median relative abundance not 0, and p < 0.05
sig_irbd_vs_control_diff_abun_species <- irbd_vs_control_diff_abun_median[[1]][irbd_vs_control_diff_abun_median[[1]]$p_value < 0.05 & irbd_vs_control_diff_abun_median[[1]]$log2fc != 0, ]
nrow(sig_irbd_vs_control_diff_abun_species)

# 8. Create a new column in the significant species data frame to indicate whether the species is more abundant in iRBD or Control based on the log2 fold-change. If the log2 fold-change is greater than 0, it indicates that the species is more abundant in iRBD; otherwise, it is more abundant in Control.
sig_irbd_vs_control_diff_abun_species$group <- ifelse(sig_irbd_vs_control_diff_abun_species$log2fc > 0, "irbd>Control", "irbd<Control")

## 8.1 Count the number of species that are more abundant in iRBD
nrow(sig_irbd_vs_control_diff_abun_species[sig_irbd_vs_control_diff_abun_species$group == "irbd>Control", ])

## 8.2 Count the number of species that are more abundant in Control
nrow(sig_irbd_vs_control_diff_abun_species[sig_irbd_vs_control_diff_abun_species$group == "irbd<Control", ])

# 9. Create the volcano plot
irbd_vs_control_fc_results <- volcano_plot(irbd_vs_control_diff_abun_median[[2]], "iRBD", "Control", -12, 12) 

# 10. Create the overall boxplot
# 10.1 Select species that are used for boxplots
positive_fold <- irbd_vs_control_fc_results[[1]][1:5, ]
negative_fold <- irbd_vs_control_fc_results[[2]][1:4, ]

# 10.2 Order species based on median relative abundance in each group and combine the ordered species names for plotting
order_of_species_pos <- positive_fold[order(positive_fold$median_g1, decreasing = T), ]$Row.names
order_of_species_neg <- negative_fold[order(negative_fold$median_g2, decreasing = T), ]$Row.names
order_of_species <- c(order_of_species_pos, order_of_species_neg)

# 10.3 Combine the results from positive and negative fold changes for boxplot visualization and filter the data to include only the selected species for plotting.
boxplot_data <- rbind(irbd_vs_control_fc_results[[1]], irbd_vs_control_fc_results[[2]])
boxplot_data_done <- boxplot_data[boxplot_data$Row.names %in% order_of_species, ]

# 10.4 Melt the data frame to long format for ggplot2 visualization, create a new column to indicate the condition (iRBD or Control) based on the variable names, and set the factor levels for species and condition to ensure proper ordering in the boxplot.
boxplot_long_df <- melt(boxplot_data_done[, 1:21], id.vars = "Row.names")
boxplot_long_df$condition <- ifelse(grepl("control", boxplot_long_df$variable), "Control", "iRBD")
colnames(boxplot_long_df)[1] <- "species"
boxplot_long_df$species <- factor(boxplot_long_df$species, levels = order_of_species)
boxplot_long_df$condition <- factor(boxplot_long_df$condition, levels = c("iRBD", "Control"))

# 10.5 Create the boxplot for the selected species
pdf(file = "./2-irbd_vs_control_diff_abundance_results/irbd_vs_control_boxplot_species.pdf", height = 8, width = 10)
plot <- ggplot(boxplot_long_df, 
                     aes(x = species, y = value, fill = condition)) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white", 
               outlier.size = 2, 
               outlier.shape = 21) +
  scale_fill_manual(values = c("#fdc848", "#2DA248")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Relative abundance (arcsine square root transformed)")
print(plot)
dev.off()

# 11. Calculate Cohen's d
# 11.1 Create a new data frame that includes only the relevant columns for calculating Cohen's d, including the 'condition' column and household columns.
cohen_df <- irbd_vs_control_relab_transformed_filtered[, c(1:590,592)]

# 11.2 Create a list of microbe column names to be used for calculating Cohen's d. This list includes all columns except for the 'household_id' and 'condition' columns, which are not relevant for the effect size calculation.
p_list <- colnames(irbd_vs_control_relab_transformed_filtered)[1:589]

# 11.3 Create a vector of microbe column names by excluding the 'household_id' and 'condition' columns from the list of all column names in the 'cohen_df' data frame.
microbe_cols <- setdiff(names(cohen_df), c("household_id", "condition"))

# 11.4 Calculate Cohen's d for each microbe
cohen_results <- lapply(microbe_cols, function(microbe) {
  
  df_wide <- cohen_df %>%
    select(household_id, condition, all_of(microbe)) %>%
    pivot_wider(names_from = condition, values_from = all_of(microbe)) %>%
    mutate(diff = irbd - irbd_control)
  
  d <- cohens_d(df_wide$diff, mu = 0, paired = FALSE)  # one-sample on diff scores
  
  data.frame(
    microbe = microbe,
    cohens_d = d$Cohens_d,
    ci_low   = d$CI_low,
    ci_high  = d$CI_high
  )}) %>%
  
  bind_rows()

# 10.5 Merge the Cohen's d results with the log2 fold-change results.
cohens_results_all <- merge(irbd_vs_control_fc_results[[3]], cohen_results, by.x = "Row.names", by.y = "microbe", all = TRUE)

# 10.6 Create a new data frame that includes the species, p-value, log2 fold-change, and Cohen's d results for reporting.
cohens_results_all <- cohens_results_all %>%
  mutate(p_effect_size = (-log10(cohens_results_all$p_value))*abs(cohens_results_all$cohens_d)) %>%
  arrange(desc(p_effect_size))
write.csv(cohens_results_all, file = "./2-irbd_vs_control_diff_abundance_results/cohensD_irbd_vs_control_species.csv", row.names = F)

# 10.7 Select species that are used for boxplots and order them based on the previously defined order of species.
sig_cohen_results <- cohen_results[cohen_results$microbe %in% order_of_species, ]
sig_cohen_results$microbe <- factor(sig_cohen_results$microbe, levels = order_of_species)

# 10.8 Create a bar plot of Cohen's d for each microbe
pdf(file = "./2-irbd_vs_control_diff_abundance_results/cohensD_irbd_vs_control_species.pdf", height = 3.5, width = 8)
ggplot(sig_cohen_results, aes(x = microbe, y = cohens_d)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  ylim(-2, 2) + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Cohen's d")
dev.off()