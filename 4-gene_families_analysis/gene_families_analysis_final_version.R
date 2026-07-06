library(readxl)
library(dplyr)
library(writexl)
library(lmerTest)
library(tidyr)

###############################################################################
folder_path <- "~/Desktop/LBD/analysis/gene_families/"
setwd(folder_path)
###############################################################################

gene_families <- as.data.frame(read_excel(path = "lbd_vs_control_selected_gene_families.xlsx", sheet = "gene_families"))
rownames(gene_families) <- gene_families$UniRef90_ID
gene_families <- gene_families[, -78]

colnames(gene_families) <- gsub("_S.*", "", colnames(gene_families))

keys <- as.data.frame(read_excel(path = "lbd_vs_control_selected_gene_families.xlsx", sheet = "keys"))

## Check the original data
colSums(gene_families)

# 6. Find the appropriate presence cutoff for bacterial species data
relab_number <- unlist(gene_families, use.name = FALSE)
relab_number_ordered <- relab_number[order(relab_number, decreasing = T)]

pdf("rank_plot_with_threshold_gene_families_lbd_vs_control.pdf", height = 6, width = 8)
plot(x=c(1:length(relab_number_ordered)), y=log10(relab_number_ordered), pch = 20, xlab = "Rank", ylab = "Relative abundance of bacterial species", main = "Rank-plot of ordered relative abundance of bacterial species", xlim = c(0, 170000))
abline(h=-7, col = "orange", lwd = 3)
abline(h=-7.2, col = "blue")
abline(h=-7.5, col = "orange")
dev.off()

## It seems 10^-7 is the good cutoff

# 7. Use prevalence cut-off to remove bacterial species that are present in very low abundance
## Number of low abundance species in each sample
cut_off <- 10^-7
apply(gene_families, 2, function(x) sum(as.numeric(x) <= cut_off, na.rm = TRUE))
apply(gene_families, 2, function(x) sum(as.numeric(x) <= 0, na.rm = TRUE))
## If use presence cutoff, how many low abundance species are needed to be removed
apply(gene_families, 2, function(x) sum(as.numeric(x) <= cut_off, na.rm = TRUE)) - apply(gene_families, 2, function(x) sum(as.numeric(x) <= 0, na.rm = TRUE))

gene_families_cutoff <- gene_families
gene_families_cutoff[gene_families_cutoff <= cut_off] <- 0

## Number of absence species in each sample
apply(gene_families_cutoff, 2, function(x) sum(as.numeric(x) == 0, na.rm = TRUE))

metadata <- read.csv(file = "~/Desktop/LBD/analysis/imputed_BMI_metadata.csv", row.names = 1)

gene_families_t <- data.frame(t(gene_families_cutoff))
gene_families_done <- merge(gene_families_t, metadata, by = "row.names")
rownames(gene_families_done) <- gene_families_done$Row.names
gene_families_done <- gene_families_done[, -1]

###############################################################################

## The following function is designed to filter out gene families (columns) with low prevalence from a given data set.
remove_low_prev_gene_families <- function(data, group){
  
  filtered_data <- data[, colSums(data > 0) >= 0.1*nrow(data)]
  result_df <- merge(filtered_data, group, by = "row.names", all = F)
  rownames(result_df) <- result_df$Row.names
  result_df <- result_df[, -1]
  
  return(result_df)
  
}

## The following function performs a transformation on a dataset and merges it with metadata and group information to produce a final structured data frame.

trs <- function(dataset, group) {
  
  transformed <- asin(sqrt(dataset[, 1:(ncol(dataset)-33)]))
  result_df <- merge(transformed, group, by = "row.names", all = F)
  rownames(result_df) <- result_df$Row.names
  result_df <- result_df[, -1]
  
  return(result_df)
  
} 

mixed_effect_gene_families <- function(filtered_data) {
  
  p_diff_relab <- c()
  
  for (i in 1:(ncol(filtered_data)-33)) {
    p_diff_relab[i] <- tryCatch({summary(lmer(filtered_data[, i] ~ condition + age + BMI + (1|household_id), data = filtered_data, REML = F))[["coefficients"]][2,5]}, error = function(e) NA) # If error occurs, assign NA and continue
  }
  
  results <- data.frame("gene_families" = names(filtered_data)[1:(ncol(filtered_data)-33)], "p_value" = p_diff_relab)
  
  results$q_value <- p.adjust(results$p_value, method = "BH")
  
  return(results)
  
}

log2_fc_gene_families_results <- function(filtered_data, group1, group2, p_val_results, file_name) {

  # Re-name the column names so it contains the group labels
  fc_data <- filtered_data[, c(1:(ncol(filtered_data)-33), ncol(filtered_data)-30)]
  fc_data$new_id <- paste(fc_data$condition, rownames(fc_data), sep = "_")
  rownames(fc_data) <- fc_data$new_id
  fc_data <- fc_data[, 1:(ncol(filtered_data)-33)]
  fc_data_t <- as.data.frame(t(fc_data))
  
  # Compute medians for each group
  fc_data_t$median_g1 <- apply(fc_data_t[, grepl(group1, colnames(fc_data_t))], 1, function(x) median(x))
  fc_data_t$median_g2 <- apply(fc_data_t[, grepl(group2, colnames(fc_data_t))], 1, function(x) median(x))
  
  # Compute log2 fold change
  pseudocount <- cut_off
  fc_data_t$log2fc <- log2((fc_data_t$median_g1 + pseudocount)/(fc_data_t$median_g2 + pseudocount))
  
  # Ensure pathway names from p_val and log2fc results are the same
  if (!all(rownames(fc_data_t) %in% p_val_results$gene_families)) {
    stop("Error: Some pathway in `fc_data_t` are missing in `p_results`.")
  }
  
  # Ensure pathway there is a column named "p_value" in the p_val results
  if (!"p_value" %in% colnames(p_val_results)) {
    stop("Error: `p_results` must contain a `p_value` column.")
  }
  
  # Merge p-values with log2fc results
  fc_data_t_p <- merge(fc_data_t, p_val_results, by.x = "row.names", by.y = "gene_families", all = F)
  fc_data_t_p_clean <- fc_data_t_p[, c(1, (ncol(fc_data_t_p)-4):ncol(fc_data_t_p))]
  colnames(fc_data_t_p_clean)[1] <- "gene_families"
  
  write.csv(fc_data_t_p_clean, file = file_name, row.names = F)
  return(list(fc_data_t_p_clean, fc_data_t_p))
  
}

###############################################################################

# LBD vs. Control 

# 1. Create the folder if not exist
folder_name <- "lbd_vs_control_gene_families"
folder <- file.path(folder_path, folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

# 2. Extract data for LBD and Control group
lbd_vs_control <- gene_families_done[gene_families_done$condition == "lbd" | gene_families_done$condition == "lbd_control", ]

# 3. Use arcsine square root transformation for the data
lbd_vs_control_transformed <- trs(lbd_vs_control, metadata)

# 4. Find the bacterial gene families that are passing the cut-off
lbd_vs_control_transformed_filtered <- remove_low_prev_gene_families(lbd_vs_control_transformed[, 1:(ncol(lbd_vs_control_transformed)-33)], metadata)

# 5. Use the mixed-effect linear model to get differentially abundant gene families 
## lbd_vs_control_gene_families_p_result <- mixed_effect_gene_families(lbd_vs_control_transformed_filtered)
## write.csv(lbd_vs_control_gene_families_p_result, file = paste0(folder, "/lbd_vs_control_gene_families_p_all.csv"), row.names = F)
lbd_vs_control_gene_families_p_result <- read.csv(file = paste0(folder, "/lbd_vs_control_gene_families_p_all.csv"))

# 6. Use mean and median of each group to calculate the log2fc

# ## 1) Re-name the row names of filtered data so it contains the group labels
# lbd_vs_control_fc_df <- lbd_vs_control_transformed_filtered[, c(1:(ncol(lbd_vs_control_transformed_filtered)-33), ncol(lbd_vs_control_transformed_filtered)-30)]
# lbd_vs_control_fc_df$new_id <- paste(lbd_vs_control_fc_df$condition, rownames(lbd_vs_control_fc_df), sep = "_")
# rownames(lbd_vs_control_fc_df) <- lbd_vs_control_fc_df$new_id
# lbd_vs_control_fc_df <- lbd_vs_control_fc_df[, 1:(ncol(lbd_vs_control_fc_df)-2)]
# lbd_vs_control_fc_t <- as.data.frame(t(lbd_vs_control_fc_df))
# 
# ## 2) Compute means for each group
# lbd_vs_control_fc_t$mean_g1 <- apply(lbd_vs_control_fc_t[, grepl("lbd_BIOME", colnames(lbd_vs_control_fc_t))], 1, function(x) mean(x))
# lbd_vs_control_fc_t$mean_g2 <- apply(lbd_vs_control_fc_t[, grepl("lbd_control", colnames(lbd_vs_control_fc_t))], 1, function(x) mean(x))
# 
# ## 3) Compute medians for each group
# lbd_vs_control_fc_t$median_g1 <- apply(lbd_vs_control_fc_t[, grepl("lbd_BIOME", colnames(lbd_vs_control_fc_t))], 1, function(x) median(x))
# lbd_vs_control_fc_t$median_g2 <- apply(lbd_vs_control_fc_t[, grepl("lbd_control", colnames(lbd_vs_control_fc_t))], 1, function(x) median(x))
# 
# ## 4) Compute log2 fold changes
# pseudocount <- cut_off
# lbd_vs_control_fc_t$log2fc_mean <- log2((lbd_vs_control_fc_t$mean_g1 + pseudocount)/(lbd_vs_control_fc_t$mean_g2 + pseudocount))
# lbd_vs_control_fc_t$log2fc_median <- log2((lbd_vs_control_fc_t$median_g1 + pseudocount)/(lbd_vs_control_fc_t$median_g2 + pseudocount))
# 
# ## 5. Ensure gene families names from p-value results and log2fc results are the same
# if (!all(rownames(lbd_vs_control_fc_t) %in% lbd_vs_control_gene_families_p_result$gene_families)) {
#   stop("Error: Some pathway in log2fc results are missing in p-value results.")
# }
# 
# ## 6. Merge p-values with log2fc results
# lbd_vs_control_fc_p <- merge(lbd_vs_control_fc_t, lbd_vs_control_gene_families_p_result, by.x = "row.names", by.y = "gene_families", all = T)
# rownames(lbd_vs_control_fc_p) <- lbd_vs_control_fc_p$Row.names
# lbd_vs_control_fc_p <- lbd_vs_control_fc_p[, -1]
# write.csv(lbd_vs_control_fc_p, file = paste0(folder, "/lbd_vs_control_gene_families_fc_p_results.csv"))

lbd_vs_control_fc_p <- read.csv(file = paste0(folder, "/lbd_vs_control_gene_families_fc_p_results.csv"))
colnames(lbd_vs_control_fc_p)[1] <- "gene_families"

lbd_vs_control_fc_p_pwy <- lbd_vs_control_fc_p %>%
  left_join(keys, by = c("gene_families" = "UniRef90_ID"))

# 7. Find significant gene families

## sig_gene_families_lbd_vs_control <- lbd_vs_control_fc_p[lbd_vs_control_fc_p$p_value < 0.05, ]
## write.csv(sig_gene_families_lbd_vs_control, file = paste0(folder, "/sig_lbd_vs_control_gene_families_fc_p_results.csv"))
sig_gene_families_lbd_vs_control <- read.csv(file = paste0(folder, "/sig_lbd_vs_control_gene_families_fc_p_results.csv"))
colnames(sig_gene_families_lbd_vs_control)[1] <- "gene_families"

# 8. Combine significant gene families with their associated pathways
sig_gene_families_lbd_vs_control_pwy <- sig_gene_families_lbd_vs_control %>%
  left_join(keys, by = c("gene_families" = "UniRef90_ID"))

# 9. Combine significant results with EC names
# EC_mapping <- read.delim(file = "map_level4ec_uniref90.txt", header = F, sep = "\t", stringsAsFactors = F)
# 
# ## Modify EC table so each uniref_id will match with one enzyme
# 
# EC_lookup_table <- EC_mapping %>%
#   pivot_longer(cols = -V1,  # All columns except the EC number
#                names_to = "source_column",
#                values_to = "UniRef90_ID") %>%
#   select(UniRef90_ID, EC_ID = V1)
# saveRDS(EC_lookup_table, file = "EC_lookup_table.rds")
EC_lookup_table <- readRDS("EC_lookup_table.rds")
# 
# sig_gene_families_lbd_vs_control_pwy_ec <- merge(sig_gene_families_lbd_vs_control_pwy, EC_lookup_table, by.x = "gene_families", by.y = "UniRef90_ID", all = F)
# write.csv(sig_gene_families_lbd_vs_control_pwy_ec, file = paste0(folder, "/sig_lbd_vs_control_gene_families_fc_p_pwy_ec_results.csv"))

lbd_vs_control_fc_p_pwy_ec <- merge(lbd_vs_control_fc_p_pwy, EC_lookup_table, by.x = "gene_families", by.y = "UniRef90_ID", all = F)
write.csv(lbd_vs_control_fc_p_pwy_ec, file = paste0(folder, "/all_lbd_vs_control_gene_families_fc_p_pwy_ec_results.csv"))

sig_gene_families_lbd_vs_control_pwy_ec <- read.csv(file = paste0(folder, "/sig_lbd_vs_control_gene_families_fc_p_pwy_ec_results.csv"), row.names = 1)

# 10. PWY-6731

pwy6731_df <- sig_gene_families_lbd_vs_control_pwy_ec[sig_gene_families_lbd_vs_control_pwy_ec$Pathway == "PWY-6731", c(1, 52:61)]

all_pwy6731 <- lbd_vs_control_fc_p_pwy_ec[lbd_vs_control_fc_p_pwy_ec$Pathway == "PWY-6731", c(1, 52:61)]
write.csv(all_pwy6731, file = "all_pwy6731.csv")


## based on mean

pwy6731_df_high_lbd <- pwy6731_df[pwy6731_df$log2fc_mean > 0, ]
write.csv(pwy6731_df_high_lbd, file = paste0(folder, "/pwy6731_high_lbd_mean.csv"))
pwy6731_df_higher_control <- pwy6731_df[pwy6731_df$log2fc_mean < 0, ]
write.csv(pwy6731_df_higher_control, file = paste0(folder, "/pwy6731_high_control_mean.csv"))

## based on median

pwy6731_df_high_lbd_median <- pwy6731_df[pwy6731_df$log2fc_median > 0, ]
write.csv(pwy6731_df_high_lbd_median, file = paste0(folder, "/pwy6731_high_lbd_median.csv"))
pwy6731_df_higher_control_median <- pwy6731_df[pwy6731_df$log2fc_median < 0, ]
write.csv(pwy6731_df_higher_control_median, file = paste0(folder, "/pwy6731_high_control_median.csv"))

# 11. PWY-7456

pwy7456_df <- sig_gene_families_lbd_vs_control_pwy_ec[sig_gene_families_lbd_vs_control_pwy_ec$Pathway == "PWY-7456", c(1, 52:61)]

## based on mean

pwy7456_df_high_lbd <- pwy7456_df[pwy7456_df$log2fc_mean > 0, ]
write.csv(pwy7456_df_high_lbd, file = paste0(folder, "/pwy7456_high_lbd_mean.csv"))
pwy7456_df_higher_control <- pwy7456_df[pwy7456_df$log2fc_mean < 0, ]
write.csv(pwy7456_df_higher_control, file = paste0(folder, "/pwy7456_high_control_mean.csv"))

## based on median

pwy7456_df_high_lbd_median <- pwy7456_df[pwy7456_df$log2fc_median > 0, ]
write.csv(pwy7456_df_high_lbd_median, file = paste0(folder, "/pwy7456_high_lbd_median.csv"))
pwy7456_df_higher_control_median <- pwy7456_df[pwy7456_df$log2fc_median < 0, ]
write.csv(pwy7456_df_higher_control_median, file = paste0(folder, "/pwy7456_high_control_median.csv"))

# 12. PWY-5030

pwy5030_df <- sig_gene_families_lbd_vs_control_pwy_ec[sig_gene_families_lbd_vs_control_pwy_ec$Pathway == "PWY-5030", c(1, 52:61)]

## based on mean

pwy5030_df_high_lbd <- pwy5030_df[pwy5030_df$log2fc_mean > 0, ]
write.csv(pwy5030_df_high_lbd, file = paste0(folder, "/pwy5030_high_lbd_mean.csv"))
pwy5030_df_higher_control <- pwy5030_df[pwy5030_df$log2fc_mean < 0, ]
write.csv(pwy5030_df_higher_control, file = paste0(folder, "/pwy5030_high_control_mean.csv"))

## based on median

pwy5030_df_high_lbd_median <- pwy5030_df[pwy5030_df$log2fc_median > 0, ]
write.csv(pwy5030_df_high_lbd_median, file = paste0(folder, "/pwy5030_high_lbd_median.csv"))
pwy5030_df_higher_control_median <- pwy5030_df[pwy5030_df$log2fc_median < 0, ]
write.csv(pwy5030_df_higher_control_median, file = paste0(folder, "/pwy5030_high_control_median.csv"))

###############################################################################

# irbd vs. Control 

# 1. Create the folder if not exist
folder_name <- "irbd_vs_control_gene_families"
folder <- file.path(folder_path, folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

# 2. Extract data for iRBD and Control group
irbd_vs_control <- gene_families_done[gene_families_done$condition == "irbd" | gene_families_done$condition == "irbd_control", ]

# 3. Use arcsine square root transformation for the data
irbd_vs_control_transformed <- trs(irbd_vs_control, metadata)

# 4. Find the bacterial gene families that are passing the cut-off
irbd_vs_control_transformed_filtered <- remove_low_prev_gene_families(irbd_vs_control_transformed[, 1:(ncol(irbd_vs_control_transformed)-33)], metadata)

# 5. Use the mixed-effect linear model to get differentially abundant gene families 
# irbd_vs_control_gene_families_p_result <- mixed_effect_gene_families(irbd_vs_control_transformed_filtered)
# write.csv(irbd_vs_control_gene_families_p_result, file = paste0(folder, "/irbd_vs_control_gene_families_p_all.csv"), row.names = F)
irbd_vs_control_gene_families_p_result <- read.csv(file = paste0(folder, "/irbd_vs_control_gene_families_p_all.csv"))

# 6. Use mean and median of each group to calculate the log2fc
# ## 1) Re-name the row names of filtered data so it contains the group labels
# irbd_vs_control_fc_df <- irbd_vs_control_transformed_filtered[, c(1:(ncol(irbd_vs_control_transformed_filtered)-33), ncol(irbd_vs_control_transformed_filtered)-30)]
# irbd_vs_control_fc_df$new_id <- paste(irbd_vs_control_fc_df$condition, rownames(irbd_vs_control_fc_df), sep = "_")
# rownames(irbd_vs_control_fc_df) <- irbd_vs_control_fc_df$new_id
# irbd_vs_control_fc_df <- irbd_vs_control_fc_df[, 1:(ncol(irbd_vs_control_fc_df)-2)]
# irbd_vs_control_fc_t <- as.data.frame(t(irbd_vs_control_fc_df))
# 
# ## 2) Compute means for each group
# irbd_vs_control_fc_t$mean_g1 <- apply(irbd_vs_control_fc_t[, grepl("irbd_BIOME", colnames(irbd_vs_control_fc_t))], 1, function(x) mean(x))
# irbd_vs_control_fc_t$mean_g2 <- apply(irbd_vs_control_fc_t[, grepl("irbd_control", colnames(irbd_vs_control_fc_t))], 1, function(x) mean(x))
# 
# ## 3) Compute medians for each group
# irbd_vs_control_fc_t$median_g1 <- apply(irbd_vs_control_fc_t[, grepl("irbd_BIOME", colnames(irbd_vs_control_fc_t))], 1, function(x) median(x))
# irbd_vs_control_fc_t$median_g2 <- apply(irbd_vs_control_fc_t[, grepl("irbd_control", colnames(irbd_vs_control_fc_t))], 1, function(x) median(x))
# 
# ## 4) Compute log2 fold changes
# pseudocount <- cut_off
# irbd_vs_control_fc_t$log2fc_mean <- log2((irbd_vs_control_fc_t$mean_g1 + pseudocount)/(irbd_vs_control_fc_t$mean_g2 + pseudocount))
# irbd_vs_control_fc_t$log2fc_median <- log2((irbd_vs_control_fc_t$median_g1 + pseudocount)/(irbd_vs_control_fc_t$median_g2 + pseudocount))
# 
# ## 5. Ensure gene families names from p-value results and log2fc results are the same
# if (!all(rownames(irbd_vs_control_fc_t) %in% irbd_vs_control_gene_families_p_result$gene_families)) {
#   stop("Error: Some pathway in log2fc results are missing in p-value results.")
# }
# 
# ## 6. Merge p-values with log2fc results
# irbd_vs_control_fc_p <- merge(irbd_vs_control_fc_t, irbd_vs_control_gene_families_p_result, by.x = "row.names", by.y = "gene_families", all = T)
# rownames(irbd_vs_control_fc_p) <- irbd_vs_control_fc_p$Row.names
# irbd_vs_control_fc_p <- irbd_vs_control_fc_p[, -1]
# write.csv(irbd_vs_control_fc_p, file = paste0(folder, "/irbd_vs_control_gene_families_fc_p_results.csv"))

irbd_vs_control_fc_p <- read.csv(file = paste0(folder, "/irbd_vs_control_gene_families_fc_p_results.csv"), row.names = 1)

# 7. Find significant gene families
# sig_gene_families_irbd_vs_control <- irbd_vs_control_fc_p[irbd_vs_control_fc_p$p_value < 0.05, ]
# write.csv(sig_gene_families_irbd_vs_control, file = paste0(folder, "/sig_irbd_vs_control_gene_families_fc_p_results.csv"))
sig_gene_families_irbd_vs_control <- read.csv(file = paste0(folder, "/sig_irbd_vs_control_gene_families_fc_p_results.csv"))
colnames(sig_gene_families_irbd_vs_control)[1] <- "gene_families"

# 8. Combine significant gene families with their associated pathways
sig_gene_families_irbd_vs_control_pwy <- sig_gene_families_irbd_vs_control %>%
  left_join(keys, by = c("gene_families" = "UniRef90_ID"))

# 9. Combine significant results with EC names
# EC_mapping <- read.delim(file = "map_level4ec_uniref90.txt", header = F, sep = "\t", stringsAsFactors = F)
# 
# ## Modify EC table so each uniref_id will match with one enzyme
# 
# EC_lookup_table <- EC_mapping %>%
#   pivot_longer(cols = -V1,  # All columns except the EC number
#                names_to = "source_column",
#                values_to = "UniRef90_ID") %>%
#   select(UniRef90_ID, EC_ID = V1)
# saveRDS(EC_lookup_table, file = "EC_lookup_table.rds")
EC_lookup_table <- readRDS("EC_lookup_table.rds")
# sig_gene_families_irbd_vs_control_pwy_ec <- merge(sig_gene_families_irbd_vs_control_pwy, EC_lookup_table, by.x = "gene_families", by.y = "UniRef90_ID", all = F)
# write.csv(sig_gene_families_irbd_vs_control_pwy_ec, file = paste0(folder, "/sig_irbd_vs_control_gene_families_fc_p_pwy_ec_results.csv"))

irbd_vs_control_fc_p_ec <- merge(irbd_vs_control_fc_p, EC_lookup_table, by.x = "row.names", by.y = "UniRef90_ID", all = F)
write.csv(irbd_vs_control_fc_p_ec, file = paste0(folder, "/irbd_vs_control_gene_families_fc_p_pwy_ec_results.csv"))

sig_gene_families_irbd_vs_control_pwy_ec <- read.csv(file = paste0(folder, "/sig_irbd_vs_control_gene_families_fc_p_pwy_ec_results.csv"), row.names = 1)

# 10. PWY-6731

pwy6731_df <- sig_gene_families_irbd_vs_control_pwy_ec[sig_gene_families_irbd_vs_control_pwy_ec$Pathway == "PWY-6731", c(1, 52:61)]

## based on mean

pwy6731_df_high_irbd <- pwy6731_df[pwy6731_df$log2fc_mean > 0, ]
write.csv(pwy6731_df_high_irbd, file = paste0(folder, "/pwy6731_high_irbd_mean.csv"))
pwy6731_df_higher_control <- pwy6731_df[pwy6731_df$log2fc_mean < 0, ]
write.csv(pwy6731_df_higher_control, file = paste0(folder, "/pwy6731_high_control_mean.csv"))

## based on median

pwy6731_df_high_irbd_median <- pwy6731_df[pwy6731_df$log2fc_median > 0, ]
write.csv(pwy6731_df_high_irbd_median, file = paste0(folder, "/pwy6731_high_irbd_median.csv"))
pwy6731_df_higher_control_median <- pwy6731_df[pwy6731_df$log2fc_median < 0, ]
write.csv(pwy6731_df_higher_control_median, file = paste0(folder, "/pwy6731_high_control_median.csv"))

# 11. PWY-7456

pwy7456_df <- sig_gene_families_irbd_vs_control_pwy_ec[sig_gene_families_irbd_vs_control_pwy_ec$Pathway == "PWY-7456", c(1, 52:61)]

## based on mean

pwy7456_df_high_irbd <- pwy7456_df[pwy7456_df$log2fc_mean > 0, ]
write.csv(pwy7456_df_high_irbd, file = paste0(folder, "/pwy7456_high_irbd_mean.csv"))
pwy7456_df_higher_control <- pwy7456_df[pwy7456_df$log2fc_mean < 0, ]
write.csv(pwy7456_df_higher_control, file = paste0(folder, "/pwy7456_high_control_mean.csv"))

## based on median

pwy7456_df_high_irbd_median <- pwy7456_df[pwy7456_df$log2fc_median > 0, ]
write.csv(pwy7456_df_high_irbd_median, file = paste0(folder, "/pwy7456_high_irbd_median.csv"))
pwy7456_df_higher_control_median <- pwy7456_df[pwy7456_df$log2fc_median < 0, ]
write.csv(pwy7456_df_higher_control_median, file = paste0(folder, "/pwy7456_high_control_median.csv"))

# 12. PWY-5030

pwy5030_df <- sig_gene_families_irbd_vs_control_pwy_ec[sig_gene_families_irbd_vs_control_pwy_ec$Pathway == "PWY-5030", c(1, 22:31)]

## based on mean

pwy5030_df_high_irbd <- pwy5030_df[pwy5030_df$log2fc_mean > 0, ]
write.csv(pwy5030_df_high_irbd, file = paste0(folder, "/pwy5030_high_irbd_mean.csv"))
pwy5030_df_higher_control <- pwy5030_df[pwy5030_df$log2fc_mean < 0, ]
write.csv(pwy5030_df_higher_control, file = paste0(folder, "/pwy5030_high_control_mean.csv"))

## based on median

pwy5030_df_high_irbd_median <- pwy5030_df[pwy5030_df$log2fc_median > 0, ]
write.csv(pwy5030_df_high_irbd_median, file = paste0(folder, "/pwy5030_high_irbd_median.csv"))
pwy5030_df_higher_control_median <- pwy5030_df[pwy5030_df$log2fc_median < 0, ]
write.csv(pwy5030_df_higher_control_median, file = paste0(folder, "/pwy5030_high_control_median.csv"))