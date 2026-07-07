library(ggplot2)
library(dplyr)
library(reshape2)
library(lmerTest)
library(ggrepel)
library(ggpubr)

###############################################################################
setwd("../3-microbial_functional_pathway_analysis")
###############################################################################

# 1. Read the pathway relative abundance data
pathway_data <- read.csv(file = "pathway_results.tsv", sep = "\t", header = T, row.names = 1)

# 2. Modify sample names
colnames(pathway_data) <- gsub("_S.*", "", colnames(pathway_data))

# 3. Read the metadata
metadata <- read.csv(file = "../0-raw_data/imputed_BMI_metadata.csv", sep = ",", header = T, row.names = 1)
pathway_clean <- pathway_data[, c(rownames(metadata))]

# 4. Change the pathway data from percentage into proportion
colSums(pathway_clean)
pathway_prop <- data.frame(apply(pathway_clean, 2, function(x) x/sum(x)))
colSums(pathway_prop)

# 5. Find the appropriate presence cutoff for pathway data
relab_number <- unlist(pathway_prop, use.name = FALSE)
relab_number_ordered <- relab_number[order(relab_number, decreasing = T)]

pdf("rank_plot_with_threshold_bacterial_pathway_done.pdf", height = 6, width = 8)
plot(x=c(1:length(relab_number_ordered)), y=log10(relab_number_ordered), pch = 20, xlab = "Rank", ylab = "Relative abundance of bacterial pathway", main = "Rank-plot of ordered relative abundance of bacterial pathway", xlim = c(0, 24000))
abline(h=-5, col = "red")
abline(h=-4.5, col = "orange", lwd = 3)
abline(h=-4, col = "orange")
dev.off()

## It seems 10^-4.5 is the good cutoff

# 6. Use prevalence cut-off to remove bacterial pathway that are present in very low abundance
## Number of low abundance pathway in each sample
cut_off <- 10^-4.5
apply(pathway_prop, 2, function(x) sum(as.numeric(x) <= cut_off, na.rm = TRUE))
apply(pathway_prop, 2, function(x) sum(as.numeric(x) <= 0, na.rm = TRUE))
## If use presence cutoff, how many low abundance pathway are needed to be removed
apply(pathway_prop, 2, function(x) sum(as.numeric(x) <= cut_off, na.rm = TRUE)) - apply(pathway_prop, 2, function(x) sum(as.numeric(x) <= 0, na.rm = TRUE))

pathway_cutoff <- pathway_prop
pathway_cutoff[pathway_cutoff < cut_off] <- 0

## Number of absence pathway in each sample
apply(pathway_cutoff, 2, function(x) sum(as.numeric(x) == 0, na.rm = TRUE))

# 8. Add metadata to the proportion data
pathway_t <- as.data.frame(t(pathway_cutoff))
pathway_all <- merge(pathway_t, metadata, by = "row.names", all = F)
rownames(pathway_all) <- pathway_all$Row.names
pathway_all <- pathway_all[, -1]

###############################################################################

## The following function performs a transformation on a dataset and merges it with metadata and group information to produce a final structured data frame.

trs <- function(dataset, metadata) {
  
  transformed <- asin(sqrt(dataset[, 1:(ncol(dataset)-33)]))
  result_df <- merge(transformed, metadata, by = "row.names", all = F)
  rownames(result_df) <- result_df$Row.names
  result_df <- result_df[, -1]
  
  return(result_df)
  
} 

# Find the bacterial pathway that are passing the prevalence cut-off

prevalence_cutoff <- function(transformed_data, metadata) {
  
  filtered_data <- transformed_data[, c(1:(ncol(transformed_data)-33))][, colSums(transformed_data[, c(1:(ncol(transformed_data)-33))] > 0) >= 0.1*nrow(transformed_data)]
  result_df <- merge(filtered_data, metadata, by = "row.names", all = F)
  rownames(result_df) <- result_df$Row.names
  result_df <- result_df[, -1]
  
  return(result_df)
  
}

# The following function performs a mixed-effects model analysis to compare the abundance of pathway between two time points. It calculates statistical significance (p-value) and trend direction (difference) for each pathway across the given time points.

mixed_effect <- function(filtered_data, file_name) {
  
  p_diff_relab <- c()
  
  for (i in 1:(ncol(filtered_data)-33)) {
    
    p_diff_relab[i] <- tryCatch({summary(lmer(filtered_data[, i] ~ condition + age + BMI + (1|household_id), data = filtered_data, REML = F))[["coefficients"]][2,5]}, error = function(e) NA)
    # If error occurs, assign NA and continue
    
  }
  
  results <- data.frame("pathway" = names(filtered_data)[1:(ncol(filtered_data)-33)], "p_value" = p_diff_relab)
  results$q_value <- p.adjust(results$p_value, method = "BH")
  
  write.csv(results, file = file_name, row.names = F)
  return(results)
  
}

# 

log2_fc_results <- function(filtered_data, group1, group2, p_val_results, file_name) {
  
  # Re-name the column names so it contains the group labels
  fc_data <- filtered_data[, c(1:(ncol(filtered_data)-33), ncol(filtered_data)-30)]
  fc_data$new_id <- paste(fc_data$condition, rownames(fc_data), sep = "_")
  rownames(fc_data) <- fc_data$new_id
  fc_data <- fc_data[, 1:(ncol(filtered_data)-33)]
  fc_data_t <- as.data.frame(t(fc_data))
  
  # Compute medians for each group
  fc_data_t$median_g1 <- apply(fc_data_t[, grepl(group1, colnames(fc_data_t))], 1, function(x) median(x))
  fc_data_t$median_g2 <- apply(fc_data_t[, grepl(group2, colnames(fc_data_t))], 1, function(x) median(x))
  
  # Compute mean for each group
  fc_data_t$mean_g1 <- apply(fc_data_t[, grepl(group1, colnames(fc_data_t))], 1, function(x) mean(x))
  fc_data_t$mean_g2 <- apply(fc_data_t[, grepl(group2, colnames(fc_data_t))], 1, function(x) mean(x))
  
  # Compute log2 fold change
  pseudocount <- cut_off
  fc_data_t$log2fc <- log2((fc_data_t$median_g1 + pseudocount)/(fc_data_t$median_g2 + pseudocount))
  
  # Compute log2 mean fold change
  pseudocount <- cut_off
  fc_data_t$log2fc_mean <- log2((fc_data_t$mean_g1 + pseudocount)/(fc_data_t$mean_g2 + pseudocount))
  
  # Compute mean difference
  fc_data_t$mean_difference <- fc_data_t$mean_g1 - fc_data_t$mean_g2
  
  # Ensure pathway names from p_val and log2fc results are the same
  if (!all(rownames(fc_data_t) %in% p_val_results$pathway)) {
    stop("Error: Some pathway in `fc_data_t` are missing in `p_results`.")
  }
  
  # Ensure pathway there is a column named "p_value" in the p_val results
  if (!"p_value" %in% colnames(p_val_results)) {
    stop("Error: `p_results` must contain a `p_value` column.")
  }
  
  # Merge p-values with log2fc results
  fc_data_t_p <- merge(fc_data_t, p_val_results, by.x = "row.names", by.y = "pathway", all = F)
  fc_data_t_p_clean <- fc_data_t_p[, c(1, (ncol(fc_data_t_p)-8):ncol(fc_data_t_p))]
  colnames(fc_data_t_p_clean)[1] <- "pathways"
  
  write.csv(fc_data_t_p_clean, file = file_name, row.names = F)
  return(list(fc_data_t_p_clean, fc_data_t_p))
  
}

## The following function processes the filtered dataset containing pathway abundance values for two conditions, computes the log2 fold-change, integrates p-values from mixed-effects linear model, and categorizes significant differences, and generates a volcano plot. 

volcano_plot <- function(fc_data_t_p, label1, label2, xlim1, xlim2) {
  
  # Define categories
  g1 <- paste0(label1, ">", label2)
  g2 <- paste0(label1, "<", label2)
  
  fc_data_t_p$category <- "All"
  fc_data_t_p$category[fc_data_t_p$log2fc > 0 & fc_data_t_p$p_value < 0.05] <- g1
  fc_data_t_p$category[fc_data_t_p$log2fc < 0 & fc_data_t_p$p_value < 0.05] <- g2
  fc_data_t_p$category <- factor(fc_data_t_p$category, levels = c("All", g1, g2))
  
  # Extract significant data and make them in an increasing order based on p_value
  positive_fold <- fc_data_t_p %>% filter(category == g1) %>% arrange(p_value)
  negative_fold <- fc_data_t_p %>% filter(category == g2) %>% arrange(p_value)
  
  write.csv(positive_fold, file = paste0("pathway_higher_in_", label1, "_than_", label2, ".csv"), row.names = F)
  write.csv(negative_fold, file = paste0("pathway_higher_in_", label2, "_than_", label1, ".csv"), row.names = F)
  
  plot <- ggplot(data = fc_data_t_p, aes(x = log2fc, y = -log10(p_value), color = category)) + 
    geom_point(size = 1, color = "lightgrey") + 
    xlim(as.numeric(xlim1), as.numeric(xlim2)) + 
    geom_hline(yintercept= -log10(0.05), color = "black", linetype = "dashed") +
    geom_vline(xintercept = log2(1), color = "black",linetype = "dashed") +
    
    ## Add points that have significant q-values in positive fold change
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
  
  pdf(file = paste0("volcano_plot_", label1, "_vs_", label2, "_pathway.pdf"), height = 7, width = 8)
  print(plot)
  dev.off()
  
  return(list(positive_fold, negative_fold, fc_data_t_p))
  
}

# Create the overall boxplot
overall_boxplot <- function(fc_results, grepl_term, group1, group2, file_name, color1, color2) {
  
  positive_fold <- fc_results[[1]][1:(min(nrow(fc_results[[1]]), 5)), ]
  negative_fold <- fc_results[[2]][1:(min(nrow(fc_results[[2]]), 5)), ]
  order_of_pathway_pos <- positive_fold[order(positive_fold$median_g1, decreasing = T), ]$Row.names
  order_of_pathway_neg <- negative_fold[order(negative_fold$median_g2, decreasing = T), ]$Row.names
  order_of_pathway <- c(order_of_pathway_pos, order_of_pathway_neg)
  boxplot_data <- rbind(fc_results[[1]], fc_results[[2]])

  boxplot_long_df <- melt(boxplot_data[, 1:(ncol(boxplot_data)-6)], id.vars = "Row.names")
  colnames(boxplot_long_df)[1] <- "pathway"
  boxplot_long_df$condition <- ifelse(grepl(grepl_term, boxplot_long_df$variable), group2, group1)
  boxplot_long_df$condition <- factor(boxplot_long_df$condition, levels = c(group1, group2))
  boxplot_long_df_done <- boxplot_long_df[boxplot_long_df$pathway %in% order_of_pathway, ]
  boxplot_long_df_done$pathway <- factor(boxplot_long_df_done$pathway, levels = order_of_pathway)
  boxplot_long_df_done$condition <- factor(boxplot_long_df_done$condition, levels = c(group1, group2))
  
  pdf(file = paste0("../3-microbial_functional_pathway_analysis/", file_name, "_boxplot_pathways.pdf"), height = 8, width = 10)
  plot <- ggplot(boxplot_long_df_done, 
                 aes(x = pathway, y = value, fill = condition)) + 
    geom_boxplot(outlier.color = "black", 
                 outlier.fill = "white", 
                 outlier.size = 2, 
                 outlier.shape = 21) +
    scale_fill_manual(values = c(color1, color2)) +
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
  
  return(boxplot_long_df_done)
  
}

# individual_boxplot generates individual boxplots for each specified pathway and saves each plot as a separate PDF file.
individual_boxplot <- function(boxplot_file, color1, color2, file_name) {
  
  for (i in 1:length(unique(boxplot_file[, 1]))) {
    plot <- ggplot(boxplot_file[boxplot_file$pathway == boxplot_file[i, 1], ], 
                   aes(x = condition, y = value, fill = condition)) + 
      geom_boxplot(outlier.shape = NA) +
      geom_point(shape = 21, size = 2, position = position_jitter(width = 0.05)) +
      scale_fill_manual(values = c(color1, color2)) +
      theme_bw() + 
      theme(
        panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
      labs(x = "", 
           y = "Relative abundance (arcsine square root transformed)",
           title = boxplot_file[i, 1])
    
    pdf(file = paste0("pathway_analysis/", file_name, "_boxplot_", boxplot_file[i, 1], ".pdf"), height = 6, width = 4)
    print(plot)
    dev.off()
  }
}




fisher_exact_test <- function(filtered_data, file_name, name1, name2) {
  
  fisher_data <- filtered_data[, c(1:(ncol(filtered_data)-33), ncol(filtered_data)-30)]
  
  p <- c()
  species <- c()
  group1 <- c()
  group2 <- c()
  
  for (i in 1:(ncol(fisher_data)-1)) {
    
    species[i] <- colnames(fisher_data)[i]
    fisher_table <- as.matrix(table(fisher_data[, i], fisher_data$condition))
    perc <- sweep(fisher_table, 2, colSums(fisher_table), FUN = "/")
    
    if (length(unique(fisher_data[, i])) < 2) {
      p[i] <- "NA"
      group1[i] <- perc[1,1]
      group2[i] <- perc[1,2]
    } else {
      p[i] <- fisher.test(fisher_table)[["p.value"]]
      group1[i] <- perc[2,1]
      group2[i] <- perc[2,2]
    }
  }
  
  prev_results <- data.frame(species, p, group1, group2)
  colnames(prev_results) <- c("pathway", "p_value", name1, name2)
  write.csv(prev_results, file = file_name, row.names = F)
  
  return(prev_results)
  
}

heatmap_z <- function(sig, 
                      label1, 
                      label2, 
                      sig_file_name, 
                      heatmap_file_name, 
                      taxonomy_table = taxonomy) {
  
  # Create 'condition' column based on comparison of the 3rd and 4th columns
  sig$condition <- ifelse(sig[, 3] > sig[, 4], label1, label2)
  
  # Split and sort by descending abundance
  sig_1 <- sig[sig$condition == label1, ]
  sig_1 <- sig_1[order(-sig_1[, 3]), ]
  
  sig_2 <- sig[sig$condition == label2, ]
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

###############################################################################

### Differential abundance analysis ###

# LBD vs. Control#

# 1. Extract data for LBD and their controls
lbd_vs_control_relab <- pathway_all[pathway_all$condition %in% c("lbd", "lbd_control"), ]

# 2. Use arcsine square root transformation for the data
lbd_vs_control_relab_transformed <- trs(lbd_vs_control_relab, metadata)

# 3. Find the bacterial pathway that are passing the prevalence cut-off
lbd_vs_control_relab_transformed_filtered <- prevalence_cutoff(lbd_vs_control_relab_transformed, metadata)

# 4. Use the mixed-effect linear model to get p values. household_id will be random effect, and age and BMI are confounding variables. 
lbd_vs_control_diff_abun_p_results <- mixed_effect(lbd_vs_control_relab_transformed_filtered, "lbd_vs_control_diff_abun_p_values.csv")
# lbd_vs_control_diff_abun_p_results <- read.csv(file = "../3-microbial_functional_pathway_analysis/lbd_vs_control/lbd_vs_control_diff_abun_p_values.csv")
sum(lbd_vs_control_diff_abun_p_results$p_value < 0.01)
sum(lbd_vs_control_diff_abun_p_results$p_value < 0.05)

# 5. Use median of each group to calculate the log2fc
lbd_vs_control_diff_abun_median <- log2_fc_results(lbd_vs_control_relab_transformed_filtered, "lbd_BIOME", "lbd_control", lbd_vs_control_diff_abun_p_results, "lbd_vs_control_diff_abun_p_fc_done.csv")

# 6. Only focus on the differentially abundant pathway that have median relative abundance not 0, and p<0.05
sig_lbd_vs_control_diff_abun <- lbd_vs_control_diff_abun_median[[1]][lbd_vs_control_diff_abun_median[[1]]$p_value < 0.05 & lbd_vs_control_diff_abun_median[[1]]$log2fc != 0, ]
nrow(sig_lbd_vs_control_diff_abun)
sig_lbd_vs_control_diff_abun$group <- ifelse(sig_lbd_vs_control_diff_abun$log2fc>0, "LBD>Control", "LBD<Control")
# write.csv(sig_lbd_vs_control_diff_abun, file = "sig_lbd_vs_control_diff_abun_pathway.csv", row.names = F)

# 7. Create the volcano plot
lbd_vs_control_fc_results <- volcano_plot(lbd_vs_control_diff_abun_median[[2]], "LBD", "Control", -10, 10) 

# 8. Create the overall boxplot
positive_fold <- lbd_vs_control_fc_results[[1]][1:5, ]
negative_fold <- lbd_vs_control_fc_results[[2]][1:5, ]
order_of_pathway_pos <- positive_fold[order(positive_fold$median_g1, decreasing = T), ]$Row.names[1:4]
order_of_pathway_neg <- negative_fold[order(negative_fold$median_g2, decreasing = T), ]$Row.names[1:5]
order_of_pathway <- c(order_of_pathway_pos, order_of_pathway_neg)
boxplot_data <- rbind(lbd_vs_control_fc_results[[1]], lbd_vs_control_fc_results[[2]])
boxplot_data_done <- boxplot_data[boxplot_data$Row.names %in% order_of_pathway, ]

boxplot_long_df <- melt(boxplot_data_done[, 1:51], id.vars = "Row.names")
# write.csv(boxplot_long_df, "../3-microbial_functional_pathway_analysis/lbd_vs_control_boxplot_long_df.csv")
boxplot_long_df$condition <- ifelse(grepl("control", boxplot_long_df$variable), "Control", "LBD")
colnames(boxplot_long_df)[1] <- "pathway"
boxplot_long_df$pathway <- factor(boxplot_long_df$pathway, levels = order_of_pathway)
boxplot_long_df$condition <- factor(boxplot_long_df$condition, levels = c("LBD", "Control"))

pwy_5505_lbd <- boxplot_long_df[grepl("PWY-5505", boxplot_long_df$pathway), ]

pdf(file = paste0("../3-microbial_functional_pathway_analysis/lbd_vs_control_boxplot_pathways.pdf"), height = 8, width = 10)
plot <- ggplot(boxplot_long_df, 
               aes(x = pathway, y = value, fill = condition)) + 
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


# 10. Cohen's d

cohen_df <- lbd_vs_control_relab_transformed_filtered[, c(1:434, 436)]

# Assume your microbe columns are everything except household_id and condition
pathway_cols <- setdiff(names(cohen_df), c("household_id", "condition"))

# Calculate Cohen's d for each microbe
cohen_results <- lapply(pathway_cols, function(microbe) {
  
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

all_results_lbd_path <- merge(lbd_vs_control_fc_results[[3]], cohen_results, by.x = "Row.names", by.y = "microbe", all = TRUE)

all_results_lbd_path <- all_results_lbd_path %>%
  mutate(p_effect_size = (-log10(all_results_lbd_path$p_value))*abs(all_results_lbd_path$cohens_d)) %>%
  arrange(desc(p_effect_size))

sig_all_results_lbd <- all_results_lbd_path %>%
  mutate(p_effect_size = (-log10(all_results_lbd_path$p_value))*abs(all_results_lbd_path$cohens_d)) %>%
  arrange(desc(p_effect_size)) %>%
  filter(p_value < 0.05 & log2fc != 0)

write.csv(all_results_lbd_path, file = "all_results_lbd_path.csv", row.names = F)


cohen_results$microbe <- factor(cohen_results$microbe, levels = order_of_pathway)

pdf(file = "cohensD_lbd_pathway.pdf", height = 5, width = 8)
ggplot(cohen_results, aes(x = microbe, y = cohens_d)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  ylim(-0.7, 0.7) + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Cohen's d")
dev.off()

pdf(file = "cohensD_ci_irbd_species.pdf", height = 3.5, width = 8)
ggplot(cohen_results, aes(x = microbe, y = cohens_d)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() + 
  # ylim(-1.2, 1.2) + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Cohen's d")
dev.off()


# 9. Create boxplot for each significant pathway
individual_boxplot(boxplot_long_df, "#ff9274", "#55b7e6", "lbd_vs_control")

# 10. Calculate the difference between mean based on the relative abundance data before arcsine square-root transformed)
sig_lbd_vs_control_pathway_relab_og <- lbd_vs_control_relab[, colnames(lbd_vs_control_relab) %in% c(order_of_pathway, "condition")]

# Re-name the column names so it contains the group labels
sig_lbd_vs_control_pathway_relab_og$new_id <- paste(sig_lbd_vs_control_pathway_relab_og$condition, rownames(sig_lbd_vs_control_pathway_relab_og), sep = "_")
rownames(sig_lbd_vs_control_pathway_relab_og) <- sig_lbd_vs_control_pathway_relab_og$new_id
sig_lbd_vs_control_pathway_relab_og <- sig_lbd_vs_control_pathway_relab_og[, 1:(length(order_of_pathway))]
sig_lbd_vs_control_pathway_relab_og_t <- as.data.frame(t(sig_lbd_vs_control_pathway_relab_og))

# Compute mean for each group
sig_lbd_vs_control_pathway_relab_og_t$mean_lbd <- apply(sig_lbd_vs_control_pathway_relab_og_t[, grepl("lbd_BIOME", colnames(sig_lbd_vs_control_pathway_relab_og_t))], 1, function(x) mean(x))
sig_lbd_vs_control_pathway_relab_og_t$mean_control <- apply(sig_lbd_vs_control_pathway_relab_og_t[, grepl("control_BIOME", colnames(sig_lbd_vs_control_pathway_relab_og_t))], 1, function(x) mean(x))

sig_lbd_vs_control_pathway_relab_og_t$mean_diff <- sig_lbd_vs_control_pathway_relab_og_t$mean_lbd - sig_lbd_vs_control_pathway_relab_og_t$mean_control
sig_lbd_vs_control_pathway_relab_og_t$pathway <- rownames(sig_lbd_vs_control_pathway_relab_og_t)

sig_lbd_vs_control_pathway_relab_og_t$pathway <- factor(sig_lbd_vs_control_pathway_relab_og_t$pathway, levels = order_of_pathway)

pdf(file = "../3-microbial_functional_pathway_analysis/lbd_vs_control/mean_diff_bar_lbd_control_v2.pdf", height = 7, width = 8)
ggplot(sig_lbd_vs_control_pathway_relab_og_t, aes(x = as.factor(pathway), y = mean_diff)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "mean difference")
dev.off()

###############################################################################

# iRBD vs. Control#

# 1. Extract data for iRBD and their controls
irbd_vs_control_relab <- pathway_all[pathway_all$condition %in% c("irbd", "irbd_control"), ]

# 2. Use arcsine square root transformation for the data
irbd_vs_control_relab_transformed <- trs(irbd_vs_control_relab, metadata)

# 3. Find the bacterial pathway that are passing the prevalence cut-off
irbd_vs_control_relab_transformed_filtered <- prevalence_cutoff(irbd_vs_control_relab_transformed, metadata)

# 4. Use the mixed-effect linear model to get p values condition variable. household_id will be random effect, and age and BMI are confounding variables. 
irbd_vs_control_diff_abun_p_results <- mixed_effect(irbd_vs_control_relab_transformed_filtered, "irbd_vs_control_diff_abun_p_values.csv")
# irbd_vs_control_diff_abun_p_results <- read.csv(file = "irbd_vs_control_diff_abun_p_values.csv")
sum(irbd_vs_control_diff_abun_p_results$p_value < 0.01)
sum(irbd_vs_control_diff_abun_p_results$p_value < 0.05)

# 5. Use median of each group to calculate the log2fc
irbd_vs_control_diff_abun_median <- log2_fc_results(irbd_vs_control_relab_transformed_filtered, "irbd_BIOME", "irbd_control", irbd_vs_control_diff_abun_p_results, "irbd_vs_control_diff_abun_p_fc_done.csv")

# 6. Only focus on the differentially abundant pathway that have median relative abundance not 0, and p<0.05
sig_irbd_vs_control_diff_abun <- irbd_vs_control_diff_abun_median[[1]][irbd_vs_control_diff_abun_median[[1]]$p_value < 0.05 & irbd_vs_control_diff_abun_median[[1]]$log2fc != 0, ]
nrow(sig_irbd_vs_control_diff_abun)
sig_irbd_vs_control_diff_abun$group <- ifelse(sig_irbd_vs_control_diff_abun$log2fc>0, "iRBD>Control", "iRBD<Control")
# write.csv(sig_irbd_vs_control_diff_abun, file = "../3-microbial_functional_pathway_analysis/irbd_vs_control/sig_irbd_vs_control_diff_abun_pathway.csv", row.names = F)

# 7. Create the volcano plot
irbd_vs_control_fc_results <- volcano_plot(irbd_vs_control_diff_abun_median[[2]], "iRBD", "Control", -10, 10) 

# 8. Create the overall boxplot
positive_fold <- irbd_vs_control_fc_results[[1]][1:5, ]
negative_fold <- irbd_vs_control_fc_results[[2]][1:5, ]
order_of_pathway_pos <- positive_fold[order(positive_fold$median_g1, decreasing = T), ]$Row.names[1:5]
order_of_pathway_neg <- negative_fold[order(negative_fold$median_g2, decreasing = T), ]$Row.names[1:5]
order_of_pathway <- c(order_of_pathway_pos, order_of_pathway_neg)
boxplot_data <- rbind(irbd_vs_control_fc_results[[1]], irbd_vs_control_fc_results[[2]])
boxplot_data_done <- boxplot_data[boxplot_data$Row.names %in% order_of_pathway, ]

boxplot_long_df <- melt(boxplot_data_done[, 1:21], id.vars = "Row.names")
# write.csv(boxplot_long_df, "../3-microbial_functional_pathway_analysis/irbd_vs_control_boxplot_long_df.csv")
boxplot_long_df$condition <- ifelse(grepl("control", boxplot_long_df$variable), "Control", "irbd")
colnames(boxplot_long_df)[1] <- "pathway"
boxplot_long_df$pathway <- factor(boxplot_long_df$pathway, levels = order_of_pathway)
boxplot_long_df$condition <- factor(boxplot_long_df$condition, levels = c("irbd", "Control"))

pwy_5505_irbd_1 <- boxplot_data[grepl("PWY-5505", boxplot_data$Row.names), ]
pwy_5505_irbd_2 <- melt(pwy_5505_irbd_1[, 1:21], id.vars = "Row.names")
pwy_5505_irbd_2$condition <- ifelse(grepl("control", pwy_5505_irbd_2$variable), "irbd_control", "irbd")
colnames(pwy_5505_irbd_2)[1] <- "pathway"

pwy_5505 <- rbind(pwy_5505_lbd, pwy_5505_irbd_2)
pwy_5505$condition <- factor(pwy_5505$condition, levels = c("LBD", "Control", "irbd", "irbd_control"))

plot <- ggplot(pwy_5505, aes(x = condition, y = value, fill = condition)) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white", 
               outlier.size = 2, 
               outlier.shape = 21) +
  # geom_point(shape = 21, size = 2, position = position_jitter(width = 0.05)) +
  scale_fill_manual(values = c("#ff9274", "#55b7e6", "#2DA248", "#fdc848")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", 
       y = "Relative abundance (arcsine square root transformed)")

pdf(file = "../3-microbial_functional_pathway_analysis/pwy5505.pdf", height = 6, width = 4)
print(plot)
dev.off()


df_plot <- read.csv(file = "../3-microbial_functional_pathway_analysis/pwy-5030.csv")
df_plot$condition <- factor(df_plot$condition, levels = c("LBD", "LBD_Control", "iRBD", "iRBD_Control"))

plot <- ggplot(df_plot, aes(x = condition, y = value, fill = condition)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 21, size = 2, position = position_jitter(width = 0.05)) +
  scale_fill_manual(values = c("#ff9274", "#55b7e6", "#2DA248", "#fdc848")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", 
       y = "Relative abundance (arcsine square root transformed)")

pdf(file = "../3-microbial_functional_pathway_analysis/pwy5030.pdf", height = 6, width = 6)
print(plot)
dev.off()

log2fc_mean_df <- sig_irbd_vs_control_diff_abun[sig_irbd_vs_control_diff_abun$pathways %in% order_of_pathway, ]
log2fc_mean_df$pathways <- factor(log2fc_mean_df$pathways, levels = order_of_pathway)

pdf(file = "../3-microbial_functional_pathway_analysis/irbd_vs_control/mean_diff_bar_irbd.pdf", height = 5, width = 8)
ggplot(log2fc_mean_df, aes(x = as.factor(pathways), y = mean_difference)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "log2fc")
dev.off()

# 10. Cohen's d

cohen_df <- irbd_vs_control_relab_transformed_filtered[, c(1:440, 442)]

# Assume your microbe columns are everything except household_id and condition
pathway_cols <- setdiff(names(cohen_df), c("household_id", "condition"))

# Calculate Cohen's d for each microbe
cohen_results <- lapply(pathway_cols, function(microbe) {
  
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

all_results_irbd_path <- merge(irbd_vs_control_fc_results[[3]], cohen_results, by.x = "Row.names", by.y = "microbe", all = TRUE)

all_results_irbd_path <- all_results_irbd_path %>%
  mutate(p_effect_size = (-log10(all_results_irbd_path$p_value))*abs(all_results_irbd_path$cohens_d)) %>%
  arrange(desc(p_effect_size))

sig_all_results_irbd <- all_results_irbd_path %>%
  mutate(p_effect_size = (-log10(all_results_irbd_path$p_value))*abs(all_results_irbd_path$cohens_d)) %>%
  arrange(desc(p_effect_size)) %>%
  filter(p_value < 0.05 & log2fc > 0)

write.csv(all_results_irbd_path, file = "all_results_irbd_path.csv", row.names = F)

cohen_results$microbe <- factor(cohen_results$microbe, levels = order_of_pathway)

pdf(file = "cohensD_irbd_pathway.pdf", height = 5, width = 8)
ggplot(cohen_results, aes(x = microbe, y = cohens_d)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  ylim(-2.2, 2.2) + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Cohen's d")
dev.off()

pdf(file = "cohensD_ci_irbd_species.pdf", height = 3.5, width = 8)
ggplot(cohen_results, aes(x = microbe, y = cohens_d)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() + 
  # ylim(-1.2, 1.2) + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Cohen's d")
dev.off()


boxplot_df <- overall_boxplot(irbd_vs_control_fc_results, "control", "iRBD", "Control", "irbd_vs_control", "#2DA248", "#fdc848")

# 9. Create boxplot for each significant pathway
individual_boxplot(boxplot_df, "#2DA248", "#fdc848", "irbd_vs_control")

# 10. Calculate the difference between mean based on the relative abundance data before arcsine square-root transformed)
sig_irbd_vs_control_pathway_relab_og <- irbd_vs_control_relab[, colnames(irbd_vs_control_relab) %in% c(order_of_pathway, "condition")]

# Re-name the column names so it contains the group labels
sig_irbd_vs_control_pathway_relab_og$new_id <- paste(sig_irbd_vs_control_pathway_relab_og$condition, rownames(sig_irbd_vs_control_pathway_relab_og), sep = "_")
rownames(sig_irbd_vs_control_pathway_relab_og) <- sig_irbd_vs_control_pathway_relab_og$new_id
sig_irbd_vs_control_pathway_relab_og <- sig_irbd_vs_control_pathway_relab_og[, 1:(length(order_of_pathway))]
sig_irbd_vs_control_pathway_relab_og_t <- as.data.frame(t(sig_irbd_vs_control_pathway_relab_og))

# Compute mean for each group
sig_irbd_vs_control_pathway_relab_og_t$mean_irbd <- apply(sig_irbd_vs_control_pathway_relab_og_t[, grepl("irbd_BIOME", colnames(sig_irbd_vs_control_pathway_relab_og_t))], 1, function(x) mean(x))
sig_irbd_vs_control_pathway_relab_og_t$mean_control <- apply(sig_irbd_vs_control_pathway_relab_og_t[, grepl("control_BIOME", colnames(sig_irbd_vs_control_pathway_relab_og_t))], 1, function(x) mean(x))

sig_irbd_vs_control_pathway_relab_og_t$mean_diff <- sig_irbd_vs_control_pathway_relab_og_t$mean_irbd - sig_irbd_vs_control_pathway_relab_og_t$mean_control
sig_irbd_vs_control_pathway_relab_og_t$pathway <- rownames(sig_irbd_vs_control_pathway_relab_og_t)

sig_irbd_vs_control_pathway_relab_og_t$pathway <- factor(sig_irbd_vs_control_pathway_relab_og_t$pathway, levels = order_of_pathway)

pdf(file = "../3-microbial_functional_pathway_analysis/irbd_vs_control/mean_diff_bar_irbd_control_v2.pdf", height = 5, width = 8)
ggplot(sig_irbd_vs_control_pathway_relab_og_t, aes(x = as.factor(pathway), y = mean_diff)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "mean difference")
dev.off()

###############################################################################

# LBD vs. iRBD #

# 1. Extract data for LBD and their controls
lbd_vs_irbd_relab <- pathway_all[pathway_all$condition %in% c("lbd", "irbd"), ]

# 2. Use arcsine square root transformation for the data
lbd_vs_irbd_relab_transformed <- trs(lbd_vs_irbd_relab, metadata)

# 3. Find the bacterial pathway that are passing the prevalence cut-off
lbd_vs_irbd_relab_transformed_filtered <- prevalence_cutoff(lbd_vs_irbd_relab_transformed, metadata)

# 4. Use the regular linear model to get p values condition variable. age and BMI are confounding variables. 
p <- c()

for (i in 1:(ncol(lbd_vs_irbd_relab_transformed_filtered)-33)) {
  
  p[i] <- summary(lm(lbd_vs_irbd_relab_transformed_filtered[, i] ~ condition + age + BMI, data = lbd_vs_irbd_relab_transformed_filtered))[["coefficients"]][2,4]
  
}

p_results <- data.frame("pathway" = names(lbd_vs_irbd_relab_transformed_filtered)[1:(ncol(lbd_vs_irbd_relab_transformed_filtered)-33)], "p_value" = p)

p_results$q_value <- p.adjust(p_results$p_value, method = "BH")

lbd_vs_irbd_relab_results <- p_results

write.csv(lbd_vs_irbd_relab_results, file = "lbd_vs_irbd_diff_abun_p_values.csv", row.names = F)

sum(lbd_vs_irbd_relab_results$p_value < 0.01)
sum(lbd_vs_irbd_relab_results$p_value < 0.05)

# 5. Use median of each group to calculate the log2fc
lbd_vs_irbd_diff_abun_median <- log2_fc_results(lbd_vs_irbd_relab_transformed_filtered, "lbd_BIOME", "irbd_BIOME", lbd_vs_irbd_relab_results, "lbd_vs_irbd_diff_abun_p_fc_done.csv")

# 6. Only focus on the differentially abundant pathway that have median relative abundance not 0, and p<0.05
sig_lbd_vs_irbd_diff_abun <- lbd_vs_irbd_diff_abun_median[[1]][lbd_vs_irbd_diff_abun_median[[1]]$p_value < 0.05 & lbd_vs_irbd_diff_abun_median[[1]]$log2fc != 0, ]
nrow(sig_lbd_vs_irbd_diff_abun)
sig_lbd_vs_irbd_diff_abun$group <- ifelse(sig_lbd_vs_irbd_diff_abun$log2fc>0, "LBD>iRBD", "LBD<iRBD")
write.csv(sig_lbd_vs_irbd_diff_abun, file = "../3-microbial_functional_pathway_analysis/lbd_vs_control/sig_lbd_vs_irbd_diff_abun.csv", row.names = F)

# 7. Create the volcano plot
lbd_vs_irbd_fc_results <- volcano_plot(lbd_vs_irbd_diff_abun_median[[2]], "LBD", "iRBD", -10, 10) 

# 8. Create the overall boxplot
positive_fold <- lbd_vs_irbd_fc_results[[1]][1:5, ]
negative_fold <- lbd_vs_irbd_fc_results[[2]][1:5, ]
order_of_pathway_pos <- positive_fold[order(positive_fold$median_g1, decreasing = T), ]$Row.names[1:5]
order_of_pathway_neg <- negative_fold[order(negative_fold$median_g2, decreasing = T), ]$Row.names[1:5]
order_of_pathway <- c(order_of_pathway_pos, order_of_pathway_neg)
boxplot_data <- rbind(lbd_vs_irbd_fc_results[[1]], lbd_vs_irbd_fc_results[[2]])
boxplot_data_done <- boxplot_data[boxplot_data$Row.names %in% order_of_pathway, ]

boxplot_long_df <- melt(boxplot_data_done[, 1:(ncol(boxplot_data_done)-6)], id.vars = "Row.names")
write.csv(boxplot_long_df, "../3-microbial_functional_pathway_analysis/irbd_vs_control_boxplot_long_df.csv")
boxplot_long_df$condition <- ifelse(grepl("lbd", boxplot_long_df$variable), "LBD", "iRBD")
colnames(boxplot_long_df)[1] <- "pathway"
boxplot_long_df$pathway <- factor(boxplot_long_df$pathway, levels = order_of_pathway)
boxplot_long_df$condition <- factor(boxplot_long_df$condition, levels = c("LBD", "iRBD"))

pdf(file = "../3-microbial_functional_pathway_analysis/lbd_vs_irbd_boxplot.pdf", height = 6, width = 10)
plot <- ggplot(boxplot_long_df, 
               aes(x = pathway, y = value, fill = condition)) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white", 
               outlier.size = 2, 
               outlier.shape = 21) +
  # Add jitter if desired:
  # geom_point(aes(fill = condition), shape = 21, size = 1.5, 
  #            position = position_jitter(width = 0.2), color = "black") +
  scale_fill_manual(values = c("#ff9274", "#fec849")) +
  # ylim(as.numeric(0), as.numeric(0.5)) + 
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Relative abundance (arcsine square-root transformed)")
print(plot)
dev.off()

df_plot <- read.csv(file = "../3-microbial_functional_pathway_analysis/pwy-5030.csv")
df_plot$condition <- factor(df_plot$condition, levels = c("LBD", "LBD_Control", "iRBD", "iRBD_Control"))

plot <- ggplot(df_plot, aes(x = condition, y = value, fill = condition)) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white", 
               outlier.size = 2, 
               outlier.shape = 21) +
  # geom_point(shape = 21, size = 2, position = position_jitter(width = 0.05)) +
  scale_fill_manual(values = c("#ff9274", "#55b7e6", "#2DA248", "#fdc848")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", 
       y = "Relative abundance (arcsine square-root transformed)")

pdf(file = "../3-microbial_functional_pathway_analysis/pwy5030.pdf", height = 6, width = 4)
print(plot)
dev.off()

# 9. Create boxplot for each significant pathway
individual_boxplot(boxplot_df, "#ff9274", "#fdc848", "lbd_vs_irbd")

# 10. Calculate the difference between mean based on the relative abundance data before arcsine square-root transformed)
sig_lbd_vs_irbd_pathway_relab_og <- lbd_vs_irbd_relab[, colnames(lbd_vs_irbd_relab) %in% c(order_of_pathway, "condition")]

# Re-name the column names so it contains the group labels
sig_lbd_vs_irbd_pathway_relab_og$new_id <- paste(sig_lbd_vs_irbd_pathway_relab_og$condition, rownames(sig_lbd_vs_irbd_pathway_relab_og), sep = "_")
rownames(sig_lbd_vs_irbd_pathway_relab_og) <- sig_lbd_vs_irbd_pathway_relab_og$new_id
sig_lbd_vs_irbd_pathway_relab_og <- sig_lbd_vs_irbd_pathway_relab_og[, 1:(length(order_of_pathway))]
sig_lbd_vs_irbd_pathway_relab_og_t <- as.data.frame(t(sig_lbd_vs_irbd_pathway_relab_og))

# Compute mean for each group
sig_lbd_vs_irbd_pathway_relab_og_t$mean_lbd <- apply(sig_lbd_vs_irbd_pathway_relab_og_t[, grepl("lbd_BIOME", colnames(sig_lbd_vs_irbd_pathway_relab_og_t))], 1, function(x) mean(x))
sig_lbd_vs_irbd_pathway_relab_og_t$mean_irbd <- apply(sig_lbd_vs_irbd_pathway_relab_og_t[, grepl("irbd_BIOME", colnames(sig_lbd_vs_irbd_pathway_relab_og_t))], 1, function(x) mean(x))

sig_lbd_vs_irbd_pathway_relab_og_t$mean_diff <- sig_lbd_vs_irbd_pathway_relab_og_t$mean_lbd - sig_lbd_vs_irbd_pathway_relab_og_t$mean_irbd
sig_lbd_vs_irbd_pathway_relab_og_t$pathway <- rownames(sig_lbd_vs_irbd_pathway_relab_og_t)

sig_lbd_vs_irbd_pathway_relab_og_t$pathway <- factor(sig_lbd_vs_irbd_pathway_relab_og_t$pathway, levels = order_of_pathway)

pdf(file = "../3-microbial_functional_pathway_analysis/lbd_vs_irbd/mean_diff_bar_lbd_vs_irbd_v2.pdf", height = 5, width = 8)
ggplot(sig_lbd_vs_irbd_pathway_relab_og_t, aes(x = as.factor(pathway), y = mean_diff)) +
  geom_bar(stat = "identity", fill = "#6abd45") +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "mean difference")
dev.off()

###############################################################################

## Prevalence Analysis ##

prevalence_data <- as.data.frame(ifelse(pathway_cutoff > 0, 1, 0))

# 7. Add metadata to the prevalence data
prevalence_t <- as.data.frame(t(prevalence_data))
prevalence_all <- merge(prevalence_t, metadata, by = "row.names", all = F)
rownames(prevalence_all) <- prevalence_all$Row.names
prevalence_all <- prevalence_all[, -1]

###############################################################################

# LBD vs. Control#

# 1. Extract data for LBD and their controls
lbd_vs_control_prev <- prevalence_all[prevalence_all$condition %in% c("lbd", "lbd_control"), ]

# 2. Find the bacterial prevalence that are passing the cut-off
lbd_vs_control_prev_filtered <- prevalence_cutoff(lbd_vs_control_prev, metadata)

# 3. Fisher's exact test
lbd_vs_control_prev_results <- fisher_exact_test(lbd_vs_control_prev_filtered, "lbd_vs_control_prevalence_results_v2.csv", "LBD", "LBD_control")

###############################################################################

# iRBD vs. Control#

# 1. Extract data for iRBD and their controls
irbd_vs_control_prev <- prevalence_all[prevalence_all$condition %in% c("irbd", "irbd_control"), ]

# 2. Find the bacterial prevalence that are passing the cut-off
irbd_vs_control_prev_filtered <- prevalence_cutoff(irbd_vs_control_prev, metadata)

# 3. Fisher's exact test
irbd_vs_control_prev_results <- fisher_exact_test(irbd_vs_control_prev_filtered, "irbd_vs_control_prevalence_results_v2.csv", "iRBD", "iRBD_control")


bbb <- read.csv(file = "../3-microbial_functional_pathway_analysis/irbd_vs_control_prevalence_results.csv", row.names = 1)
bbb$qval <- p.adjust(bbb$p_value, method = "BH")

aaa <- read.csv(file = "../3-microbial_functional_pathway_analysis/lbd_vs_control_prevalence_results.csv", row.names = 1)
aaa$qval <- p.adjust(aaa$p_value, method = "BH")

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


###############################################################################

lbd_vs_control_prev_pathway_sig <- lbd_vs_control_prev_results[lbd_vs_control_prev_results$p_value < 0.05, ]
names(lbd_vs_control_prev_pathway_sig)[2] <- "lbd_vs_control_pval"
# write.csv(lbd_vs_control_prev_pathway_sig, file = "lbd_vs_control_prev_pathway_sig.csv")
lbd_vs_control_prev_pathway_sig <- read.csv(file = "../3-microbial_functional_pathway_analysis/lbd_vs_control_prev_pathway_sig.csv", row.names = 1)


irbd_vs_control_prev_pathway_sig <- irbd_vs_control_prev_results[irbd_vs_control_prev_results$p_value < 0.05, ]
# write.csv(irbd_vs_control_prev_pathway_sig, file = "irbd_vs_control_prev_pathway_sig.csv")

###############################################################################

### Differential prevalent species and pathways together ###

lbd_vs_control_prev_species_sig <- read.csv(file = "../3-microbial_functional_pathway_analysis/sig_lbd_vs_control_prev_new.csv")
lbd_vs_control_prev_species_sig <- lbd_vs_control_prev_species_sig[, 1:4]
colnames(lbd_vs_control_prev_species_sig)[1] <- "overall"
lbd_vs_control_prev_species_sig$c <- "species"
lbd_vs_control_prev_species_sig$group <- ifelse(lbd_vs_control_prev_species_sig$LBD>lbd_vs_control_prev_species_sig$LBD_control, "LBD>Control", "LBD<Control")
colnames(lbd_vs_control_prev_pathway_sig)[1] <- "overall"
lbd_vs_control_prev_pathway_sig$c <- "pathways"
lbd_vs_control_prev_pathway_sig$group <- ifelse(lbd_vs_control_prev_pathway_sig$LBD>lbd_vs_control_prev_pathway_sig$LBD_control, "LBD>Control", "LBD<Control")

lbd_vs_control_all <- rbind(lbd_vs_control_prev_species_sig, lbd_vs_control_prev_pathway_sig)
write.csv(lbd_vs_control_all, "lbd_vs_control_species_pathway_prevalence.csv", row.names = F)

pdf(file = "heatmap_all_lbd_vs_control_v2.pdf", width = 8, height = 5)
library(grid)

p <- pheatmap(data.matrix(lbd_vs_control_all[3:4]), 
              color = colorRampPalette(c("#f7f3e8", "#b1182d"))(100), 
              border_color = "white",
              cluster_rows = F, 
              cluster_cols = F, 
              cellwidth = 25, cellheight = 15,
              scale = "none", 
              main = "",
              angle_col = 0,
              labels_row = lbd_vs_control_all$overall,
              breaks = seq(0, 1, length.out = 101),
              legend_breaks = seq(0, 1, by = 0.2),
              legend_labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))
print(p)
dev.off()

irbd_vs_control_prev_species_sig <- read.csv(file = "../3-microbial_functional_pathway_analysis/irbd_vs_control/sig_results_irbd_vs_control.csv")
irbd_vs_control_prev_species_sig <- irbd_vs_control_prev_species_sig[, 1:4]
colnames(irbd_vs_control_prev_species_sig)[1] <- "overall"
irbd_vs_control_prev_species_sig$c <- "species"
irbd_vs_control_prev_species_sig$group <- ifelse(irbd_vs_control_prev_species_sig$iRBD>irbd_vs_control_prev_species_sig$iRBD_control, "iRBD>Control", "iRBD<Control")
colnames(irbd_vs_control_prev_pathway_sig)[1] <- "overall"
irbd_vs_control_prev_pathway_sig$c <- "pathways"
irbd_vs_control_prev_pathway_sig$group <- ifelse(irbd_vs_control_prev_pathway_sig$iRBD>irbd_vs_control_prev_pathway_sig$iRBD_control, "iRBD>Control", "iRBD<Control")

irbd_vs_control_all <- rbind(irbd_vs_control_prev_species_sig, irbd_vs_control_prev_pathway_sig)
write.csv(irbd_vs_control_all, "irbd_vs_control_species_pathway_prevalence.csv", row.names = F)

pdf(file = "heatmap_all_irbd_vs_control.pdf", width = 4, height = 4)
pheatmap(t(data.matrix(irbd_vs_control_all[3:4])), 
         color = colorRampPalette(c("#f7f3e8", "#b1182d"))(100), 
         border_color = "white",
         cluster_rows = F, 
         cluster_cols = F, 
         cellwidth = 25, cellheight = 15,
         scale = "none", 
         main = "",
         angle_col = 90,
         labels_col = irbd_vs_control_all$overall)
dev.off()

###############################################################################

### Gradual trend from control to iRBD to LBD ###

pathway_name <- lbd_vs_control_prev_pathway_sig$overall

lbd_df <- lbd_vs_control_prev_results[lbd_vs_control_prev_results$pathway %in% pathway_name, ]
colnames(lbd_df)[2] <- "lbd_p_value"
irbd_df <- irbd_vs_control_prev_results[irbd_vs_control_prev_results$pathway %in% pathway_name, ]
colnames(irbd_df)[2] <- "irbd_p_value"

lbd_irbd_df <- merge(lbd_df, irbd_df, by = "pathway", all = T)
pathway_lbd_irbd_df <- lbd_irbd_df[, c(1, 3, 4, 6, 7)]
colnames(pathway_lbd_irbd_df)[1] <- "overall"

species_lbd_irbd_df <- read.csv(file = "../3-microbial_functional_pathway_analysis/lbd_and_irbd_prev_species.csv")
colnames(species_lbd_irbd_df)[1] <- "overall"

all_lbd_irbd_df <- rbind(species_lbd_irbd_df, pathway_lbd_irbd_df)

pdf(file = "heatmap_figure_4.pdf", width = 6, height = 5)
pheatmap(data.matrix(all_lbd_irbd_df[2:5]), 
         color = colorRampPalette(c("#f7f3e8", "#b1182d"))(100), 
         breaks = seq(0, 1, length.out = 101),  # 0 to 1 scale with 100 color steps
         border_color = "white",
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         cellwidth = 25, cellheight = 15,
         scale = "none", 
         main = "",
         angle_col = 90,
         labels_row = all_lbd_irbd_df$overall)

dev.off()

###############################################################################

pdf(file = "PWY-6922.pdf", width = 4, height = 5.5)
ggplot(filter(selected_df_long, pathways == "PWY-6922: L-N&delta;-acetylornithine biosynthesis"), 
       aes(x = condition, y = percentage, fill = condition)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Prevalence (%)") +
  ylim(0,70) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.y = element_line(linetype = "dashed", color = "gray70", linewidth = 0.1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank()) + 
  scale_fill_manual(values = c("lbd" = "#d36480", 
                               "lbd_control" = "#e4c475",
                               "irbd" = "#5191b4",
                               "irbd_control" = "#9ac9db")) +
  rremove("legend")
dev.off()
