library(lmerTest)
library(vegan)
library(ade4)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(RColorBrewer)

####################################################################
folder_path <- "./1-alpha_beta_diversity"
setwd(folder_path)
####################################################################

# 1. Read the species relative abundance data
species_data <- read.table(file = "../0-raw_data/metaphlan_results_new.tsv", sep = "\t", header = T)

# 2. Filter the species-level data
species <- species_data[grepl("s__", species_data$clade_name) & !grepl("t__", species_data$clade_name), ]

# 3. Modify taxonomy names
rownames(species) <- species$clade_name
new_rownames <- gsub(".*s__", "", rownames(species))
rownames(species) <- new_rownames
species <- species[, -1]

# 4. Modify sample names
col_names <- names(species)
new_colnames <- gsub("metaphlan_|_S.*", "", col_names)
colnames(species) <- new_colnames

# 5. Read the metadata
metadata <- read.csv(file = "../0-raw_data/imputed_BMI_metadata.csv", sep = ",", header = T, row.names = 1)

# 6. Filter the species data to include only samples present in the metadata
species_clean <- species[, rownames(metadata)]

# 7. Change the taxonomic data from percentage into proportion
colSums(species_clean)
species_prop <- data.frame(apply(species_clean, 2, function(x) x/sum(x)))
colSums(species_prop)

# # 8. Find the appropriate presence cutoff for bacterial species data
# relab_number <- unlist(species_prop, use.name = FALSE)
# relab_number_ordered <- relab_number[order(relab_number, decreasing = T)]
# 
# pdf("rank_plot_with_threshold_bacterial_species_new.pdf", height = 6, width = 8)
# plot(x=c(1:length(relab_number_ordered)), y=log10(relab_number_ordered), pch = 20, xlab = "Rank", ylab = "Relative abundance of bacterial species", main = "Rank-plot of ordered relative abundance of bacterial species", xlim = c(0, 14000))
# abline(h=-5, col = "red")
# abline(h=-4, col = "blue")
# abline(h=-4.7, col = "orange")
# dev.off()

## It seems 10^-4.7 is the good cutoff

# 9. Use prevalence cut-off to remove bacterial species that are present in very low abundance
## Number of low abundance species in each sample
cut_off <- 10^-4.7

## If use presence cutoff, how many low abundance species are needed to be removed
apply(species_prop, 2, function(x) sum(as.numeric(x) <= cut_off, na.rm = TRUE)) - apply(species_prop, 2, function(x) sum(as.numeric(x) <= 0, na.rm = TRUE))

species_cutoff <- species_prop
species_cutoff[species_cutoff <= cut_off] <- 0

## Number of absence species in each sample
apply(species_cutoff, 2, function(x) sum(as.numeric(x) == 0, na.rm = TRUE))

# 10. Add metadata to the transformed data
species_t <- as.data.frame(t(species_cutoff))
species_all <- merge(species_t, metadata, by = "row.names", all = F)
rownames(species_all) <- species_all$Row.names
species_all <- species_all[, -1]

####################################################################

# LBD vs. Control

####################################################################

# 1. Create the folder if not exist
folder_name <- "lbd_vs_control_results"
folder <- file.path(".", folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

# 2. Select the samples with LBD and control from the species_all data frame
lbd_vs_control_data <- species_all[species_all$condition %in% c("lbd", "lbd_control"), ]

# 3. Calculate the Shannon index and species richness for each group
shannon <- diversity(lbd_vs_control_data[, c(1:(ncol(lbd_vs_control_data)-33))], index = "shannon", MARGIN = 1)
simpson <- diversity(lbd_vs_control_data[, c(1:(ncol(lbd_vs_control_data)-33))], index = "simpson", MARGIN = 1)
invsimpson <- diversity(lbd_vs_control_data[, c(1:(ncol(lbd_vs_control_data)-33))], index = "invsimpson", MARGIN = 1)
richness <- apply(lbd_vs_control_data[, c(1:(ncol(lbd_vs_control_data)-33))], 1, function(x) sum(x>0))

## Create a new data frame to store the diversity indices along with metadata
lbd_vs_control_diversity <- data.frame("condition" = lbd_vs_control_data$condition, "household_id" = lbd_vs_control_data$household_id, "age" = lbd_vs_control_data$age, "BMI" = lbd_vs_control_data$BMI, "sex" = lbd_vs_control_data$sex, "richness" = richness, "shannon" = shannon, "simpson" = simpson, "invsimpson" = invsimpson)

write.csv(lbd_vs_control_diversity, paste(folder, "/alpha_diversity_res_lbd_vs_control.csv", sep = ""), row.names = T)

# 4. Calculate p-values using mixed-effects linear models
shannon_p <- summary(lmer(shannon ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

simpson_p <- summary(lmer(simpson ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

invsimpson_p <- summary(lmer(invsimpson ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

richness_p <- summary(lmer(richness ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

mixed_effect_p_lbd_vs_control <- data.frame(diversity = c("shannon", "simpson", "invsimpson", "richness"), value = c(shannon_p, simpson_p, invsimpson_p, richness_p))

write.csv(mixed_effect_p_lbd_vs_control, paste0(folder, "/mixed_effect_p_lbd_vs_control.csv"))

# 5. Create long_df for boxplots
lbd_vs_control_alpha_diversity_long_df <- melt(lbd_vs_control_diversity, id.vars = "condition", measure.vars = c("shannon", "richness"))

# 6. PERMANOVA and PCoA analysis
# 6.1 Use arcsine square root transformation for relative abundance data
lbd_vs_control_transformed <- asin(sqrt(lbd_vs_control_data[, 1:(ncol(lbd_vs_control_data)-33)]))

# 6.2 Merge the transformed data with metadata for PERMANOVA analysis
lbd_vs_control_permanova <- merge(lbd_vs_control_transformed, metadata, by = "row.names", all = F)
rownames(lbd_vs_control_permanova) <- lbd_vs_control_permanova$Row.names
lbd_vs_control_permanova <- lbd_vs_control_permanova[, -1]

# 6.3 Subset the data for LBD samples only to check age distribution
lbd <- lbd_vs_control_permanova[lbd_vs_control_permanova$condition == "lbd", ]
summary(lbd$age)

# 6.4 Test for association between age, sex, BMI, and disease condition
## age: Wilcoxon rank-sum test
wilcox.test(age~condition, lbd_vs_control_permanova)
## sex: Fisher's exact test
print(table(lbd_vs_control_permanova$sex, lbd_vs_control_permanova$condition))
fisher.test(table(lbd_vs_control_permanova$sex, lbd_vs_control_permanova$condition))
## BMI: Mixed-effects linear model
summary(lmer(BMI ~ condition + (1|household_id), data = lbd_vs_control_permanova, REML = F))[["coefficients"]]

# 6.5 Perform PERMANOVA
## Remember to set random seed to ensure reproducibility of permutation-based p-values

## Condition only
set.seed(10)
adonis2(
  lbd_vs_control_permanova[, 1:(ncol(lbd_vs_control_permanova)-33)] ~ condition,
  data = lbd_vs_control_permanova,
  permutations = 999,
  strata = lbd_vs_control_permanova$household_id
)

## Age only
set.seed(10)
adonis2(
  lbd_vs_control_permanova[, 1:(ncol(lbd_vs_control_permanova)-33)] ~ age,
  data = lbd_vs_control_permanova,
  permutations = 999,
  strata = lbd_vs_control_permanova$household_id
)

## Sex only
set.seed(10)
adonis2(
  lbd_vs_control_permanova[, 1:(ncol(lbd_vs_control_permanova)-33)] ~ sex,
  data = lbd_vs_control_permanova,
  permutations = 999,
  strata = lbd_vs_control_permanova$household_id
)

## BMI only
set.seed(10)
adonis2(
  lbd_vs_control_permanova[, 1:(ncol(lbd_vs_control_permanova)-33)] ~ BMI,
  data = lbd_vs_control_permanova,
  permutations = 999,
  strata = lbd_vs_control_permanova$household_id
)

# 6.6 Perform PERMANOVA
## Calculate distance matrix
distance_matrix <- vegdist(lbd_vs_control_permanova[, 1:(ncol(lbd_vs_control_permanova)-33)], method = "bray")

## Perform PCoA analysis
pcoa_all <- dudi.pco(distance_matrix, scannf = FALSE, nf = 3)

## Extract the eigenvalues from the result of a PCoA
evals <- eigenvals(pcoa_all)
Variance <- evals / sum(evals)
Variance1 <- 100 * signif(Variance[1], 2)
Variance2 <- 100 * signif(Variance[2], 2)
Variance3 <- 100 * signif(Variance[3], 2)

## Extract first two PCoA axis values and store them into a dataframe
pc_plot_data <- data.frame(pcoa_all[["li"]][["A1"]], pcoa_all[["li"]][["A2"]])

## Add label to each point
pc_plot_data$Group <- as.factor(lbd_vs_control_permanova$condition)

## Rename the column name
colnames(pc_plot_data) <- c("pc1", "pc2", "Group")

## Calculate the centroid of each group
centroid_pc1 <- aggregate(pc1 ~ Group, data = pc_plot_data, FUN = mean)
centroid_pc2 <- aggregate(pc2 ~ Group, data = pc_plot_data, FUN = mean)
centroid <- merge(centroid_pc1, centroid_pc2, by = "Group")
colnames(centroid) <- c("Group", "c_pc1", "c_pc2")

## Merge into a new plot data frame
new_pc_pot <- merge(pc_plot_data, centroid, by = "Group")

## Generate the pcoa plot
pdf(file = paste(folder, "/pcoa_lbd_vs_control.pdf", sep = ""), width = 7, height = 4.5)
pcoa_plot <- ggplot(new_pc_pot, aes(x = c_pc1, y = c_pc2, color = Group)) + 
  scale_colour_manual(values = c("#ff9274", "#55b7e6")) + 
  theme(legend.title = element_blank()) + 
  labs(x = paste("PCoA1 (", Variance1, "%)", sep=""), 
       y = paste("PCoA2 (", Variance2, "%)", sep="")) + 
  ## Add the segments between centroid and each tip point
  geom_segment(aes(x = c_pc1, y = c_pc2, xend = pc1, yend = pc2)) + 
  # Add points for pc1 and pc2 with custom size
  geom_point(aes(x = pc1, y = pc2), size = 3, shape = 20) +
  # Add 95% confidence ellipses around groups
  stat_ellipse(data = new_pc_pot[, 1:3], aes(x = pc1, y = pc2, group = Group), level = 0.95) +
  theme_bw() + 
  theme(legend.title=element_blank(),
        panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
        panel.grid.minor = element_blank())
print(pcoa_plot)
dev.off()

####################################################################

# iRBD vs. Control

####################################################################

# 1. Create the folder if not exist
folder_name <- "irbd_vs_control_results"
folder <- file.path(".", folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

# 2. Select the samples with iRBD and control from the species_all data frame
irbd_vs_control_data <- species_all[species_all$condition %in% c("irbd", "irbd_control"), ]

# 3. Calculate the Shannon index and species richness for each group
shannon_diversity <- diversity(irbd_vs_control_data[, c(1:(ncol(irbd_vs_control_data)-33))], index = "shannon", MARGIN = 1)
simpson_diversity <- diversity(irbd_vs_control_data[, c(1:(ncol(irbd_vs_control_data)-33))], index = "simpson", MARGIN = 1)
invsimpson_diversity <- diversity(irbd_vs_control_data[, c(1:(ncol(irbd_vs_control_data)-33))], index = "invsimpson", MARGIN = 1)
richness <- apply(irbd_vs_control_data[, c(1:(ncol(irbd_vs_control_data)-33))], 1, function(x) sum(x>0))

## Create a new data frame to store the diversity indices along with metadata
irbd_vs_control_diversity <- data.frame("condition" = irbd_vs_control_data$condition, "household_id" = irbd_vs_control_data$household_id, "age" = irbd_vs_control_data$age, "BMI" = irbd_vs_control_data$BMI, "sex" = irbd_vs_control_data$sex, "richness" = richness, "shannon" = shannon_diversity, "simpson" = simpson_diversity, "invsimpson" = invsimpson_diversity)

write.csv(irbd_vs_control_diversity, paste(folder, "/alpha_diversity_results_irbd_vs_control.csv", sep = ""), row.names = T)

# 4. Calculate p-values using mixed-effects linear models
shannon_p <- summary(lmer(shannon ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p <- summary(lmer(simpson ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p <- summary(lmer(invsimpson ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p <- summary(lmer(richness ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

mixed_effect_p_irbd_vs_control <- data.frame(diversity = c("shannon", "simpson", "invsimpson", "richness"), value = c(shannon_p, simpson_p, invsimpson_p, richness_p))

write.csv(mixed_effect_p_irbd_vs_control, paste0(folder, "/mixed_effect_p_irbd_vs_control.csv"))

# 5. Create long_df for boxplots
irbd_vs_control_alpha_diversity_long_df <- melt(irbd_vs_control_diversity, id.vars = "condition", measure.vars = c("shannon", "richness"))

# 6. PERMANOVA and PCoA analysis
# 6.1 Use arcsine square root transformation for relative abundance data
irbd_vs_control_transformed <- asin(sqrt(irbd_vs_control_data[, 1:(ncol(irbd_vs_control_data)-33)]))
irbd_vs_control_permanova <- merge(irbd_vs_control_transformed, metadata, by = "row.names", all = F)
rownames(irbd_vs_control_permanova) <- irbd_vs_control_permanova$Row.names
irbd_vs_control_permanova <- irbd_vs_control_permanova[, -1]

# 6.2 Subset the data for iRBD samples only to check age distribution
irbd <- irbd_vs_control_permanova[irbd_vs_control_permanova$condition == "irbd", ]
summary(irbd$age)

# 6.3 Test for association between age, sex, BMI, and disease condition
## age: Wilcoxon rank-sum test
wilcox.test(age~condition, irbd_vs_control_permanova)
## sex: Fisher's exact test
print(table(irbd_vs_control_permanova$sex, irbd_vs_control_permanova$condition))
fisher.test(table(irbd_vs_control_permanova$sex, irbd_vs_control_permanova$condition))
## BMI: Mixed-effects linear model
summary(lmer(BMI ~ condition + (1|household_id), data = irbd_vs_control_permanova, REML = F))[["coefficients"]]

# 6.4 Perform PERMANOVA
## Remember to set random seed to ensure reproducibility of permutation-based p-values

## Condition only
set.seed(10)
adonis2(
  irbd_vs_control_permanova[, 1:(ncol(irbd_vs_control_permanova)-33)] ~ condition,
  data = irbd_vs_control_permanova,
  permutations = 999,
  strata = irbd_vs_control_permanova$household_id
)

## Age only
set.seed(10)
adonis2(
  irbd_vs_control_permanova[, 1:(ncol(irbd_vs_control_permanova)-33)] ~ age,
  data = irbd_vs_control_permanova,
  permutations = 999,
  strata = irbd_vs_control_permanova$household_id
)

## Sex only
set.seed(10)
adonis2(
  irbd_vs_control_permanova[, 1:(ncol(irbd_vs_control_permanova)-33)] ~ sex,
  data = irbd_vs_control_permanova,
  permutations = 999,
  strata = irbd_vs_control_permanova$household_id
)

## BMI only
set.seed(10)
adonis2(
  irbd_vs_control_permanova[, 1:(ncol(irbd_vs_control_permanova)-33)] ~ BMI,
  data = irbd_vs_control_permanova,
  permutations = 999,
  strata = irbd_vs_control_permanova$household_id
)

# 6.5 Perform PERMANOVA
## Calculate distance matrix
distance_matrix <- vegdist(irbd_vs_control_permanova[, 1:(ncol(irbd_vs_control_permanova)-33)], method = "bray")

## Perform PCoA analysis
pcoa_all <- dudi.pco(distance_matrix, scannf = FALSE, nf = 3)

## Extract the eigenvalues from the result of a PCoA
evals <- eigenvals(pcoa_all)
Variance <- evals / sum(evals)
Variance1 <- 100 * signif(Variance[1], 2)
Variance2 <- 100 * signif(Variance[2], 2)
Variance3 <- 100 * signif(Variance[3], 2)

## Extract first two PCoA axis values and store them into a dataframe
pc_plot_data <- data.frame(pcoa_all[["li"]][["A1"]], pcoa_all[["li"]][["A2"]])

## Add label to each point
pc_plot_data$Group <- as.factor(irbd_vs_control_permanova$condition)

## Rename the column name
colnames(pc_plot_data) <- c("pc1", "pc2", "Group")

## Calculate the centroid of each group
centroid_pc1 <- aggregate(pc1 ~ Group, data = pc_plot_data, FUN = mean)
centroid_pc2 <- aggregate(pc2 ~ Group, data = pc_plot_data, FUN = mean)
centroid <- merge(centroid_pc1, centroid_pc2, by = "Group")
colnames(centroid) <- c("Group", "c_pc1", "c_pc2")

## Merge into a new plot data frame
new_pc_pot <- merge(pc_plot_data, centroid, by = "Group")

## Generate the pcoa plot
pdf(file = paste(folder, "/pcoa_irbd_vs_control.pdf", sep = ""), width = 7, height = 4.5)
pcoa_plot <- ggplot(new_pc_pot, aes(x = c_pc1, y = c_pc2, color = Group)) + 
  scale_colour_manual(values = c("#fdc848", "#2DA248")) + 
  theme(legend.title = element_blank()) + 
  labs(x = paste("PCoA1 (", Variance1, "%)", sep=""), 
       y = paste("PCoA2 (", Variance2, "%)", sep="")) + 
  ## Add the segments between centroid and each tip point
  geom_segment(aes(x = c_pc1, y = c_pc2, xend = pc1, yend = pc2)) + 
  # Add points for pc1 and pc2 with custom size
  geom_point(aes(x = pc1, y = pc2), size = 3, shape = 20) +
  # Add 95% confidence ellipses around groups
  stat_ellipse(data = new_pc_pot[, 1:3], aes(x = pc1, y = pc2, group = Group), level = 0.95) +
  theme_bw() + 
  theme(legend.title=element_blank(),
        panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
        panel.grid.minor = element_blank())
print(pcoa_plot)
dev.off()

####################################################################

# Boxplot for Shannon index and richness for all 4 groups

####################################################################

# 1. Combine the long data frames for LBD vs Control and iRBD vs Control
all_shannon_richness <- rbind(lbd_vs_control_alpha_diversity_long_df, irbd_vs_control_alpha_diversity_long_df)

# 2. Rename the condition values to more descriptive names
all_shannon_richness_done <- all_shannon_richness %>%
  mutate(condition = case_when(
    condition == "lbd" ~ "LBD",
    condition == "irbd" ~ "iRBD",
    condition == "lbd_control" ~ "LBD_Control",
    condition == "irbd_control" ~ "iRBD_Control",
    TRUE ~ condition  # leave other values unchanged
  ))

# 3. Set the order of the condition factor levels for plotting
all_shannon_richness_done$condition <- factor(all_shannon_richness_done$condition, levels = c("LBD", "LBD_Control", "iRBD", "iRBD_Control"))

# 4. Create boxplots for Shannon index
all_shannon <- all_shannon_richness_done[all_shannon_richness_done$variable == "shannon", ]

pdf(file = "shannon_4_groups.pdf", height = 5.5, width = 3)
alpha_plot <- ggplot(all_shannon, 
                     aes(x = factor(condition), y = value, fill = factor(condition))) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 21, size = 1, position = position_jitter(width = 0.1)) +
  ylim(2.4, 5.4) +
  scale_fill_manual(values = c("#ff9274", "#55b7e6", "#fdc848", "#2DA248")) +
  theme_bw() + 
  theme(legend.title=element_blank(), 
        panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "", y = "Shannon Index") + 
  rremove('legend')
print(alpha_plot)
dev.off()

# 5. Create boxplots for species richness
all_richness <- all_shannon_richness_done[all_shannon_richness_done$variable == "richness", ]

pdf(file = "richness_4_groups.pdf", height = 5.5, width = 3.2)
alpha_plot <- ggplot(all_richness, 
                     aes(x = factor(condition), y = value, fill = factor(condition))) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape = 21, size = 1, position = position_jitter(width = 0.1)) +
  ylim(65, 330) +
  scale_fill_manual(values = c("#ff9274", "#55b7e6", "#2DA248", "#fdc848")) +
  theme_bw() + 
  theme(legend.title=element_blank(), 
        panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "", y = "Species Richness") + 
  rremove('legend')
print(alpha_plot)
dev.off()

####################################################################

# Distribution of bacterial families for LBD vs Control and iRBD vs Control

####################################################################

# 1. Prepare data for LBD vs Control and iRBD vs Control
lbd_vs_control_data <- species_all[species_all$condition %in% c("lbd", "lbd_control"), ]
irbd_vs_control_data <- species_all[species_all$condition %in% c("irbd", "irbd_control"), ]

# 2. Read taxonomy information
taxonomy <- read.csv(file = "../0-raw_data/taxonomy_info.csv")

# 3. LBD vs Control
# 3.1 Merge the taxonomy information with the species data
lbd_vs_control_t <- as.data.frame(t(lbd_vs_control_data[, 1:(ncol(lbd_vs_control_data)-33)]))
lbd_vs_control_tax <- merge(lbd_vs_control_t, taxonomy, by.x = "row.names", by.y = "species", all = F)

# 3.2 Select the relevant columns for family-level analysis
lbd_vs_control_family <- lbd_vs_control_tax[, c(1:51, 56)]
rownames(lbd_vs_control_family) <- lbd_vs_control_family$Row.names
lbd_vs_control_family <- lbd_vs_control_family[, -1]

# 3.3 Group by family and sum the abundance across samples
lbd_vs_control_family_new <- lbd_vs_control_family %>% group_by(family) %>% summarize(across(starts_with("BIOME"), sum, .names = "{.col}"))

# 3.4 Remove families with zero abundance across all samples
lbd_vs_control_family_done <- lbd_vs_control_family_new[rowSums(lbd_vs_control_family_new[, 2:51]) != 0, ]

# 3.5 Convert data to long format
df_long_lbd_vs_control <- lbd_vs_control_family_done %>%
  pivot_longer(cols = -family, names_to = "sample", values_to = "abundance")

# 3.6 For each sample, keep family with abundance >= 0.005, sum others
df_cleaned_lbd_vs_control <- df_long_lbd_vs_control %>%
  # Step 1: Mark families below threshold as "Others"
  mutate(family = if_else(abundance >= 0.05, family, "Others")) %>%
  # Step 2: Merge "Others" and *_unclassified into "OO"
  mutate(family = if_else(family == "Others" | grepl("_unclassified$", family), "OO", family)) %>%
  # Step 3: Sum abundance for each group in each sample
  group_by(sample, family) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# 3.7 Convert the cleaned data to wide format for further analysis
df_final_lbd_vs_control <- df_cleaned_lbd_vs_control %>%
  pivot_wider(names_from = sample, values_from = abundance, values_fill = 0)

# 3.8 Calculate the mean abundance for each family across all samples
df_final_lbd_vs_control$mean <- rowMeans(df_final_lbd_vs_control[, 2:51])

# 3.9 Merge the cleaned data with sample lookup information for plotting
sample_lookup_lbd_vs_control <- lbd_vs_control_data[, 1405:1406, drop = FALSE]
lbd_vs_control_family_combined <- merge(df_cleaned_lbd_vs_control, sample_lookup_lbd_vs_control, by.x = "sample", by.y = "row.names", all = F)

# 4. iRBD vs Control
# 4.1 Merge the taxonomy information with the species data
irbd_vs_control_t <- as.data.frame(t(irbd_vs_control_data[, 1:(ncol(irbd_vs_control_data)-33)]))
irbd_vs_control_tax <- merge(irbd_vs_control_t, taxonomy, by.x = "row.names", by.y = "species", all = F)

# 4.2 Select the relevant columns for family-level analysis
irbd_vs_control_family <- irbd_vs_control_tax[, c(1:21, 26)]
rownames(irbd_vs_control_family) <- irbd_vs_control_family$Row.names
irbd_vs_control_family <- irbd_vs_control_family[, -1]

# 4.3 Group by family and sum the abundance across samples
irbd_vs_control_family_new <- irbd_vs_control_family %>% group_by(family) %>% summarize(across(starts_with("BIOME"), sum, .names = "{.col}"))

# 4.4 Remove families with zero abundance across all samples
irbd_vs_control_family_done <- irbd_vs_control_family_new[rowSums(irbd_vs_control_family_new[, 2:21]) != 0, ]

# 4.5 Convert data to long format
df_long_irbd_vs_control <- irbd_vs_control_family_done %>%
  pivot_longer(cols = -family, names_to = "sample", values_to = "abundance")

# 4.6 For each sample, keep family with abundance >= 0.005, sum others
df_cleaned_irbd_vs_control <- df_long_irbd_vs_control %>%
  # Step 1: Mark families below threshold as "Others"
  mutate(family = if_else(abundance >= 0.05, family, "Others")) %>%
  # Step 2: Merge "Others" and *_unclassified into "OO"
  mutate(family = if_else(family == "Others" | grepl("_unclassified$", family), "OO", family)) %>%
  # Step 3: Sum abundance for each group in each sample
  group_by(sample, family) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

# 4.7 Convert the cleaned data to wide format for further analysis
df_final_irbd_vs_control <- df_cleaned_irbd_vs_control %>%
  pivot_wider(names_from = sample, values_from = abundance, values_fill = 0)

# 4.8 Calculate the mean abundance for each family across all samples
df_final_irbd_vs_control$mean <- rowMeans(df_final_irbd_vs_control[, 2:21])

# 4.9 Merge with sample lookup table  
sample_lookup_irbd_vs_control <- irbd_vs_control_data[, 1405:1406, drop = FALSE]
irbd_vs_control_family_combined <- merge(df_cleaned_irbd_vs_control, sample_lookup_irbd_vs_control, by.x = "sample", by.y = "row.names", all = F)

# 5. Combine the family-level data for LBD vs Control and iRBD vs Control
unique(lbd_vs_control_family_combined$family)
unique(irbd_vs_control_family_combined$family)
union(unique(lbd_vs_control_family_combined$family), 
      unique(irbd_vs_control_family_combined$family))

# 6. Combine the family-level data for all 4 groups and set the order of families for plotting
all_4_family <- rbind(lbd_vs_control_family_combined, irbd_vs_control_family_combined)
all_4_family$family <- factor(all_4_family$family, levels = c("FGB79294", "FGB3054", df_final_lbd_vs_control$family[order(df_final_lbd_vs_control$mean)][1:25], "Lachnospiraceae", "OO"))

# 7. Set the order of conditions for plotting
all_4_family$condition <- factor(all_4_family$condition, levels = c("lbd", "lbd_control", "irbd", "irbd_control"))

# 8. Plot stacked bar plot for family-level distribution across all 4 groups
# 8.1 Set the color palette for families
base_colors <- brewer.pal(12, "Set3")
my_colors <- colorRampPalette(base_colors)(28)
set.seed(2112)
my_colors <- sample(my_colors)

# 8.2 Create the stacked bar plot and save it as a PDF
pdf(file = "stacked_bar_plot_family_all.pdf", width = 20.5, height = 4)
ggplot(all_4_family, aes(x = new_sample_name, y = abundance, fill = family)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
  # scale_fill_brewer(palette = "Paired") +
  scale_fill_manual(values = c(my_colors, "lightblue")) +
  labs(x = "", y = "Relative Abundance") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.y = element_line(linetype = "dashed", color = "gray70", linewidth = 0.1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()

# 9. Create a wide-format data frame for all 4 groups and save it as a CSV file
all_4_family_wide_df <- pivot_wider(all_4_family[, 2:4], names_from = "new_sample_name", values_from = "abundance", values_fill = 0)
write.csv(all_4_family_wide_df, "distribution_family_relative_abundance.csv", row.names = F)