library(lmerTest)
library(vegan)
library(ade4)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(RColorBrewer)

###############################################################################
folder_path <- "~/Desktop/LBD/analysis/alpha_beta_diversity"
setwd(folder_path)
###############################################################################

# 1. Read the species relative abundance data
species_data <- read.table(file = "~/Desktop/LBD/analysis/metaphlan_results_new.tsv", sep = "\t", header = T)
species <- species_data[grepl("s__", species_data$clade_name) & !grepl("t__", species_data$clade_name), ]

# 2. Modify taxonomy names
rownames(species) <- species$clade_name
new_rownames <- gsub(".*s__", "", rownames(species))
rownames(species) <- new_rownames
species <- species[, -1]

# 3. Modify sample names
col_names <- names(species)
new_colnames <- gsub("metaphlan_|_S.*", "", col_names)
colnames(species) <- new_colnames

# 4. Read the metadata
metadata <- read.csv(file = "~/Desktop/LBD/analysis/imputed_BMI_metadata.csv", sep = ",", header = T, row.names = 1)
species_clean <- species[, rownames(metadata)]

# 5. Change the taxonomic data from percentage into proportion
colSums(species_clean)
species_prop <- data.frame(apply(species_clean, 2, function(x) x/sum(x)))
colSums(species_prop)

# # 6. Find the appropriate presence cutoff for bacterial species data
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

# 7. Use prevalence cut-off to remove bacterial species that are present in very low abundance
## Number of low abundance species in each sample
cut_off <- 10^-4.7

## If use presence cutoff, how many low abundance species are needed to be removed
apply(species_prop, 2, function(x) sum(as.numeric(x) <= cut_off, na.rm = TRUE)) - apply(species_prop, 2, function(x) sum(as.numeric(x) <= 0, na.rm = TRUE))

species_cutoff <- species_prop
species_cutoff[species_cutoff <= cut_off] <- 0

## Number of absence species in each sample
apply(species_cutoff, 2, function(x) sum(as.numeric(x) == 0, na.rm = TRUE))

# 8. Add metadata to the transformed data
species_t <- as.data.frame(t(species_cutoff))
species_all <- merge(species_t, metadata, by = "row.names", all = F)
rownames(species_all) <- species_all$Row.names
species_all <- species_all[, -1]

###############################################################################

# LBD vs. Control

# 1. Create the folder if not exist
folder_name <- "lbd_vs_control_results"
folder <- file.path(folder_path, folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

lbd_vs_control_data <- species_all[species_all$condition %in% c("lbd", "lbd_control"), ]

# 2. Calculate the Shannon index and species richness for each group, then use a mixed-effects linear model to compute p-values.

shannon <- diversity(lbd_vs_control_data[, c(1:(ncol(lbd_vs_control_data)-33))], index = "shannon", MARGIN = 1)
simpson <- diversity(lbd_vs_control_data[, c(1:(ncol(lbd_vs_control_data)-33))], index = "simpson", MARGIN = 1)
invsimpson <- diversity(lbd_vs_control_data[, c(1:(ncol(lbd_vs_control_data)-33))], index = "invsimpson", MARGIN = 1)
richness <- apply(lbd_vs_control_data[, c(1:(ncol(lbd_vs_control_data)-33))], 1, function(x) sum(x>0))

lbd_vs_control_diversity <- data.frame("condition" = lbd_vs_control_data$condition, "household_id" = lbd_vs_control_data$household_id, "age" = lbd_vs_control_data$age, "BMI" = lbd_vs_control_data$BMI, "sex" = lbd_vs_control_data$sex, "richness" = richness, "shannon" = shannon, "simpson" = simpson, "invsimpson" = invsimpson)

write.csv(lbd_vs_control_diversity, paste(folder, "/alpha_diversity_results_lbd_vs_control.csv", sep = ""), row.names = T)

# lbd_vs_control_diversity <- read.csv(file = "~/Desktop/LBD/analysis/alpha_beta_diversity/lbd_vs_control_results/alpha_diversity_results_lbd_vs_control.csv", row.names = 1)

## Calculate p-values
shannon_p <- summary(lmer(shannon ~ age + condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

simpson_p <- summary(lmer(simpson ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

invsimpson_p <- summary(lmer(invsimpson ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

richness_p <- summary(lmer(richness ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

mixed_effect_p_lbd_vs_control <- data.frame(diversity = c("shannon", "simpson", "invsimpson", "richness"), value = c(shannon_p, simpson_p, invsimpson_p, richness_p))

write.csv(mixed_effect_p_lbd_vs_control, paste0(folder, "/mixed_effect_p_lbd_vs_control.csv"))

## long_df for boxplots
lbd_vs_control_alpha_diversity_long_df <- melt(lbd_vs_control_diversity, id.vars = "condition", measure.vars = c("shannon", "richness"))

## Option: p-value from Wilcoxon rank-sum test
shannon_wilcox <- with(lbd_vs_control_diversity, 
                       wilcox.test(shannon[condition == "lbd"], 
                                   shannon[condition == "lbd_control"], 
                                   paired = FALSE)$p.value)

simpson_wilcox <- with(lbd_vs_control_diversity, 
                       wilcox.test(simpson[condition == "lbd"], 
                                   simpson[condition == "lbd_control"], 
                                   paired = FALSE)$p.value)

invsimpson_wilcox <- with(lbd_vs_control_diversity, 
                          wilcox.test(invsimpson[condition == "lbd"], 
                                      invsimpson[condition == "lbd_control"], 
                                      paired = FALSE)$p.value)

richness_wilcox <- with(lbd_vs_control_diversity, 
                        wilcox.test(richness[condition == "lbd"], 
                                    richness[condition == "lbd_control"], 
                                    paired = FALSE)$p.value)

wilcoxon_p_lbd_vs_control <- data.frame(diversity = c("shannon", "simpson", "invsimpson", "richness"), value = c(shannon_wilcox, simpson_wilcox, invsimpson_wilcox, richness_wilcox))

write.csv(wilcoxon_p_lbd_vs_control, paste0(folder, "/wilcoxon_p_lbd_vs_control.csv"))

## Optional: Apply mixed-effects linear model to account for additional covariates.

shannon_p <- summary(lmer(shannon ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
shannon_p_age <- summary(lmer(shannon ~ condition + age + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
shannon_p_sex <- summary(lmer(shannon ~ condition + sex + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
shannon_p_BMI <- summary(lmer(shannon ~ condition + BMI + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
shannon_p_all <- summary(lmer(shannon ~ condition + age + sex + BMI + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
shannon_p_list <- c(shannon_p, shannon_p_age, shannon_p_sex, shannon_p_BMI, shannon_p_all)

simpson_p <- summary(lmer(simpson ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p_age <- summary(lmer(simpson ~ condition + age + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p_sex <- summary(lmer(simpson ~ condition + sex + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p_BMI <- summary(lmer(simpson ~ condition + BMI + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p_all <- summary(lmer(simpson ~ condition + age + sex + BMI + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p_list <- c(simpson_p, simpson_p_age, simpson_p_sex, simpson_p_BMI, simpson_p_all)

invsimpson_p <- summary(lmer(invsimpson ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p_age <- summary(lmer(invsimpson ~ condition + age + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p_sex <- summary(lmer(invsimpson ~ condition + sex + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p_BMI <- summary(lmer(invsimpson ~ condition + BMI + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p_all <- summary(lmer(invsimpson ~ condition + age + sex + BMI + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p_list <- c(invsimpson_p, invsimpson_p_age, invsimpson_p_sex, invsimpson_p_BMI, invsimpson_p_all)

richness_p <- summary(lmer(richness ~ condition + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p_age <- summary(lmer(richness ~ condition + age + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p_sex <- summary(lmer(richness ~ condition + sex + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p_BMI <- summary(lmer(richness ~ condition + BMI + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p_all <- summary(lmer(richness ~ condition + age + sex + BMI + (1|household_id), data = lbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p_list <- c(richness_p, richness_p_age, richness_p_sex, richness_p_BMI, richness_p_all)

mixed_effects_all <- data.frame("model" = c("condition_only", "condition_age", "condition_sex", "condition_BMI", "condition_all"), shannon_p_list, simpson_p_list, invsimpson_p_list, richness_p_list) 

write.csv(mixed_effects_all, paste0(folder, "/mixed_effect_p_all_lbd_vs_control.csv"))

# 3. Create a Principle Coordinate Analysis plot

## Use arcsine square root transformation for relative abundance data
lbd_vs_control_transformed <- asin(sqrt(lbd_vs_control_data[, 1:(ncol(lbd_vs_control_data)-33)]))

lbd_vs_control_permanova <- merge(lbd_vs_control_transformed, metadata, by = "row.names", all = F)
rownames(lbd_vs_control_permanova) <- lbd_vs_control_permanova$Row.names
lbd_vs_control_permanova <- lbd_vs_control_permanova[, -1]

## 
lbd <- lbd_vs_control_permanova[lbd_vs_control_permanova$condition == "lbd", ]
summary(lbd$age)

## Test for association between age, sex, BMI, and disease condition
wilcox.test(age~condition, lbd_vs_control_permanova)

print(table(lbd_vs_control_permanova$sex, lbd_vs_control_permanova$condition))
fisher.test(table(lbd_vs_control_permanova$sex, lbd_vs_control_permanova$condition))

summary(lmer(BMI ~ condition + (1|household_id), data = lbd_vs_control_permanova, REML = F))[["coefficients"]]

## Perform PERMANOVA

# Set random seed to ensure reproducibility of permutation-based p-values
set.seed(10)

# Perform PERMANOVA using marginal testing (`by = "margin"`) to assess the individual contribution of each variable 
# (age, sex, BMI, condition) to microbiome beta diversity, adjusting for all other variables.
# Bray-Curtis dissimilarity is calculated from arcsine square root transformed relative abundance data.
# Stratification by household_id accounts for paired or clustered samples (e.g., matched case-control design).
set.seed(10)
results_all_condition <- adonis2(
  lbd_vs_control_permanova[, 1:(ncol(lbd_vs_control_permanova)-33)] ~ age + sex + BMI + condition,
  data = lbd_vs_control_permanova,
  permutations = 999,
  strata = lbd_vs_control_permanova$household_id,
  by = "margin"
)
# Save PERMANOVA output to text file
capture.output(results_all_condition, file = paste0(folder, "/PERMANOVA_results_all_marginal_lbd_vs_control.txt"))

# Repeatability ensured by setting the same seed before each individual model
set.seed(10)
# Test effect of condition only (univariate model)
results_condition_only <- adonis2(
  lbd_vs_control_permanova[, 1:(ncol(lbd_vs_control_permanova)-33)] ~ condition,
  data = lbd_vs_control_permanova,
  permutations = 999,
  strata = lbd_vs_control_permanova$household_id
)
capture.output(results_condition_only, file = paste0(folder, "/PERMANOVA_results_condition_only_lbd_vs_control.txt"))

set.seed(10)
# Test effect of age only
results_age_only <- adonis2(
  lbd_vs_control_permanova[, 1:(ncol(lbd_vs_control_permanova)-33)] ~ age,
  data = lbd_vs_control_permanova,
  permutations = 999,
  strata = lbd_vs_control_permanova$household_id
)
capture.output(results_age_only, file = paste0(folder, "/PERMANOVA_results_age_only_lbd_vs_control.txt"))

set.seed(10)
# Test effect of BMI only
results_BMI_only <- adonis2(
  lbd_vs_control_permanova[, 1:(ncol(lbd_vs_control_permanova)-33)] ~ BMI,
  data = lbd_vs_control_permanova,
  permutations = 999,
  strata = lbd_vs_control_permanova$household_id
)
capture.output(results_BMI_only, file = paste0(folder, "/PERMANOVA_results_BMI_only_lbd_vs_control.txt"))

set.seed(10)
# Test effect of sex only
results_sex_only <- adonis2(
  lbd_vs_control_permanova[, 1:(ncol(lbd_vs_control_permanova)-33)] ~ sex,
  data = lbd_vs_control_permanova,
  permutations = 999,
  strata = lbd_vs_control_permanova$household_id
)
capture.output(results_sex_only, file = paste0(folder, "/PERMANOVA_results_sex_only_lbd_vs_control.txt"))

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
pdf(file = paste(folder, "/pcoa_lbd_vs_control_v2.pdf", sep = ""), width = 7.5, height = 4.5)
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

###############################################################################

# iRBD vs. Control

# 1. Create the folder if not exist
folder_name <- "irbd_vs_control_results"
folder <- file.path(folder_path, folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

irbd_vs_control_data <- species_all[species_all$condition %in% c("irbd", "irbd_control"), ]

# 2. Calculate the Shannon index and species richness for each group, then use a mixed-effects linear model to compute p-values.

shannon_diversity <- diversity(irbd_vs_control_data[, c(1:(ncol(irbd_vs_control_data)-33))], index = "shannon", MARGIN = 1)
simpson_diversity <- diversity(irbd_vs_control_data[, c(1:(ncol(irbd_vs_control_data)-33))], index = "simpson", MARGIN = 1)
invsimpson_diversity <- diversity(irbd_vs_control_data[, c(1:(ncol(irbd_vs_control_data)-33))], index = "invsimpson", MARGIN = 1)
richness <- apply(irbd_vs_control_data[, c(1:(ncol(irbd_vs_control_data)-33))], 1, function(x) sum(x>0))

irbd_vs_control_diversity <- data.frame("condition" = irbd_vs_control_data$condition, "household_id" = irbd_vs_control_data$household_id, "age" = irbd_vs_control_data$age, "BMI" = irbd_vs_control_data$BMI, "sex" = irbd_vs_control_data$sex, "richness" = richness, "shannon" = shannon_diversity, "simpson" = simpson_diversity, "invsimpson" = invsimpson_diversity)

write.csv(irbd_vs_control_diversity, paste(folder, "/alpha_diversity_results_irbd_vs_control.csv", sep = ""), row.names = T)

# irbd_vs_control_diversity <- read.csv(file = "~/Desktop/LBD/analysis/alpha_beta_diversity/irbd_vs_control_results/alpha_diversity_results_irbd_vs_control.csv", row.names = 1)

## Calculate p-values
shannon_p <- summary(lmer(shannon ~ age + condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

simpson_p <- summary(lmer(simpson ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

invsimpson_p <- summary(lmer(invsimpson ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

richness_p <- summary(lmer(richness ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]

mixed_effect_p_irbd_vs_control <- data.frame(diversity = c("shannon", "simpson", "invsimpson", "richness"), value = c(shannon_p, simpson_p, invsimpson_p, richness_p))

write.csv(mixed_effect_p_irbd_vs_control, paste0(folder, "/mixed_effect_p_irbd_vs_control.csv"))

## long_df for boxplots
irbd_vs_control_alpha_diversity_long_df <- melt(irbd_vs_control_diversity, id.vars = "condition", measure.vars = c("shannon", "richness"))

## Option: p-value from Wilcoxon rank-sum test
shannon_wilcox <- with(irbd_vs_control_diversity, 
                       wilcox.test(shannon[condition == "irbd"], 
                                   shannon[condition == "irbd_control"], 
                                   paired = FALSE)$p.value)

simpson_wilcox <- with(irbd_vs_control_diversity, 
                       wilcox.test(simpson[condition == "irbd"], 
                                   simpson[condition == "irbd_control"], 
                                   paired = FALSE)$p.value)

invsimpson_wilcox <- with(irbd_vs_control_diversity, 
                          wilcox.test(invsimpson[condition == "irbd"], 
                                      invsimpson[condition == "irbd_control"], 
                                      paired = FALSE)$p.value)

richness_wilcox <- with(irbd_vs_control_diversity, 
                        wilcox.test(richness[condition == "irbd"], 
                                    richness[condition == "irbd_control"], 
                                    paired = FALSE)$p.value)

wilcoxon_p_irbd_vs_control <- data.frame(diversity = c("shannon", "simpson", "invsimpson", "richness"), value = c(shannon_wilcox, simpson_wilcox, invsimpson_wilcox, richness_wilcox))

write.csv(wilcoxon_p_irbd_vs_control, paste0(folder, "/wilcoxon_p_irbd_vs_control.csv"))

## Optional: Apply mixed-effects linear model to account for additional covariates.

shannon_p <- summary(lmer(shannon ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
shannon_p_age <- summary(lmer(shannon ~ condition + age + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
shannon_p_sex <- summary(lmer(shannon ~ condition + sex + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
shannon_p_BMI <- summary(lmer(shannon ~ condition + BMI + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
shannon_p_all <- summary(lmer(shannon ~ condition + age + sex + BMI + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
shannon_p_list <- c(shannon_p, shannon_p_age, shannon_p_sex, shannon_p_BMI, shannon_p_all)

simpson_p <- summary(lmer(simpson ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p_age <- summary(lmer(simpson ~ condition + age + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p_sex <- summary(lmer(simpson ~ condition + sex + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p_BMI <- summary(lmer(simpson ~ condition + BMI + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p_all <- summary(lmer(simpson ~ condition + age + sex + BMI + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
simpson_p_list <- c(simpson_p, simpson_p_age, simpson_p_sex, simpson_p_BMI, simpson_p_all)

invsimpson_p <- summary(lmer(invsimpson ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p_age <- summary(lmer(invsimpson ~ condition + age + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p_sex <- summary(lmer(invsimpson ~ condition + sex + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p_BMI <- summary(lmer(invsimpson ~ condition + BMI + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p_all <- summary(lmer(invsimpson ~ condition + age + sex + BMI + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
invsimpson_p_list <- c(invsimpson_p, invsimpson_p_age, invsimpson_p_sex, invsimpson_p_BMI, invsimpson_p_all)

richness_p <- summary(lmer(richness ~ condition + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p_age <- summary(lmer(richness ~ condition + age + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p_sex <- summary(lmer(richness ~ condition + sex + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p_BMI <- summary(lmer(richness ~ condition + BMI + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p_all <- summary(lmer(richness ~ condition + age + sex + BMI + (1|household_id), data = irbd_vs_control_diversity, REML = F))[["coefficients"]][2,5]
richness_p_list <- c(richness_p, richness_p_age, richness_p_sex, richness_p_BMI, richness_p_all)

mixed_effects_all <- data.frame("model" = c("condition_only", "condition_age", "condition_sex", "condition_BMI", "condition_all"), shannon_p_list, simpson_p_list, invsimpson_p_list, richness_p_list) 

write.csv(mixed_effects_all, paste0(folder, "/mixed_effect_p_all_irbd_vs_control.csv"))

# 3. Create a Principle Coordinate Analysis plot

## Use arcsine square root transformation for relative abundance data
irbd_vs_control_transformed <- asin(sqrt(irbd_vs_control_data[, 1:(ncol(irbd_vs_control_data)-33)]))
irbd_vs_control_permanova <- merge(irbd_vs_control_transformed, metadata, by = "row.names", all = F)
rownames(irbd_vs_control_permanova) <- irbd_vs_control_permanova$Row.names
irbd_vs_control_permanova <- irbd_vs_control_permanova[, -1]

irbd <- irbd_vs_control_permanova[irbd_vs_control_permanova$condition == "irbd", ]
summary(irbd$age)

## Test for association between age, sex, BMI, and disease condition
wilcox.test(age~condition, irbd_vs_control_permanova)

print(table(irbd_vs_control_permanova$sex, irbd_vs_control_permanova$condition))
fisher.test(table(irbd_vs_control_permanova$sex, irbd_vs_control_permanova$condition))

summary(lmer(BMI ~ condition + (1|household_id), data = irbd_vs_control_permanova, REML = F))[["coefficients"]]

## Perform PERMANOVA
set.seed(10)
adonis2(
  irbd_vs_control_permanova[, 1:(ncol(irbd_vs_control_permanova)-33)] ~ BMI,
  data = irbd_vs_control_permanova,
  permutations = 999,
  strata = irbd_vs_control_permanova$household_id
)



## PERMANOVA results for condition only, with household_id as random effect
set.seed(10) ## Since the permutation calculation uses random seeds, setting the seed is necessary to ensure reproducible results.
results_condition_only <- adonis2(irbd_vs_control_permanova[, 1:(ncol(irbd_vs_control_permanova)-32)] ~ condition, data = irbd_vs_control_permanova, permutations = 999, strata = irbd_vs_control_permanova$household_id, by = "terms")
capture.output(results_condition_only, file = paste0(folder, "/PERMANOVA_results_condition_only_irbd_vs_control.txt"))

## PERMANOVA results for condition and age, with household_id as random effect
set.seed(10)
results_age_condition <- adonis2(irbd_vs_control_permanova[, 1:(ncol(irbd_vs_control_permanova)-32)] ~ age + condition, data = irbd_vs_control_permanova, permutations = 999, strata = irbd_vs_control_permanova$household_id, by = "terms")
capture.output(results_age_condition, file = paste0(folder, "/PERMANOVA_results_age_condition_irbd_vs_control.txt"))

## PERMANOVA results for condition and BMI, with household_id as random effect
set.seed(10)
results_BMI_condition <- adonis2(irbd_vs_control_permanova[, 1:(ncol(irbd_vs_control_permanova)-32)] ~ BMI + condition, data = irbd_vs_control_permanova, permutations = 999, strata = irbd_vs_control_permanova$household_id, by = "terms")
capture.output(results_BMI_condition, file = paste0(folder, "/PERMANOVA_results_BMI_condition_irbd_vs_control.txt"))

## PERMANOVA results for condition, age, and sex and BMI, with household_id as random effect
set.seed(10)
results_age_BMI_condition <- adonis2(irbd_vs_control_permanova[, 1:(ncol(irbd_vs_control_permanova)-32)] ~ age + BMI + condition, data = irbd_vs_control_permanova, permutations = 999, strata = irbd_vs_control_permanova$household_id, by = "terms")
capture.output(results_age_BMI_condition, file = paste0(folder, "/PERMANOVA_results_age_BMI_condition_irbd_vs_control.txt"))

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
pdf(file = paste(folder, "/pcoa_irbd_vs_control_v2.pdf", sep = ""), width = 7.5, height = 4.5)
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

###############################################################################

# LBD vs. iRBD

# 1. Create the folder if not exist
folder_name <- "lbd_vs_irbd_results"
folder <- file.path(folder_path, folder_name)

# Check if the folder already exists
if (!file.exists(folder)) {
  # If it doesn't exist, create the folder
  dir.create(folder)
}

lbd_vs_irbd_data <- species_all[species_all$condition %in% c("lbd", "irbd"), ]

# 2. Calculate the Shannon index and species richness for each group, then use a normal linear model to compute p-values.

shannon <- diversity(lbd_vs_irbd_data[, c(1:(ncol(lbd_vs_irbd_data)-32))], index = "shannon", MARGIN = 1)
simpson <- diversity(lbd_vs_irbd_data[, c(1:(ncol(lbd_vs_irbd_data)-32))], index = "simpson", MARGIN = 1)
invsimpson <- diversity(lbd_vs_irbd_data[, c(1:(ncol(lbd_vs_irbd_data)-32))], index = "invsimpson", MARGIN = 1)
richness <- apply(lbd_vs_irbd_data[, c(1:(ncol(lbd_vs_irbd_data)-32))], 1, function(x) sum(x>0))

lbd_vs_irbd_diversity <- data.frame("condition" = lbd_vs_irbd_data$condition, "household_id" = lbd_vs_irbd_data$household_id, "age" = lbd_vs_irbd_data$age, "BMI" = lbd_vs_irbd_data$BMI, "sex" = lbd_vs_irbd_data$sex, "richness" = richness, "shannon" = shannon, "simpson" = simpson, "invsimpson" = invsimpson)

write.csv(lbd_vs_irbd_diversity, paste(folder, "/alpha_diversity_results_lbd_vs_irbd.csv", sep = ""), row.names = T)

# lbd_vs_irbd_diversity <- read.csv(file = "~/Desktop/LBD/analysis/alpha_beta_diversity/lbd_vs_irbd_results/alpha_diversity_results_lbd_vs_irbd.csv", row.names = 1)

## Calculate p-values
shannon_p <- summary(lm(shannon ~ condition, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]

simpson_p <- summary(lm(simpson ~ condition, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]

invsimpson_p <- summary(lm(invsimpson ~ condition, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]

richness_p <- summary(lm(richness ~ condition, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]

lm_p_lbd_vs_control <- data.frame(diversity = c("shannon", "simpson", "invsimpson", "richness"), value = c(shannon_p, simpson_p, invsimpson_p, richness_p))

write.csv(lm_p_lbd_vs_control, paste0(folder, "/lm_p_lbd_vs_control.csv"))

## long_df for boxplots
lbd_vs_irbd_alpha_diversity_long_df <- melt(lbd_vs_irbd_diversity, id.vars = "condition", measure.vars = c("shannon", "richness"))

## Option: p-value from Wilcoxon rank-sum test
shannon_wilcox <- with(lbd_vs_irbd_diversity, 
                       wilcox.test(shannon[condition == "lbd"], 
                                   shannon[condition == "irbd"], 
                                   paired = FALSE)$p.value)

simpson_wilcox <- with(lbd_vs_irbd_diversity, 
                       wilcox.test(simpson[condition == "lbd"], 
                                   simpson[condition == "irbd"], 
                                   paired = FALSE)$p.value)

invsimpson_wilcox <- with(lbd_vs_irbd_diversity, 
                          wilcox.test(invsimpson[condition == "lbd"], 
                                      invsimpson[condition == "irbd"], 
                                      paired = FALSE)$p.value)

richness_wilcox <- with(lbd_vs_irbd_diversity, 
                        wilcox.test(richness[condition == "lbd"], 
                                    richness[condition == "irbd"], 
                                    paired = FALSE)$p.value)

wilcoxon_p_lbd_vs_irbd <- data.frame(diversity = c("shannon", "simpson", "invsimpson", "richness"), value = c(shannon_wilcox, simpson_wilcox, invsimpson_wilcox, richness_wilcox))

write.csv(wilcoxon_p_lbd_vs_irbd, paste0(folder, "/wilcoxon_p_lbd_vs_irbd.csv"))

## Optional: Apply mixed-effects linear model to account for additional covariates.

shannon_p <- summary(lm(shannon ~ condition, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
shannon_p_age <- summary(lm(shannon ~ condition + age, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
shannon_p_sex <- summary(lm(shannon ~ condition + sex, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
shannon_p_BMI <- summary(lm(shannon ~ condition + BMI, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
shannon_p_all <- summary(lm(shannon ~ condition + age + sex + BMI, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
shannon_p_list <- c(shannon_p, shannon_p_age, shannon_p_sex, shannon_p_BMI, shannon_p_all)

simpson_p <- summary(lm(simpson ~ condition, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
simpson_p_age <- summary(lm(simpson ~ condition + age, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
simpson_p_sex <- summary(lm(simpson ~ condition + sex, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
simpson_p_BMI <- summary(lm(simpson ~ condition + BMI, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
simpson_p_all <- summary(lm(simpson ~ condition + age + sex + BMI, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
simpson_p_list <- c(simpson_p, simpson_p_age, simpson_p_sex, simpson_p_BMI, simpson_p_all)

invsimpson_p <- summary(lm(invsimpson ~ condition, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
invsimpson_p_age <- summary(lm(invsimpson ~ condition + age, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
invsimpson_p_sex <- summary(lm(invsimpson ~ condition + sex, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
invsimpson_p_BMI <- summary(lm(invsimpson ~ condition + BMI, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
invsimpson_p_all <- summary(lm(invsimpson ~ condition + age + sex + BMI, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
invsimpson_p_list <- c(invsimpson_p, invsimpson_p_age, invsimpson_p_sex, invsimpson_p_BMI, invsimpson_p_all)

richness_p <- summary(lm(richness ~ condition, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
richness_p_age <- summary(lm(richness ~ condition + age, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
richness_p_sex <- summary(lm(richness ~ condition + sex, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
richness_p_BMI <- summary(lm(richness ~ condition + BMI, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
richness_p_all <- summary(lm(richness ~ condition + age + sex + BMI, data = lbd_vs_irbd_diversity))[["coefficients"]][2,4]
richness_p_list <- c(richness_p, richness_p_age, richness_p_sex, richness_p_BMI, richness_p_all)

lm_all <- data.frame("model" = c("condition_only", "condition_age", "condition_sex", "condition_BMI", "condition_all"), shannon_p_list, simpson_p_list, invsimpson_p_list, richness_p_list) 

write.csv(lm_all, paste0(folder, "/lm_p_all_lbd_vs_irbd.csv"))

# 3. Create a Principle Coordinate Analysis plot

## Use arcsine square root transformation for relative abundance data
lbd_vs_irbd_transformed <- asin(sqrt(lbd_vs_irbd_data[, 1:(ncol(lbd_vs_irbd_data)-32)]))

lbd_vs_irbd_permanova <- merge(lbd_vs_irbd_transformed, metadata, by = "row.names", all = F)
rownames(lbd_vs_irbd_permanova) <- lbd_vs_irbd_permanova$Row.names
lbd_vs_irbd_permanova <- lbd_vs_irbd_permanova[, -1]

## Test for association between age, sex, BMI, and disease condition
wilcox.test(age~condition, lbd_vs_irbd_permanova)

print(table(lbd_vs_irbd_permanova$sex, lbd_vs_irbd_permanova$condition))
fisher.test(table(lbd_vs_irbd_permanova$sex, lbd_vs_irbd_permanova$condition))

wilcox.test(BMI~condition, lbd_vs_irbd_permanova)

## Perform PERMANOVA
## PERMANOVA results for condition only, with household_id as random effect
set.seed(10) ## Since the permutation calculation uses random seeds, setting the seed is necessary to ensure reproducible results.
results_condition_only <- adonis2(lbd_vs_irbd_permanova[, 1:(ncol(lbd_vs_irbd_permanova)-32)] ~ condition, data = lbd_vs_irbd_permanova, permutations = 999, by = "terms")
capture.output(results_condition_only, file = paste0(folder, "/PERMANOVA_results_condition_only_lbd_vs_irbd.txt"))

## PERMANOVA results for condition and age, with household_id as random effect
set.seed(10)
results_age_condition <- adonis2(lbd_vs_irbd_permanova[, 1:(ncol(lbd_vs_irbd_permanova)-32)] ~ age + condition, data = lbd_vs_irbd_permanova, permutations = 999, by = "terms")
capture.output(results_age_condition, file = paste0(folder, "/PERMANOVA_results_age_condition_lbd_vs_irbd.txt"))

## PERMANOVA results for condition and sex, with household_id as random effect
set.seed(10)
results_sex_condition <- adonis2(lbd_vs_irbd_permanova[, 1:(ncol(lbd_vs_irbd_permanova)-32)] ~ sex + condition, data = lbd_vs_irbd_permanova, permutations = 999, by = "terms")
capture.output(results_sex_condition, file = paste0(folder, "/PERMANOVA_results_sex_condition_lbd_vs_irbd.txt"))

## PERMANOVA results for condition and BMI, with household_id as random effect
set.seed(10)
results_BMI_condition <- adonis2(lbd_vs_irbd_permanova[, 1:(ncol(lbd_vs_irbd_permanova)-32)] ~ BMI + condition, data = lbd_vs_irbd_permanova, permutations = 999, by = "terms")
capture.output(results_BMI_condition, file = paste0(folder, "/PERMANOVA_results_BMI_condition_lbd_vs_irbd.txt"))

## PERMANOVA results for condition, age, and sex and BMI, with household_id as random effect
set.seed(10)
results_age_sex_BMI_condition <- adonis2(lbd_vs_irbd_permanova[, 1:(ncol(lbd_vs_irbd_permanova)-32)] ~ age + sex + BMI + condition, data = lbd_vs_irbd_permanova, permutations = 999, by = "terms")
capture.output(results_age_sex_BMI_condition, file = paste0(folder, "/PERMANOVA_results_age_sex_BMI_condition_lbd_vs_irbd.txt"))

## Calculate distance matrix
distance_matrix <- vegdist(lbd_vs_irbd_permanova[, 1:(ncol(lbd_vs_irbd_permanova)-32)], method = "bray")

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
pc_plot_data$Group <- as.factor(lbd_vs_irbd_permanova$condition)
pc_plot_data$Group <- factor(pc_plot_data$Group, levels = c("lbd", "irbd"))

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
pdf(file = paste(folder, "/pcoa_lbd_vs_irbd.pdf", sep = ""), width = 5.65, height = 4.5)
pcoa_plot <- ggplot(new_pc_pot, aes(x = c_pc1, y = c_pc2, color = Group)) + 
  scale_colour_manual(values = c("#ff9274", "#2DA248")) + 
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

###############################################################################

## Boxplot for Shannon index and richness

all_shannon_richness <- rbind(lbd_vs_control_alpha_diversity_long_df, irbd_vs_control_alpha_diversity_long_df)

all_shannon_richness_done <- all_shannon_richness %>%
  mutate(condition = case_when(
    condition == "lbd" ~ "LBD",
    condition == "irbd" ~ "iRBD",
    condition == "lbd_control" ~ "LBD_Control",
    condition == "irbd_control" ~ "iRBD_Control",
    TRUE ~ condition  # leave other values unchanged
  ))

all_shannon_richness_done$condition <- factor(all_shannon_richness_done$condition, levels = c("LBD", "LBD_Control", "iRBD", "iRBD_Control"))

all_shannon <- all_shannon_richness_done[all_shannon_richness_done$variable == "shannon", ]

pdf(file = paste0("~/Desktop/LBD/analysis/alpha_beta_diversity/shannon_4_groups_solid_clean.pdf"), height = 5.5, width = 3)
alpha_plot <- ggplot(all_shannon, 
                     aes(x = factor(condition), y = value, fill = factor(condition))) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white",   # outline color
               outlier.size = 2,         # point size
               outlier.shape = 21) +
  geom_point(shape = 21, size = 1, position = position_jitter(width = 0.1)) +
  ylim(2.4, 5.4) +
  scale_fill_manual(values = c("#ff9274", "#55b7e6", "#2DA248", "#fdc848")) +
  theme_bw() + 
  theme(legend.title=element_blank(), 
        panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "", y = "Shannon Index") + 
  rremove('legend')
print(alpha_plot)
dev.off()

all_richness <- all_shannon_richness_done[all_shannon_richness_done$variable == "richness", ]

pdf(file = paste0("~/Desktop/LBD/analysis/alpha_beta_diversity/richness_4_groups_solid_clean.pdf"), height = 5.5, width = 3.2)
alpha_plot <- ggplot(all_richness, 
                     aes(x = factor(condition), y = value, fill = factor(condition))) + 
  geom_boxplot(outlier.color = "black", 
               outlier.fill = "white",   # outline color
               outlier.size = 2,         # point size
               outlier.shape = 21) +
  # geom_point(shape = 21, size = 1, position = position_jitter(width = 0.1)) +
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

###############################################################################

## stacked bar plots to demonstrate relative abundance of phylum and family groups of each samples

lbd_vs_control_data <- species_all[species_all$condition %in% c("lbd", "lbd_control"), ]
irbd_vs_control_data <- species_all[species_all$condition %in% c("irbd", "irbd_control"), ]

taxonomy <- read.csv(file = "~/Desktop/LBD/analysis/taxonomy_info.csv")
lbd_vs_control_t <- as.data.frame(t(lbd_vs_control_data[, 1:(ncol(lbd_vs_control_data)-33)]))
lbd_vs_control_tax <- merge(lbd_vs_control_t, taxonomy, by.x = "row.names", by.y = "species", all = F)

###############################################################################

lbd_vs_control_family <- lbd_vs_control_tax[, c(1:51, 56)]
rownames(lbd_vs_control_family) <- lbd_vs_control_family$Row.names
lbd_vs_control_family <- lbd_vs_control_family[, -1]
lbd_vs_control_family_new <- lbd_vs_control_family %>% group_by(family) %>% summarize(across(starts_with("BIOME"), sum, .names = "{.col}"))
lbd_vs_control_family_done <- lbd_vs_control_family_new[rowSums(lbd_vs_control_family_new[, 2:51]) != 0, ]

# Step 1: Convert data to long format
df_long <- lbd_vs_control_family_done %>%
  pivot_longer(cols = -family, names_to = "sample", values_to = "abundance")

# Step 2: For each sample, keep family with abundance >= 0.005, sum others
df_cleaned <- df_long %>%
  # Step 1: Mark families below threshold as "Others"
  mutate(family = if_else(abundance >= 0.05, family, "Others")) %>%
  # Step 2: Merge "Others" and *_unclassified into "OO"
  mutate(family = if_else(family == "Others" | grepl("_unclassified$", family), "OO", family)) %>%
  # Step 3: Sum abundance for each group in each sample
  group_by(sample, family) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

df_final_lbd_vs_control <- df_cleaned %>%
  pivot_wider(names_from = sample, values_from = abundance, values_fill = 0)

df_final_lbd_vs_control$mean <- rowMeans(df_final_lbd_vs_control[, 2:51])

sample_lookup <- lbd_vs_control_data[, 1405:1406, drop = FALSE]
lbd_vs_control_family_combined <- merge(df_cleaned, sample_lookup, by.x = "sample", by.y = "row.names", all = F)

# Start with max Set3 colors
base_colors <- brewer.pal(12, "Set3")
# Expand using colorRampPalette
my_colors <- colorRampPalette(base_colors)(27)
set.seed(2134)
my_colors <- sample(my_colors)

pdf(file = "stacked_bar_plot_family.pdf", width = 12, height = 4)
ggplot(lbd_vs_control_family_combined, aes(x = new_sample_name, y = abundance, fill = family)) +
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

###############################################################################

irbd_vs_control_t <- as.data.frame(t(irbd_vs_control_data[, 1:(ncol(irbd_vs_control_data)-33)]))
irbd_vs_control_tax <- merge(irbd_vs_control_t, taxonomy, by.x = "row.names", by.y = "species", all = F)

###############################################################################

irbd_vs_control_family <- irbd_vs_control_tax[, c(1:21, 26)]
rownames(irbd_vs_control_family) <- irbd_vs_control_family$Row.names
irbd_vs_control_family <- irbd_vs_control_family[, -1]
irbd_vs_control_family_new <- irbd_vs_control_family %>% group_by(family) %>% summarize(across(starts_with("BIOME"), sum, .names = "{.col}"))
irbd_vs_control_family_done <- irbd_vs_control_family_new[rowSums(irbd_vs_control_family_new[, 2:21]) != 0, ]

# Step 1: Convert data to long format
df_long <- irbd_vs_control_family_done %>%
  pivot_longer(cols = -family, names_to = "sample", values_to = "abundance")

# Step 2: For each sample, keep family with abundance >= 0.005, sum others
df_cleaned <- df_long %>%
  # Step 1: Mark families below threshold as "Others"
  mutate(family = if_else(abundance >= 0.05, family, "Others")) %>%
  # Step 2: Merge "Others" and *_unclassified into "OO"
  mutate(family = if_else(family == "Others" | grepl("_unclassified$", family), "OO", family)) %>%
  # Step 3: Sum abundance for each group in each sample
  group_by(sample, family) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

df_final <- df_cleaned %>%
  pivot_wider(names_from = sample, values_from = abundance, values_fill = 0)

df_final$mean <- rowMeans(df_final[, 2:21])

sample_lookup <- irbd_vs_control_data[, 1405:1406, drop = FALSE]
irbd_vs_control_family_combined <- merge(df_cleaned, sample_lookup, by.x = "sample", by.y = "row.names", all = F)

# irbd_vs_control_family_combined$family <- factor(irbd_vs_control_family_combined$family, levels = c(df_final$family[order(df_final$mean)]))
# 
# # Step 1: Get Lachnospiraceae abundance per sample
# Lachnospiraceae_order <- irbd_vs_control_family_combined %>%
#   filter(family == "Lachnospiraceae") %>%
#   arrange(desc(abundance)) %>%
#   pull(new_sample_name)
# 
# # Step 2: Reorder sample factor levels by decreasing Firmicutes abundance
# irbd_vs_control_family_combined$new_sample_name <- factor(irbd_vs_control_family_combined$new_sample_name, levels = Lachnospiraceae_order)

# Start with max Set3 colors
base_colors <- brewer.pal(12, "Set3")
# Expand using colorRampPalette
my_colors <- colorRampPalette(base_colors)(16)
set.seed(1111)
my_colors <- sample(my_colors)

pdf(file = "stacked_bar_plot_family_irbd_vs_control.pdf", width = 6, height = 4)
ggplot(irbd_vs_control_family_combined, aes(x = sample, y = abundance, fill = family)) +
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

###############################################################################

unique(lbd_vs_control_family_combined$family)
unique(irbd_vs_control_family_combined$family)
union(unique(lbd_vs_control_family_combined$family), 
      unique(irbd_vs_control_family_combined$family))

all_4_family <- rbind(lbd_vs_control_family_combined, irbd_vs_control_family_combined)
all_4_family$family <- factor(all_4_family$family, levels = c("FGB79294", "FGB3054", df_final_lbd_vs_control$family[order(df_final_lbd_vs_control$mean)][1:25], "Lachnospiraceae", "OO"))

# # Step 1: Get Firmicutes abundance per sample
# lachnospiraceae_order <- all_4_family %>%
#   filter(family == "Lachnospiraceae") %>%
#   arrange(desc(abundance)) %>%
#   pull(new_sample_name)
# 
# # Step 2: Reorder sample factor levels by decreasing Firmicutes abundance
# all_4_family$new_sample_name <- factor(all_4_family$new_sample_name, levels = lachnospiraceae_order)

all_4_family$condition <- factor(all_4_family$condition, levels = c("lbd", "lbd_control", "irbd", "irbd_control"))


all_4_family_wide_df <- pivot_wider(all_4_family[, 2:4], names_from = "new_sample_name", values_from = "abundance", values_fill = 0)
write.csv(all_4_family_wide_df, paste0(folder_path, "/all_relative_abundance.csv"), row.names = F)

sum_row <- colSums(all_4_family_wide_df[all_4_family_wide_df$family %in% c("Bacteroidaceae",
                                                                           "Lachnospiraceae",
                                                                           "Oscillospiraceae"), -1])

all_4_family_wide_df <- rbind(all_4_family_wide_df, c("Sum_Bact_Lach_Oscillo", sum_row))


Lachnospiraceae_LBD <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Lachnospiraceae", grepl("^LBD-\\d+$", colnames(all_4_family_wide_df))])
summary(Lachnospiraceae_LBD)

Lachnospiraceae_LBD_control <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Lachnospiraceae", grepl("LBD-Control", colnames(all_4_family_wide_df))])
summary(Lachnospiraceae_LBD_control)

Lachnospiraceae_iRBD <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Lachnospiraceae", grepl("^iRBD-\\d+$", colnames(all_4_family_wide_df))])
summary(Lachnospiraceae_iRBD)

Lachnospiraceae_iRBD_control <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Lachnospiraceae", grepl("iRBD-Control", colnames(all_4_family_wide_df))])
summary(Lachnospiraceae_iRBD_control)

OO <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "OO", 2:ncol(all_4_family_wide_df)])
summary(OO)

Bacteroidaceae <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Bacteroidaceae", 2:ncol(all_4_family_wide_df)])
summary(Bacteroidaceae)

Bacteroidaceae_LBD <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Bacteroidaceae", grepl("^LBD-\\d+$", colnames(all_4_family_wide_df))])
summary(Bacteroidaceae_LBD)

Bacteroidaceae_LBD_control <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Bacteroidaceae", grepl("LBD-Control", colnames(all_4_family_wide_df))])
summary(Bacteroidaceae_LBD_control)

Bacteroidaceae_iRBD <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Bacteroidaceae", grepl("^iRBD-\\d+$", colnames(all_4_family_wide_df))])
summary(Bacteroidaceae_iRBD)

Bacteroidaceae_iRBD_control <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Bacteroidaceae", grepl("iRBD-Control", colnames(all_4_family_wide_df))])
summary(Bacteroidaceae_iRBD_control)

Oscillospiraceae

Oscillospiraceae_LBD <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Oscillospiraceae", grepl("^LBD-\\d+$", colnames(all_4_family_wide_df))])
summary(Oscillospiraceae_LBD)

Oscillospiraceae_LBD_control <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Oscillospiraceae", grepl("LBD-Control", colnames(all_4_family_wide_df))])
summary(Oscillospiraceae_LBD_control)

Oscillospiraceae_iRBD <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Oscillospiraceae", grepl("^iRBD-\\d+$", colnames(all_4_family_wide_df))])
summary(Oscillospiraceae_iRBD)

Oscillospiraceae_iRBD_control <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Oscillospiraceae", grepl("iRBD-Control", colnames(all_4_family_wide_df))])
summary(Oscillospiraceae_iRBD_control)

Rikenellaceae

Rikenellaceae <- as.numeric(all_4_family_wide_df[all_4_family_wide_df$family == "Rikenellaceae", 2:ncol(all_4_family_wide_df)])
summary(Rikenellaceae)

# Start with max Set3 colors
base_colors <- brewer.pal(12, "Set3")
# Expand using colorRampPalette
my_colors <- colorRampPalette(base_colors)(28)
set.seed(2112)
my_colors <- sample(my_colors)

pdf(file = "~/Desktop/LBD/analysis/alpha_beta_diversity/stacked_bar_plot_family_all_new_NEW.pdf", width = 20.5, height = 4)
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
