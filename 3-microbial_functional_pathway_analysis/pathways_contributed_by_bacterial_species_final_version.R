library(dplyr)
library(tidyr)
library(reshape)
library(ggplot2)
library(purrr)

###############################################################################
setwd("../3-microbial_functional_pathway_analysis")
###############################################################################

pathway_taxonomy <- read.delim("./pathway_analysis/pathway_taxonomy_new.tsv", header = TRUE, stringsAsFactors = FALSE)
rownames(pathway_taxonomy) <- pathway_taxonomy$X..Pathway
pathway_taxonomy <- pathway_taxonomy[, -1]
colnames(pathway_taxonomy) <- gsub("_S.*", "", colnames(pathway_taxonomy))

colSums(pathway_taxonomy)

pathway_taxonomy_t <- data.frame(t(pathway_taxonomy))
pathway_taxonomy_trans <- asin(sqrt(pathway_taxonomy_t))
metadata <- read.csv(file = "./metadata/imputed_BMI_metadata.csv", row.names = 1)
pathway_taxonomy_t_1 <- merge(pathway_taxonomy_t, metadata, by = "row.names")
rownames(pathway_taxonomy_t_1) <- pathway_taxonomy_t_1$Row.names
pathway_taxonomy_t_1 <- pathway_taxonomy_t_1[, -1]
# new_names <- paste(pathway_taxonomy_t_1$condition, rownames(pathway_taxonomy_t_1), sep = "_")
# rownames(pathway_taxonomy_t_1) <- new_names
pathway_taxonomy_t_2 <- pathway_taxonomy_t_1[, 1:17128]
pathway_taxonomy_t_3 <- data.frame(t(pathway_taxonomy_t_2))
rownames(pathway_taxonomy_t_3) <- rownames(pathway_taxonomy)

pathway_taxonomy_1 <- pathway_taxonomy_t_3 %>% 
  mutate(temp = row.names(pathway_taxonomy_t_3)) %>% 
  separate(temp, into = c("pathway_ID", "pathway_name"), sep = ": ")

pathway_taxonomy_done <- pathway_taxonomy_1 %>% 
  mutate(temp = pathway_name) %>% 
  separate(temp, into = c("pathway_name_done", "taxonomy"), sep = "\\|")

###############################################################################

lbd_vs_control_pathway_taxonomy <- pathway_taxonomy_done[, grepl("BIOME|pathway|taxonomy", colnames(pathway_taxonomy_done))]

## PWY-5030
lbd_vs_control_df <- lbd_vs_control_pathway_taxonomy[lbd_vs_control_pathway_taxonomy$pathway_ID == "PWY-5030", ]
lbd_vs_control_df_1 <- lbd_vs_control_df[rowSums(lbd_vs_control_df[, 1:50]) != 0, ]
lbd_vs_control_df_2 <- lbd_vs_control_df_1[, c(1:70, 74)]
rownames(lbd_vs_control_df_2) <- lbd_vs_control_df_2$taxonomy
lbd_vs_control_df_2$taxonomy <- NULL

lbd_vs_control_df_3 <- data.frame(t(lbd_vs_control_df_2))
lbd_vs_control_df_4 <- merge(lbd_vs_control_df_3, metadata, by = "row.names", all = F)

lbd_vs_control_df_5 <- lbd_vs_control_df_4[, c(2:27, 30)] %>% pivot_longer(cols = -new_sample_name, names_to = "taxonomy", values_to = "relab")

lbd_vs_control_df_6 <- lbd_vs_control_df_5 %>%
  mutate(condition = case_when(
    grepl("LBD-[0-9]+", new_sample_name) ~ "LBD", 
    grepl("iRBD-[0-9]+", new_sample_name) ~ "iRBD", 
    grepl("LBD-Control", new_sample_name) ~ "LBD_Control", 
    grepl("iRBD-Control", new_sample_name) ~ "iRBD_Control", 
    TRUE ~ new_sample_name  # leave other values unchanged
  ))

lbd_vs_control_df_6$condition <- factor(lbd_vs_control_df_6$condition, levels = c("LBD", "LBD_Control", "iRBD", "iRBD_Control"))

# Compute total per sample and set order within each condition
category_order <- lbd_vs_control_df_6 %>%
  group_by(condition, new_sample_name) %>%
  summarise(total = sum(relab), .groups = "drop") %>%
  arrange(condition, desc(total)) %>%
  mutate(order = row_number())

# Set the category as a factor based on the desired order
lbd_vs_control_df_6$new_sample_name <- factor(lbd_vs_control_df_6$new_sample_name, levels = category_order$new_sample_name)

library(RColorBrewer)
# Start with max Set3 colors
base_colors <- brewer.pal(12, "Set3")
# Expand using colorRampPalette
my_colors <- colorRampPalette(base_colors)(27)
set.seed(2112)
my_colors <- sample(my_colors)

pdf(file = "all_bacterial_contribution.pdf", width = 24, height = 5)
plot <- ggplot(lbd_vs_control_df_6, aes(x = new_sample_name, y = relab, fill = taxonomy)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(my_colors, "lightblue")) +
  labs(x = "", y = "Relative Abundance") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.y = element_line(linetype = "dashed", color = "gray70", linewidth = 0.1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank())
print(plot)
dev.off()

paired_test_df <- lbd_vs_control_df_4[, c(2:29, 31)]
paired_test_lbd_df_done <- paired_test_df[grepl("lbd", paired_test_df$condition), ]

# Select the columns only containing species abundances (excluding metadata columns like IDs and groups)

species_cols <- setdiff(names(paired_test_lbd_df_done), c("household_id", "condition"))

# Function to run Wilcoxon test for each species
result <- lapply(unique(paired_test_lbd_df_done$household_id), function(hh) {
  subset_df <- paired_test_lbd_df_done[paired_test_lbd_df_done$household_id == hh, ]
  if (nrow(subset_df) == 2) { # Only proceed if both case and control are present
    case_vector <- as.numeric(subset_df[subset_df$condition == "lbd", species_cols])
    control_vector <- as.numeric(subset_df[subset_df$condition == "lbd_control", species_cols])
    # Wilcoxon test: paired=TRUE
    test <- wilcox.test(case_vector, control_vector, paired = TRUE)
    return(data.frame(household_id = hh, p.value = test$p.value))
  } else {
    return(data.frame(household_id = hh, p.value = NA))
  }
})
wilcox_results_lbd <- do.call(rbind, result)

t1 <- paired_test_df_done[paired_test_df_done$household_id == "h21",]
wilcox.test(as.numeric(t1[1, 1:27]), as.numeric(t1[2, 1:27]), paired = T)

t2 <- paired_test_df_done[paired_test_df_done$household_id == "h25",]
wilcox.test(as.numeric(t2[1, 1:27]), as.numeric(t2[2, 1:27]), paired = T)

t3 <- paired_test_df_done[paired_test_df_done$household_id == "h7",]
wilcox.test(as.numeric(t3[1, 1:27]), as.numeric(t3[2, 1:27]), paired = T)


paired_test_irbd_df_done <- paired_test_df[grepl("irbd", paired_test_df$condition), ]

# Select the columns only containing species abundances (excluding metadata columns like IDs and groups)

species_cols <- setdiff(names(paired_test_irbd_df_done), c("household_id", "condition"))

# Function to run Wilcoxon test for each species
result <- lapply(unique(paired_test_irbd_df_done$household_id), function(hh) {
  subset_df <- paired_test_irbd_df_done[paired_test_irbd_df_done$household_id == hh, ]
  if (nrow(subset_df) == 2) { # Only proceed if both case and control are present
    case_vector <- as.numeric(subset_df[subset_df$condition == "irbd", species_cols])
    control_vector <- as.numeric(subset_df[subset_df$condition == "irbd_control", species_cols])
    # Wilcoxon test: paired=TRUE
    test <- wilcox.test(case_vector, control_vector, paired = TRUE)
    return(data.frame(household_id = hh, p.value = test$p.value))
  } else {
    return(data.frame(household_id = hh, p.value = NA))
  }
})
wilcox_results_irbd <- do.call(rbind, result)

t1 <- paired_test_irbd_df_done[paired_test_irbd_df_done$household_id == "h32",]
wilcox.test(as.numeric(t1[1, 1:27]), as.numeric(t1[2, 1:27]), paired = T)

t2 <- paired_test_df_done[paired_test_df_done$household_id == "h25",]
wilcox.test(as.numeric(t2[1, 1:27]), as.numeric(t2[2, 1:27]), paired = T)

t3 <- paired_test_df_done[paired_test_df_done$household_id == "h7",]
wilcox.test(as.numeric(t3[1, 1:27]), as.numeric(t3[2, 1:27]), paired = T)


paired_test_irbd_df_done$total_abundance <- rowSums(paired_test_irbd_df_done[, 1:27])
irbd_group <- paired_test_irbd_df_done$total_abundance[paired_test_irbd_df_done$condition == "irbd"]
irbd_control_group <- paired_test_irbd_df_done$total_abundance[paired_test_irbd_df_done$condition == "irbd_control"]

wilcox.test(lbd_group, lbd_control_group)
wilcox.test(irbd_group, irbd_control_group)
wilcox.test(lbd_group, irbd_group)

ks.test(lbd_group, lbd_control_group)
ks.test(irbd_group, irbd_control_group)
ks.test(lbd_group, irbd_group)

install.packages("kSamples")  # if you haven't installed it yet
library(kSamples)

# Suppose you already have LBD_total and Control_total as in your previous example
ad.test(lbd_group, lbd_control_group)
ad.test(irbd_group, irbd_control_group)
ad.test(lbd_group, irbd_group)

install.packages("cramer")  # or install.packages("goftest")
library(cramer)

# Assume LBD_total and Control_total are your total abundances for each sample in each group
cramer.test(lbd_group, lbd_control_group)
cramer.test(irbd_group, irbd_control_group)
cramer.test(lbd_group, irbd_group)


## PWY-6731
lbd_vs_control_df_6731 <- lbd_vs_control_pathway_taxonomy[lbd_vs_control_pathway_taxonomy$pathway_ID == "PWY-6731", ]

lbd_vs_control_df_6731_1 <- lbd_vs_control_df_6731[rowSums(lbd_vs_control_df_6731[, 1:50]) != 0, ]

lbd_vs_control_df_6731_2 <- lbd_vs_control_df_6731_1[, c(1:70, 74)]
rownames(lbd_vs_control_df_6731_2) <- lbd_vs_control_df_6731_2$taxonomy
lbd_vs_control_df_6731_2$taxonomy <- NULL

lbd_vs_control_df_6731_3 <- data.frame(t(lbd_vs_control_df_6731_2))
lbd_vs_control_df_6731_4 <- merge(lbd_vs_control_df_6731_3, metadata, by = "row.names", all = F)

lbd_vs_control_df_6731_5 <- lbd_vs_control_df_6731_4[, c(2:22, 24)] %>% pivot_longer(cols = -new_sample_name, names_to = "taxonomy", values_to = "relab")

lbd_vs_control_df_6731_6 <- lbd_vs_control_df_6731_5 %>%
  mutate(condition = case_when(
    grepl("LBD-[0-9]+", new_sample_name) ~ "LBD", 
    grepl("iRBD-[0-9]+", new_sample_name) ~ "iRBD", 
    grepl("LBD-Control", new_sample_name) ~ "LBD_Control", 
    grepl("iRBD-Control", new_sample_name) ~ "iRBD_Control", 
    TRUE ~ new_sample_name  # leave other values unchanged
  ))

lbd_vs_control_df_6731_6$condition <- factor(lbd_vs_control_df_6731_6$condition, levels = c("LBD", "LBD_Control", "iRBD", "iRBD_Control"))

# Compute total per sample and set order within each condition
category_order <- lbd_vs_control_df_6731_6 %>%
  group_by(condition, new_sample_name) %>%
  summarise(total = sum(relab), .groups = "drop") %>%
  arrange(condition, desc(total)) %>%
  mutate(order = row_number())

# Set the category as a factor based on the desired order
lbd_vs_control_df_6731_6$new_sample_name <- factor(lbd_vs_control_df_6731_6$new_sample_name, levels = category_order$new_sample_name)

library(RColorBrewer)
# Start with max Set3 colors
base_colors <- brewer.pal(12, "Set3")
# Expand using colorRampPalette
my_colors <- colorRampPalette(base_colors)(length(unique(lbd_vs_control_df_6731_6$taxonomy)))
set.seed(218)
my_colors <- sample(my_colors)

pdf(file = "bacterial_contribution_pwy_6731_new.pdf", width = 24, height = 4.5)
plot <- ggplot(lbd_vs_control_df_6731_6, aes(x = new_sample_name, y = relab, fill = taxonomy)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = my_colors) +
  labs(x = "", y = "Relative Abundance") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.y = element_line(linetype = "dashed", color = "gray70", linewidth = 0.1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank())
print(plot)
dev.off()

paired_test_df <- lbd_vs_control_df_6731_4[, c(2:23, 25)]

paired_test_lbd_df_done <- paired_test_df[grepl("lbd", paired_test_df$condition), ]

# Select the columns only containing species abundances (excluding metadata columns like IDs and groups)

species_cols <- setdiff(names(paired_test_lbd_df_done), c("household_id", "condition"))

# Function to run Wilcoxon test for each species
result <- lapply(unique(paired_test_lbd_df_done$household_id), function(hh) {
  subset_df <- paired_test_lbd_df_done[paired_test_lbd_df_done$household_id == hh, ]
  if (nrow(subset_df) == 2) { # Only proceed if both case and control are present
    case_vector <- as.numeric(subset_df[subset_df$condition == "lbd", species_cols])
    control_vector <- as.numeric(subset_df[subset_df$condition == "lbd_control", species_cols])
    # Wilcoxon test: paired=TRUE
    test <- wilcox.test(case_vector, control_vector, paired = TRUE)
    return(data.frame(household_id = hh, p.value = test$p.value))
  } else {
    return(data.frame(household_id = hh, p.value = NA))
  }
})
wilcox_results_lbd <- do.call(rbind, result)

t1 <- paired_test_lbd_df_done[paired_test_lbd_df_done$household_id == "h1",]
wilcox.test(as.numeric(t1[1, 1:21]), as.numeric(t1[2, 1:21]), paired = T)

paired_test_lbd_df_done$total_abundance <- rowSums(paired_test_lbd_df_done[, 1:21])

# Function to run Wilcoxon test for each species
result <- lapply(unique(paired_test_lbd_df_done$household_id), function(hh) {
  subset_df <- paired_test_lbd_df_done[paired_test_lbd_df_done$household_id == hh, ]
  if (nrow(subset_df) == 2) { # Only proceed if both case and control are present
    case_total <- subset_df[subset_df$condition == "lbd", "total_abundance"]
    control_total <- subset_df[subset_df$condition == "lbd_control", "total_abundance"]
    log2_fc <- log2((case_total+0.00001)/(control_total+0.00001))
    return(data.frame(household_id = hh, fold_change = log2_fc, control = control_total))
  } else {
    return(data.frame(household_id = hh, fold_change = NA, control = control_total))
  }
})
fc_lbd <- do.call(rbind, result)
fc_lbd_ordered <- fc_lbd[order(fc_lbd$fold_change, decreasing = T), ]
fc_lbd_ordered$household_id <- factor(fc_lbd_ordered$household_id, levels = fc_lbd_ordered$household_id)
fc_lbd_ordered$log10control <- -log10(fc_lbd_ordered$control)

# Plot
pdf(file = "log2fc_6731_household.pdf", height = 3, width = 8)


# Create the lollipop plot
ggplot(fc_lbd_ordered, aes(x = household_id, y = fold_change)) +
  geom_segment(aes(x = household_id, xend = household_id, y = 0, yend = fold_change), 
               color = "#b2b3a5", size = 1.5) +
  geom_point(aes(size = 1/log10control), shape = 21, fill = "#caaec5", color = "black") +
  scale_size_continuous(name = "Total bacteria contribution of pathway in Controls", 
                        breaks = c(1/-log10(0.0005), 1/-log10(0.001), 1/-log10(0.002), 1/-log10(0.003)),
                        labels = c("0.0005", "0.001", "0.002", "0.003")) +
  theme_bw() +
  labs(x = "Household ID", y = "Fold Change") +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()

result <- lapply(unique(paired_test_lbd_df_done$household_id), function(hh) {
  subset_df <- paired_test_lbd_df_done[paired_test_lbd_df_done$condition == cc, ]
  
  case_total <- subset_df[subset_df$condition == "lbd", "total_abundance"]
  log2_fc <- log2((case_total+0.00001)/(control_total+0.00001))
  return(data.frame(household_id = hh, fold_change = log2_fc))
  
})

lbd <- paired_test_lbd_df_done[paired_test_lbd_df_done$condition == "lbd", ]
lbd[lbd == 0] <- 0.00001
lbd_geo_mean <- data.frame(apply(lbd[, 1:21], 2, function(x) exp(mean(log(x)))))
colnames(lbd_geo_mean)[1] <- "geo_mean"
lbd_geo_mean$condition <- "LBD"
lbd_geo_mean$species <- rownames(lbd_geo_mean)

lbd_control <- paired_test_lbd_df_done[paired_test_lbd_df_done$condition == "lbd_control", ]
lbd_control[lbd_control == 0] <- 0.00001
lbd_control_geo_mean <- data.frame(apply(lbd_control[, 1:21], 2, function(x) exp(mean(log(x)))))
colnames(lbd_control_geo_mean)[1] <- "geo_mean"
lbd_control_geo_mean$condition <- "LBD_control"
lbd_control_geo_mean$species <- rownames(lbd_control_geo_mean)

graph_df <- rbind(lbd_geo_mean, lbd_control_geo_mean)


graph_df_order <- data.frame(species = colnames(paired_test_lbd_df_done)[1:21], 
                       LBD_geo_mean = lbd_geo_mean$geo_mean,
                       LBD_control_geo_mean = lbd_control_geo_mean$geo_mean)
graph_df_order <- graph_df_order[order(graph_df_order$LBD_geo_mean, decreasing = T), ]
graph_df$species <- factor(graph_df$species, levels = graph_df_order$species)

pdf(file = "geomean_species_6731.pdf", width = 7.5, height = 5.5)
plot <- ggplot(graph_df, 
               aes(x = species, y = geo_mean, fill = condition)) + 
  geom_col(position = "dodge", color = "black") +
  scale_fill_manual(values = c("#ff9274", "#55b7e6")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray70", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Geometric mean of pathway relative abundance contributed by different species")
print(plot)
dev.off()




lbd_group <- paired_test_lbd_df_done$total_abundance[paired_test_lbd_df_done$condition == "lbd"]
lbd_control_group <- paired_test_lbd_df_done$total_abundance[paired_test_lbd_df_done$condition == "lbd_control"]

# Suppose you already have LBD_total and Control_total as in your previous example
ad.test(lbd_group, lbd_control_group)
ad.test(irbd_group, irbd_control_group)
ad.test(lbd_group, irbd_group)

## PWY-7456
lbd_vs_control_df_7456 <- lbd_vs_control_pathway_taxonomy[lbd_vs_control_pathway_taxonomy$pathway_ID == "PWY-7456", ]

lbd_vs_control_df_7456_1 <- lbd_vs_control_df_7456[rowSums(lbd_vs_control_df_7456[, 1:50]) != 0, ]

lbd_vs_control_df_7456_2 <- lbd_vs_control_df_7456_1[, c(1:70, 74)]
rownames(lbd_vs_control_df_7456_2) <- lbd_vs_control_df_7456_2$taxonomy
lbd_vs_control_df_7456_2$taxonomy <- NULL

lbd_vs_control_df_7456_3 <- data.frame(t(lbd_vs_control_df_7456_2))
lbd_vs_control_df_7456_4 <- merge(lbd_vs_control_df_7456_3, metadata, by = "row.names", all = F)

lbd_vs_control_df_7456_5 <- lbd_vs_control_df_7456_4[, c(2:6, 8)] %>% pivot_longer(cols = -new_sample_name, names_to = "taxonomy", values_to = "relab")

lbd_vs_control_df_7456_6 <- lbd_vs_control_df_7456_5 %>%
  mutate(condition = case_when(
    grepl("LBD-[0-9]+", new_sample_name) ~ "LBD", 
    grepl("iRBD-[0-9]+", new_sample_name) ~ "iRBD", 
    grepl("LBD-Control", new_sample_name) ~ "LBD_Control", 
    grepl("iRBD-Control", new_sample_name) ~ "iRBD_Control", 
    TRUE ~ new_sample_name  # leave other values unchanged
  ))

lbd_vs_control_df_7456_6$condition <- factor(lbd_vs_control_df_7456_6$condition, levels = c("LBD", "LBD_Control", "iRBD", "iRBD_Control"))

# Compute total per sample and set order within each condition
category_order <- lbd_vs_control_df_7456_6 %>%
  group_by(condition, new_sample_name) %>%
  summarise(total = sum(relab), .groups = "drop") %>%
  arrange(condition, desc(total)) %>%
  mutate(order = row_number())

# Set the category as a factor based on the desired order
lbd_vs_control_df_7456_6$new_sample_name <- factor(lbd_vs_control_df_7456_6$new_sample_name, levels = category_order$new_sample_name)

library(RColorBrewer)
# Start with max Set3 colors
base_colors <- brewer.pal(5, "Set3")
# Expand using colorRampPalette
# my_colors <- colorRampPalette(base_colors)(length(unique(lbd_vs_control_df_7456_6$taxonomy)))
# set.seed(1111)
# my_colors <- sample(my_colors)

pdf(file = "bacterial_contribution_pwy_7456.pdf", width = 24, height = 5)
plot <- ggplot(lbd_vs_control_df_7456_6, aes(x = new_sample_name, y = relab, fill = taxonomy)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = base_colors) +
  labs(x = "", y = "Relative Abundance") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.y = element_line(linetype = "dashed", color = "gray70", linewidth = 0.1),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank())
print(plot)
dev.off()

paired_test_df <- lbd_vs_control_df_7456_4[, c(2:7, 9)]
paired_test_lbd_df_done <- paired_test_df[grepl("lbd", paired_test_df$condition), ]

# Select the columns only containing species abundances (excluding metadata columns like IDs and groups)

species_cols <- setdiff(names(paired_test_lbd_df_done), c("household_id", "condition"))

# Function to run Wilcoxon test for each species
result <- lapply(unique(paired_test_lbd_df_done$household_id), function(hh) {
  subset_df <- paired_test_lbd_df_done[paired_test_lbd_df_done$household_id == hh, ]
  if (nrow(subset_df) == 2) { # Only proceed if both case and control are present
    case_vector <- as.numeric(subset_df[subset_df$condition == "lbd", species_cols])
    control_vector <- as.numeric(subset_df[subset_df$condition == "lbd_control", species_cols])
    # Wilcoxon test: paired=TRUE
    test <- wilcox.test(case_vector, control_vector, paired = TRUE)
    return(data.frame(household_id = hh, p.value = test$p.value))
  } else {
    return(data.frame(household_id = hh, p.value = NA))
  }
})
wilcox_results_lbd <- do.call(rbind, result)

t1 <- paired_test_lbd_df_done[paired_test_lbd_df_done$household_id == "h12",]
wilcox.test(as.numeric(t1[1, 1:5]), as.numeric(t1[2, 1:5]), paired = T)

paired_test_lbd_df_done$total_abundance <- rowSums(paired_test_lbd_df_done[, 1:5])
lbd_group <- paired_test_lbd_df_done$total_abundance[paired_test_lbd_df_done$condition == "lbd"]
lbd_control_group <- paired_test_lbd_df_done$total_abundance[paired_test_lbd_df_done$condition == "lbd_control"]

# Suppose you already have LBD_total and Control_total as in your previous example
ad.test(lbd_group, lbd_control_group)
ad.test(irbd_group, irbd_control_group)
ad.test(lbd_group, irbd_group)

