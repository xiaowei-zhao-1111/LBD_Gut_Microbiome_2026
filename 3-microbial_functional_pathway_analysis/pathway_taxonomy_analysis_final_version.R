library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(RColorBrewer)

###############################################################################
setwd("../3-microbial_functional_pathway_analysis")
###############################################################################

pathway_taxonomy <- read.delim("./pathway_analysis/pathway_taxonomy_new.tsv", header = TRUE, stringsAsFactors = FALSE)
rownames(pathway_taxonomy) <- pathway_taxonomy$X..Pathway
pathway_taxonomy <- pathway_taxonomy[, -1]
colnames(pathway_taxonomy) <- gsub("_S.*", "", colnames(pathway_taxonomy))

colSums(pathway_taxonomy)

pathway_taxonomy_t <- data.frame(t(pathway_taxonomy))
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

lbd_irbd_sig <- read.csv(file = "./pathway_analysis/lbd_vs_irbd/sig_lbd_vs_irbd_diff_abun.csv")

sig_pathways <- lbd_irbd_sig$pathways
sig_pathways_done <- sub(":.*", "", sig_pathways)

pathway_taxonomy_done_lbd_irbd <- pathway_taxonomy_done[pathway_taxonomy_done$pathway_ID %in% sig_pathways_done, ]

pathway_taxonomy_done_1 <- pathway_taxonomy_done_lbd_irbd[rowSums(pathway_taxonomy_done_lbd_irbd[, 1:70]) != 0, ]
pathway_taxonomy_done_2 <- pathway_taxonomy_done_1[, c(1:70)]
pathway_taxonomy_done_3 <- data.frame(t(pathway_taxonomy_done_2))
pathway_taxonomy_done_4 <- merge(pathway_taxonomy_done_3, metadata, by = "row.names", all = F)
pathway_taxonomy_done_5 <- pathway_taxonomy_done_4[pathway_taxonomy_done_4$condition %in% c("lbd", "irbd"), ]
pathway_taxonomy_done_5$new_name <- paste0(pathway_taxonomy_done_5$condition, "_", pathway_taxonomy_done_5$Row.names)
rownames(pathway_taxonomy_done_5) <- pathway_taxonomy_done_5$new_name
pathway_taxonomy_done_6 <- pathway_taxonomy_done_5[, c(2:429)]
pathway_taxonomy_done_7 <- as.data.frame(t(pathway_taxonomy_done_6))
pathway_taxonomy_done_7[] <- lapply(pathway_taxonomy_done_7, as.numeric)

pathway_taxonomy_done_7$mean_lbd <- apply(pathway_taxonomy_done_7[, grepl("lbd", colnames(pathway_taxonomy_done_7))], 1, function(x) mean(x))
pathway_taxonomy_done_7$mean_irbd <- apply(pathway_taxonomy_done_7[, grepl("irbd", colnames(pathway_taxonomy_done_7))], 1, function(x) mean(x))

pathway_taxonomy_done_7$median_lbd <- apply(pathway_taxonomy_done_7[, grepl("lbd", colnames(pathway_taxonomy_done_7))], 1, function(x) median(x))
pathway_taxonomy_done_7$median_irbd <- apply(pathway_taxonomy_done_7[, grepl("irbd", colnames(pathway_taxonomy_done_7))], 1, function(x) median(x))

## based on mean
pathway_taxonomy_df_mean <- pathway_taxonomy_done_7[pathway_taxonomy_done_7$mean_lbd != 0 | pathway_taxonomy_done_7$mean_irbd != 0, ]

pathway_taxonomy_df_mean$log2fc_mean <- log2((pathway_taxonomy_df_mean$mean_lbd+0.00001)/(pathway_taxonomy_df_mean$mean_irbd+0.00001))

lbd_col <- grep("lbd", colnames(pathway_taxonomy_df_mean[, 1:35]), value = TRUE)
irbd_col <- grep("irbd", colnames(pathway_taxonomy_df_mean[, 1:35]), value = TRUE)

pathway_taxonomy_df_mean$wx_p_value <- apply(pathway_taxonomy_df_mean[, 1:35], 1, function(x) {
  lbd_value <- as.numeric(x[lbd_col])
  irbd_value <- as.numeric(x[irbd_col])
  test <- wilcox.test(lbd_value, irbd_value)
  return(test$p.value)
})

pathway_taxonomy_df_mean$q_value <- p.adjust(pathway_taxonomy_df_mean$wx_p_value, method = "BH")
sig_log2fc_mean_df <- pathway_taxonomy_df_mean[pathway_taxonomy_df_mean$wx_p_value < 0.05 & abs(pathway_taxonomy_df_mean$log2fc_mean) > 1, 40:41, drop = FALSE]

sig_log2fc_mean_df_1 <- sig_log2fc_mean_df %>% 
  mutate(temp = row.names(sig_log2fc_mean_df)) %>% 
  separate(temp, into = c("pathway_ID", "taxonomy"), sep = "\\.g_|\\.II\\.u")

sig_log2fc_mean_df_2 <- sig_log2fc_mean_df_1[order(sig_log2fc_mean_df_1$log2fc_mean, decreasing = F), ]

sig_log2fc_mean_df_2$new <- rownames(sig_log2fc_mean_df_2)
all_pathways <- rownames(sig_log2fc_mean_df_2)

sig_log2fc_mean_df_2$new <- factor(sig_log2fc_mean_df_2$new, levels = all_pathways)

summary(sig_log2fc_mean_df_2$wx_p_value)

base_colors <- brewer.pal(3, "Set3")

pdf(file = "lbd_irbd_pathway_taxonomy_mean_new_v2.pdf", height = 5, width = 13)
alpha_plot <- ggplot(sig_log2fc_mean_df_2, 
                     aes(x = log2fc_mean, y = new, fill = factor(taxonomy), size = -log10(wx_p_value))) + 
  geom_point(shape = 21, color = "black") + 
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_fill_manual(values = base_colors) +
  scale_size_continuous(
    name = "P-value",
    breaks = c(-log10(0.015), -log10(0.01), -log10(0.005), -log10(0.001)),
    labels = c("P = 0.015", "P = 0.01", "P = 0.005", "P = 0.001")
  ) + 
  xlim(-4, 4) + 
  theme_bw() + 
  theme(legend.title=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "log2 fold-change in mean pathway relative abundance between LBD and iRBD", y = "")
print(alpha_plot)
dev.off()

   ###############################################################################

## based on median
pathway_taxonomy_df_median <- pathway_taxonomy_done_7[pathway_taxonomy_done_7$median_lbd != 0 | pathway_taxonomy_done_7$median_irbd != 0, ]

pathway_taxonomy_df_median$log2fc_median <- log2((pathway_taxonomy_df_median$median_lbd+0.00001)/(pathway_taxonomy_df_median$median_irbd+0.00001))

lbd_col <- grep("lbd", colnames(pathway_taxonomy_df_median[, 1:35]), value = TRUE)
irbd_col <- grep("irbd", colnames(pathway_taxonomy_df_median[, 1:35]), value = TRUE)

pathway_taxonomy_df_median$wx_p_value <- apply(pathway_taxonomy_df_median[, 1:35], 1, function(x) {
  lbd_value <- as.numeric(x[lbd_col])
  irbd_value <- as.numeric(x[irbd_col])
  test <- wilcox.test(lbd_value, irbd_value)
  return(test$p.value)
})

pathway_taxonomy_df_median$q_value <- p.adjust(pathway_taxonomy_df_median$wx_p_value, method = "BH")
sig_log2fc_median_df <- pathway_taxonomy_df_median[pathway_taxonomy_df_median$wx_p_value < 0.05, 40:41, drop = FALSE]

sig_log2fc_median_df_1 <- sig_log2fc_median_df %>% 
  mutate(temp = row.names(sig_log2fc_median_df)) %>% 
  separate(temp, into = c("pathway_ID", "taxonomy"), sep = "\\.g_|\\.II\\.u")

sig_log2fc_median_df_2 <- sig_log2fc_median_df_1[order(sig_log2fc_median_df_1$log2fc_median, decreasing = F), ]

sig_log2fc_median_df_2$new <- rownames(sig_log2fc_median_df_2)
all_pathways <- rownames(sig_log2fc_median_df_2)

sig_log2fc_median_df_2$new <- factor(sig_log2fc_median_df_2$new, levels = all_pathways)

unique(sig_log2fc_median_df_2$taxonomy)

base_colors <- brewer.pal(4, "Set3")

pdf(file = "lbd_irbd_pathway_taxonomy_median.pdf", height = 7, width = 13)
alpha_plot <- ggplot(sig_log2fc_median_df_2, 
                     aes(x = log2fc_median, y = new, fill = factor(taxonomy), size = -log10(wx_p_value))) + 
  geom_point(shape = 21, color = "black") + 
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
  scale_fill_manual(values = base_colors) +
  scale_size_continuous(
    name = "P-value",
    breaks = c(-log10(0.04), -log10(0.02), -log10(0.01), -log10(0.005), -log10(0.001)),
    labels = c("P = 0.04", "P = 0.02", "P = 0.01", "P = 0.005", "P = 0.001")
  ) + 
  xlim(-6, 6) + 
  theme_bw() + 
  theme(legend.title=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  labs(x = "log2 fold-change in median pathway relative abundance between LBD and iRBD", y = "")
print(alpha_plot)
dev.off()



