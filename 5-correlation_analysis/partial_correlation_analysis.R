library(ggplot2)
library(reshape)
library(tidyr)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ppcor)

###############################################################################
setwd("../5-correlation_analysis")
###############################################################################

# 1. Read the species and pathway relative abundance data and modify the sample names. 
species_data <- read.table(file = "./metaphlan_results_new.tsv", sep = "\t", header = T, row.names = 1)
species <- species_data[grepl("s__", rownames(species_data)) & !grepl("t__", rownames(species_data)), ]
rownames(species) <- gsub(".*s__", "", rownames(species))
colnames(species) <- gsub("metaphlan_|_S.*", "", names(species))

pathway_data <- read.csv(file = "./pathway_results.tsv", sep = "\t", header = T, row.names = 1)
colnames(pathway_data) <- gsub("_S.*", "", colnames(pathway_data))
pathway <- pathway_data

# 2. Read the metadata
metadata <- read.csv(file = "./0-raw_data/imputed_BMI_metadata.csv", sep = ",", header = T, row.names = 1)
species_clean <- species[, c(rownames(metadata))]
pathway_clean <- pathway[, c(rownames(metadata))]

# 3. Change the taxonomic data from percentage into proportion
colSums(species_clean)
species_prop <- data.frame(apply(species_clean, 2, function(x) x/sum(x)))
colSums(species_prop)

colSums(pathway_clean)
pathway_prop <- data.frame(apply(pathway_clean, 2, function(x) x/sum(x)))
colSums(pathway_prop)

# 4. Use prevalence cut-off (10^-4.5) to remove bacterial species and pathways that are present in very low abundance
cut_off_s <- 10^-4.7
species_cutoff <- species_prop
species_cutoff[species_cutoff <= cut_off_s] <- 0

cut_off_p <- 10^-4.5
pathway_cutoff <- pathway_prop
pathway_cutoff[pathway_cutoff < cut_off_p] <- 0

# 5. Add metadata to the proportion data
species_t <- as.data.frame(t(species_cutoff))
species_all <- merge(species_t, metadata, by = "row.names", all = F)
rownames(species_all) <- species_all$Row.names
species_done <- species_all[, -1]

pathway_t <- as.data.frame(t(pathway_cutoff))
pathway_all <- merge(pathway_t, metadata, by = "row.names", all = F)
rownames(pathway_all) <- pathway_all$Row.names
pathway_done <- pathway_all[, -1]

###############################################################################

# ── STAGE 1: fast screen via parametric pcor.test ─────────────────────────────
run_partial_cor_screen <- function(df, feature, clinical_measurement) {
  
  df1 <- df[, c(feature, "age", "BMI", clinical_measurement)]
  df2 <- na.omit(df1)
  if (nrow(df2) < 4) return(c(rho = NA, pval = NA))
  
  Z_full <- tryCatch(
    model.matrix(~ age + BMI, data = df2)[, -1, drop = FALSE],
    error = function(e) NULL)
  if (is.null(Z_full)) return(c(rho = NA, pval = NA))
  
  res <- tryCatch(
    pcor.test(x = df2[[1]], y = df2[[clinical_measurement]],
              z = Z_full, method = "spearman"),
    error = function(e) NULL)
  if (is.null(res) || is.na(res$estimate)) return(c(rho = NA, pval = NA))
  
  return(c(rho = res$estimate, pval = res$p.value))
}


# ── STAGE 2: resampling null + LOO stability (only on Stage 1 hits) ────────────
run_partial_cor_validate <- function(df, feature, clinical_measurement, n_perm = 1000) {
  
  df1 <- df[, c(feature, "age", "BMI", clinical_measurement)]
  df2 <- na.omit(df1)
  if (nrow(df2) < 4) return(c(pval_resample = NA, ci_lower = NA, ci_upper = NA, stability = NA))
  
  n_total <- nrow(df2)
  
  # real rho (recomputed cleanly on complete cases)
  Z_full <- tryCatch(
    model.matrix(~ age + BMI, data = df2)[, -1, drop = FALSE],
    error = function(e) NULL)
  if (is.null(Z_full)) return(c(pval_resample = NA, ci_lower = NA, ci_upper = NA, stability = NA))
  
  real_res <- tryCatch(
    pcor.test(x = df2[[1]], y = df2[[clinical_measurement]],
              z = Z_full, method = "spearman"),
    error = function(e) NULL)
  if (is.null(real_res) || is.na(real_res$estimate)) {
    return(c(pval_resample = NA, ci_lower = NA, ci_upper = NA, stability = NA))
  }
  real_rho <- real_res$estimate
  
  # resampling null distribution
  set.seed(2026)
  null_rho <- c()
  for (i in seq_len(n_perm)) {
    size   <- sample(ceiling(0.8 * n_total):(n_total - 1), 1)
    idx    <- sample(seq_len(n_total), size, replace = FALSE)
    df_sub <- df2[idx, ]
    
    if (var(df_sub[[1]], na.rm = TRUE) == 0) next
    
    df_sub[[clinical_measurement]] <- sample(df_sub[[clinical_measurement]])
    
    Z_sub <- tryCatch(
      model.matrix(~ age + BMI, data = df_sub)[, -1, drop = FALSE],
      error = function(e) NULL)
    if (is.null(Z_sub)) next
    
    res_null <- tryCatch(
      pcor.test(x = df_sub[[1]], y = df_sub[[clinical_measurement]],
                z = Z_sub, method = "spearman"),
      error = function(e) NULL)
    if (is.null(res_null) || is.na(res_null$estimate)) next
    
    null_rho <- c(null_rho, res_null$estimate)
  }
  
  pval_resample <- if (length(null_rho) > 0) mean(abs(null_rho) >= abs(real_rho)) else NA
  
  # LOO stability
  sig_count   <- 0
  total_count <- 0
  all_rho     <- c()
  
  for (iter in seq_len(n_total)) {
    df_sub <- df2[-iter, ]
    
    if (var(df_sub[[1]], na.rm = TRUE) == 0) next
    
    Z_sub <- tryCatch(
      model.matrix(~ age + BMI, data = df_sub)[, -1, drop = FALSE],
      error = function(e) NULL)
    if (is.null(Z_sub)) next
    
    res_loo <- tryCatch(
      pcor.test(x = df_sub[[1]], y = df_sub[[clinical_measurement]],
                z = Z_sub, method = "spearman"),
      error = function(e) NULL)
    if (is.null(res_loo) || is.na(res_loo$estimate)) next
    
    total_count <- total_count + 1
    all_rho     <- c(all_rho, res_loo$estimate)
    
    if (sign(res_loo$estimate) == sign(real_rho) && abs(res_loo$estimate) >= 0.4) {
      sig_count <- sig_count + 1
    }
    cat(sprintf("  Feature: %-40s | LOO iter: %d/%d | rho: %.3f\n",
                feature, iter, n_total, res_loo$estimate))
  }
  
  stability <- if (total_count > 0) sig_count / total_count else NA
  ci_lower  <- unname(quantile(all_rho, 0.025, na.rm = TRUE))
  ci_upper  <- unname(quantile(all_rho, 0.975, na.rm = TRUE))
  
  return(c(pval_resample = pval_resample,
           ci_lower      = ci_lower,
           ci_upper      = ci_upper,
           stability     = stability))
}


# ── WRAPPER: two-stage pipeline ────────────────────────────────────────────────
partial_cor_subsample <- function(data, data_type, clinical_measurement, group, threshold, n_perm = 1000) {
  
  df <- data %>% filter(condition == group)
  
  filtered_df <- df[, c(1:(ncol(df) - 32))][, colSums(df[, c(1:(ncol(df) - 32))] > 0) >= threshold * nrow(df)]
  
  features <- colnames(filtered_df)
  n_feat   <- length(features)
  
  # ── open log file ────────────────────────────────────────────────────────────
  log_file <- paste0("partial_cor_log_", clinical_measurement, "_", group, "_", data_type, ".txt")
  log_con  <- file(log_file, open = "wt")
  
  log <- function(...) {
    msg <- paste0(...)
    cat(msg, "\n")
    cat(msg, "\n", file = log_con)
  }
  
  on.exit(close(log_con), add = TRUE)
  
  log("=== STAGE 1: parametric screen across ", n_feat, " features ===")
  log("clinical_measurement : ", clinical_measurement)
  log("group                : ", group)
  log("data_type            : ", data_type)
  log("prevalence threshold : ", threshold)
  log("n_perm (Stage 2)     : ", n_perm)
  log("timestamp            : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  log(strrep("-", 60))
  
  # ── Stage 1: parametric screen ──────────────────────────────────────────────
  screen_rho  <- numeric(n_feat)
  screen_pval <- numeric(n_feat)
  
  for (i in seq_len(n_feat)) {
    log(sprintf("[%d/%d] %s", i, n_feat, features[i]))
    res            <- run_partial_cor_screen(df, features[i], clinical_measurement)
    screen_rho[i]  <- res[["rho"]]
    screen_pval[i] <- res[["pval"]]
  }
  
  screen_result <- data.frame(
    feature = features,
    rho     = screen_rho,
    pval    = screen_pval,
    qval    = p.adjust(screen_pval, method = "BH")
  )
  
  # save ALL screen results
  screen_result_order <- screen_result[order(screen_result$pval), ]
  write.csv(screen_result_order,
            file = paste0("partial_cor_screen_", clinical_measurement, "_", group, "_", data_type, ".csv"),
            row.names = FALSE)
  
  hits <- screen_result[!is.na(screen_result$pval) & screen_result$pval < 0.05, ]
  
  log(strrep("-", 60))
  log(sprintf("Stage 1 complete: %d / %d features pass (raw p < 0.05)", nrow(hits), n_feat))
  log(sprintf("Full screen results saved to: partial_cor_screen_%s_%s_%s.csv", clinical_measurement, group, data_type))
  log(strrep("-", 60))
  
  if (nrow(hits) == 0) {
    log("No features passed Stage 1. Returning screen results only.")
    return(screen_result_order)
  }
  
  # ── Stage 2: resampling + LOO on hits only ───────────────────────────────────
  log(sprintf("\n=== STAGE 2: resampling + LOO on %d hits ===", nrow(hits)))
  
  val_pval      <- numeric(nrow(hits))
  val_ci_lower  <- numeric(nrow(hits))
  val_ci_upper  <- numeric(nrow(hits))
  val_stability <- numeric(nrow(hits))
  
  for (i in seq_len(nrow(hits))) {
    feat <- hits$feature[i]
    log(sprintf("\n[%d/%d] %s", i, nrow(hits), feat))
    res              <- run_partial_cor_validate(df, feat, clinical_measurement, n_perm = n_perm)
    val_pval[i]      <- res[["pval_resample"]]
    val_ci_lower[i]  <- res[["ci_lower"]]
    val_ci_upper[i]  <- res[["ci_upper"]]
    val_stability[i] <- res[["stability"]]
  }
  
  result <- data.frame(
    feature         = hits$feature,
    rho             = hits$rho,
    pval_parametric = hits$pval,
    qval_parametric = hits$qval,
    pval_resample   = val_pval,
    qval_resample   = p.adjust(val_pval, method = "BH"),
    ci_lower        = val_ci_lower,
    ci_upper        = val_ci_upper,
    stability       = val_stability
  )
  
  result_order <- result[order(result$pval_resample), ]
  
  log(strrep("-", 60))
  log("Stage 2 complete. Results saved to:")
  log(sprintf("  Screen (all)      : partial_cor_screen_%s_%s_%s.csv", clinical_measurement, group, data_type))
  log(sprintf("  Validated (hits)  : partial_cor_res_%s_%s_%s.csv",    clinical_measurement, group, data_type))
  log(sprintf("  Log               : %s", log_file))
  log(sprintf("  Timestamp         : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  write.csv(result_order,
            file = paste0("partial_cor_res_", clinical_measurement, "_", group, "_", data_type, ".csv"),
            row.names = FALSE)
  
  return(result_order)
}


## LBD: species: CDR-SB
res_PC_cdr_species_lbd_new <- partial_cor_subsample(
  data = species_done,
  data_type = "species",
  clinical_measurement = "CDR_SB",
  group = "lbd",
  threshold = 0.5,
  n_perm = 1000)

## LBD: species: MOCA
res_PC_moca_species_lbd_new <- partial_cor_subsample(
  data = species_done,
  data_type = "species",
  clinical_measurement = "MOCA",
  group = "lbd",
  threshold = 0.5,
  n_perm = 1000)

## LBD: species: STMS
res_PC_stms_species_lbd_new <- partial_cor_subsample(
  data = species_done,
  data_type = "species",
  clinical_measurement = "STMS",
  group = "lbd",
  threshold = 0.5,
  n_perm = 1000)

## LBD: species: UPDRS3
res_PC_updrs_species_lbd_new <- partial_cor_subsample(
  data = species_done,
  data_type = "species",
  clinical_measurement = "UPDRS3",
  group = "lbd",
  threshold = 0.5,
  n_perm = 1000)

## LBD: pathway: CDR-SB
res_PC_cdr_pathway_lbd_new <- partial_cor_subsample(
  data = pathway_done,
  data_type = "pathway",
  clinical_measurement = "CDR_SB",
  group = "lbd",
  threshold = 0.5,
  n_perm = 1000)

## LBD: pathway: MOCA
res_PC_moca_pathway_lbd_new <- partial_cor_subsample(
  data = pathway_done,
  data_type = "pathway",
  clinical_measurement = "MOCA",
  group = "lbd",
  threshold = 0.5,
  n_perm = 1000)

## LBD: pathway: STMS
res_PC_stms_pathway_lbd_new <- partial_cor_subsample(
  data = pathway_done,
  data_type = "pathway",
  clinical_measurement = "STMS",
  group = "lbd",
  threshold = 0.5,
  n_perm = 1000)

## LBD: pathway: UPDRS
res_PC_updrs_pathway_lbd_new <- partial_cor_subsample(
  data = pathway_done,
  data_type = "pathway",
  clinical_measurement = "UPDRS3",
  group = "lbd",
  threshold = 0.5,
  n_perm = 1000)

## iRBD: species: CDR-SB
res_PC_cdr_species_irbd_new <- partial_cor_subsample(
  data = species_done,
  data_type = "species",
  clinical_measurement = "CDR_SB",
  group = "irbd",
  threshold = 0.5,
  n_perm = 1000)

## iRBD: species: MOCA
res_PC_moca_species_irbd_new <- partial_cor_subsample(
  data = species_done,
  data_type = "species",
  clinical_measurement = "MOCA",
  group = "irbd",
  threshold = 0.5,
  n_perm = 1000)

## iRBD: species: STMS
res_PC_stms_species_irbd_new <- partial_cor_subsample(
  data = species_done,
  data_type = "species",
  clinical_measurement = "STMS",
  group = "irbd",
  threshold = 0.5,
  n_perm = 1000)

## iRBD: species: UPDRS3
res_PC_updrs_species_irbd_new <- partial_cor_subsample(
  data = species_done,
  data_type = "species",
  clinical_measurement = "UPDRS3",
  group = "irbd",
  threshold = 0.5,
  n_perm = 1000)

## iRBD: pathway: CDR-SB
res_PC_cdr_pathway_irbd_new <- partial_cor_subsample(
  data = pathway_done,
  data_type = "pathway",
  clinical_measurement = "CDR_SB",
  group = "irbd",
  threshold = 0.5,
  n_perm = 1000)

## iRBD: pathway: MOCA
res_PC_moca_pathway_irbd_new <- partial_cor_subsample(
  data = pathway_done,
  data_type = "pathway",
  clinical_measurement = "MOCA",
  group = "irbd",
  threshold = 0.5,
  n_perm = 1000)

## iRBD: pathway: STMS
res_PC_stms_pathway_irbd_new <- partial_cor_subsample(
  data = pathway_done,
  data_type = "pathway",
  clinical_measurement = "STMS",
  group = "irbd",
  threshold = 0.5,
  n_perm = 1000)

## iRBD: pathway: UPDRS3
res_PC_updrs_pathway_irbd_new <- partial_cor_subsample(
  data = pathway_done,
  data_type = "pathway",
  clinical_measurement = "UPDRS3",
  group = "irbd",
  threshold = 0.5,
  n_perm = 1000)

## LBD group:

sig_cdr_sp_lbd <- res_PC_cdr_species_lbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "CDR-SB")
sig_moca_sp_lbd <- res_PC_moca_species_lbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "MoCA")
sig_stms_sp_lbd <- res_PC_stms_species_lbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "STMS")

lbd_species_all <- rbind(sig_cdr_sp_lbd, sig_moca_sp_lbd, sig_stms_sp_lbd)

lbd_species_all$feature <- sub("_", " ", lbd_species_all$feature)
lbd_species_all_ordered <- lbd_species_all[order(lbd_species_all$feature), ]

sig_cdr_pwy_lbd <- res_PC_cdr_pathway_lbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "CDR-SB")
sig_moca_pwy_lbd <- res_PC_moca_pathway_lbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "MoCA")
sig_stms_pwy_lbd <- res_PC_stms_pathway_lbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "STMS")

lbd_pathway_all <- rbind(sig_cdr_pwy_lbd, sig_moca_pwy_lbd, sig_stms_pwy_lbd)

lbd_pathway_all$feature <- sub(".*: ", "", lbd_pathway_all$feature)
lbd_pathway_all_ordered <- lbd_pathway_all[order(lbd_pathway_all$feature), ]

lbd_all <- rbind(lbd_species_all_ordered, lbd_pathway_all_ordered)

lbd_all$log10p <- -log10(lbd_all$pval_parametric)
lbd_all$group <- factor(lbd_all$group, levels = c("CDR-SB", "MoCA", "STMS"))
lbd_all$feature <- factor(lbd_all$feature, levels = rev(unique(lbd_all$feature)))

lbd_all$overall <- "lbd"

sig_cdr_sp_irbd <- res_PC_cdr_species_irbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "CDR-SB")
sig_moca_sp_irbd <- res_PC_moca_species_irbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "MoCA")
sig_stms_sp_irbd <- res_PC_stms_species_irbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "STMS")

irbd_species_all <- rbind(sig_cdr_sp_irbd, sig_moca_sp_irbd, sig_stms_sp_irbd)

irbd_species_all$feature <- sub("_", " ", irbd_species_all$feature)
irbd_species_all_ordered <- irbd_species_all[order(irbd_species_all$feature), ]

sig_cdr_pwy_irbd <- res_PC_cdr_pathway_irbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "CDR-SB")
sig_moca_pwy_irbd <- res_PC_moca_pathway_irbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "MoCA")
sig_stms_pwy_irbd <- res_PC_stms_pathway_irbd_new %>% filter(pval_resample < 0.05) %>% mutate(group = "STMS")

irbd_pathway_all <- rbind(sig_cdr_pwy_irbd, sig_moca_pwy_irbd, sig_stms_pwy_irbd)

irbd_pathway_all$feature <- sub(".*: ", "", irbd_pathway_all$feature)
irbd_pathway_all_ordered <- irbd_pathway_all[order(irbd_pathway_all$feature), ]

irbd_all <- rbind(irbd_species_all_ordered, irbd_pathway_all_ordered)

irbd_all$log10p <- -log10(irbd_all$pval_parametric)
irbd_all$group <- factor(irbd_all$group, levels = c("CDR-SB", "MoCA", "STMS"))
irbd_all$feature <- factor(irbd_all$feature, levels = rev(unique(irbd_all$feature)))

irbd_all$overall <- "irbd"


lbd_irbd_all <- rbind(lbd_all, irbd_all)
lbd_irbd_all$overall <- factor(lbd_irbd_all$overall, levels = c("lbd", "irbd"))

#lbd_irbd_all$overall_label <- factor(lbd_irbd_all$overall_label, levels = unique(lbd_irbd_all$overall_label))

summary(lbd_irbd_all$pval_parametric)

pdf(file = "lbd_irbd_correlation_v13.pdf", height = 5, width = 11)
ggplot(lbd_irbd_all, aes(x = group, y = feature)) +
  geom_point(aes(size = log10p, fill = rho), 
             shape = 21, color = "black", stroke = 1) +  # stroke controls border thickness
  scale_fill_gradient2(
    low = "#55b7e6", mid = "white", high = "#ff9274", midpoint = 0,
    name = "Spearman Rho", 
    limits = c(-1, 1)
  ) +
  scale_size_continuous(
    name = "P-value",
    breaks = c(-log10(0.0348), -log10(0.01), -log10(0.005), -log10(0.00058)),  # Define specific sizes to show in legend
    labels = c("P = 0.0348", "P = 0.01", "P = 0.005", "P = 0.00058")) + 
  facet_wrap(~ overall, scales = "free_y") + 
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.key.height = unit(0.6, "cm"),
    legend.key.width = unit(0.6, "cm")
  ) + 
  scale_y_discrete(position = "right")
dev.off()

# Save current RNG state
rng_state <- .Random.seed
saveRDS(rng_state, "rng_state.rds")

# Later, to restore it before re-running:
.Random.seed <- readRDS("rng_state.rds")