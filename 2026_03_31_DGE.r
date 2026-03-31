suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(stringr)
  library(tibble)
})

indir  <- "/scratch/user/uqwsuwak/sc/pseudobulk_lm_vs_nonlm_by_patient_AMP"
outdir <- file.path(indir, "DESeq2_results_AMP")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

min_patients_per_group <- 2
min_genes_total_count  <- 10

safe_name <- function(x) {
  x <- gsub("[ /()+\\-]+", "_", x)
  x <- gsub("[^A-Za-z0-9_.]", "", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

make_volcano <- function(res_df, title_txt, outfile) {
  plot_df <- res_df %>%
    mutate(
      neglog10padj = ifelse(is.na(padj) | padj <= 0, NA, -log10(padj)),
      sig = ifelse(!is.na(padj) & padj < 0.05, "padj < 0.05", "ns"),
      label = ifelse(!is.na(padj) & padj < 0.05, gene, NA)
    )
  
  p <- ggplot(plot_df, aes(x = log2FoldChange, y = neglog10padj)) +
    geom_point(aes(color = sig), alpha = 0.8, size = 1.8, na.rm = TRUE) +
    geom_text_repel(
      data = subset(plot_df, !is.na(label)),
      aes(label = label),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.alpha = 0.5,
      na.rm = TRUE
    ) +
    scale_color_manual(values = c("ns" = "grey70", "padj < 0.05" = "#D55E00")) +
    labs(
      title = title_txt,
      x = "log2FoldChange (LM vs nonLM)",
      y = "-log10(padj)",
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )
  
  ggsave(outfile, p, width = 8, height = 6, dpi = 300)
}

count_files <- list.files(indir, pattern = "^counts__.*\\.csv$", full.names = TRUE)

if (length(count_files) == 0) {
  stop("No counts files found in: ", indir)
}

run_summary <- list()

for (cf in count_files) {
  
  base_name <- basename(cf)
  meta_file <- file.path(indir, str_replace(base_name, "^counts__", "metadata__"))
  
  if (!file.exists(meta_file)) {
    message("Skipping: metadata file missing for ", base_name)
    next
  }
  
  message("\nProcessing: ", base_name)
  
  counts <- read.csv(cf, row.names = 1, check.names = FALSE)
  meta   <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  
  required_meta_cols <- c("patient_id", "pathotype2")
  missing_cols <- setdiff(required_meta_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    message("Skipping: metadata missing columns for ", base_name, " -> ", paste(missing_cols, collapse = ", "))
    next
  }
  
  # parse location and cell type from filename
  # counts__ST__CD8_GZMK.csv -> ST / CD8_GZMK
  m <- str_match(base_name, "^counts__([^_]+)__(.*)\\.csv$")
  if (is.na(m[1, 2]) || is.na(m[1, 3])) {
    message("Skipping: could not parse filename ", base_name)
    next
  }
  
  location_val <- m[1, 2]
  celltype_val <- m[1, 3]
  
  # clean metadata columns
  if (!"location" %in% colnames(meta)) {
    meta$location <- location_val
  }
  if (!"ct_combined2" %in% colnames(meta)) {
    meta$ct_combined2 <- celltype_val
  }
  
  meta$location <- trimws(meta$location)
  meta$ct_combined2 <- trimws(meta$ct_combined2)
  
  # replace blank / NA with filename-derived values
  meta$location[is.na(meta$location) | meta$location == ""] <- location_val
  meta$ct_combined2[is.na(meta$ct_combined2) | meta$ct_combined2 == ""] <- celltype_val
  
  # one row per patient
  meta <- meta %>% distinct(patient_id, .keep_all = TRUE)
  
  if (!all(colnames(counts) %in% meta$patient_id)) {
    message("Skipping: some count columns not found in metadata for ", base_name)
    next
  }
  
  meta <- meta[match(colnames(counts), meta$patient_id), , drop = FALSE]
  
  if (!all(meta$patient_id == colnames(counts))) {
    message("Skipping: patient_id order mismatch for ", base_name)
    next
  }
  
  keep <- meta$pathotype2 %in% c("LM", "nonLM")
  meta <- meta[keep, , drop = FALSE]
  counts <- counts[, keep, drop = FALSE]
  
  if (ncol(counts) == 0 || nrow(counts) == 0) {
    message("Skipping: empty counts after filtering for ", celltype_val, " / ", location_val)
    next
  }
  
  group_table <- table(meta$pathotype2)
  
  if (!all(c("LM", "nonLM") %in% names(group_table))) {
    message("Skipping: one group missing for ", celltype_val, " / ", location_val)
    next
  }
  
  if (group_table["LM"] < min_patients_per_group || group_table["nonLM"] < min_patients_per_group) {
    message("Skipping: not enough patients per group for ", celltype_val, " / ", location_val)
    next
  }
  
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"
  
  keep_genes <- rowSums(counts) >= min_genes_total_count
  counts_filt <- counts[keep_genes, , drop = FALSE]
  
  if (nrow(counts_filt) == 0) {
    message("Skipping: no genes left after filtering for ", celltype_val, " / ", location_val)
    next
  }
  
  meta$pathotype2 <- factor(meta$pathotype2, levels = c("nonLM", "LM"))
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_filt,
    colData   = meta,
    design    = ~ pathotype2
  )
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("pathotype2", "LM", "nonLM"))
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
  
  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene")
  
  prefix <- paste0(safe_name(location_val), "__", safe_name(celltype_val))
  
  write.csv(
    res_df,
    file = file.path(outdir, paste0("DEG_results__", prefix, ".csv")),
    row.names = FALSE
  )
  
  write.csv(
    norm_counts,
    file = file.path(outdir, paste0("normalized_counts__", prefix, ".csv")),
    row.names = FALSE
  )
  
  make_volcano(
    res_df    = res_df,
    title_txt = paste0(location_val, " | ", celltype_val, " | LM vs nonLM"),
    outfile   = file.path(outdir, paste0("volcano__", prefix, ".png"))
  )
  
  run_summary[[length(run_summary) + 1]] <- data.frame(
    location = location_val,
    cell_type = celltype_val,
    n_patients = ncol(counts_filt),
    n_genes_tested = nrow(counts_filt),
    LM_patients = as.integer(group_table["LM"]),
    nonLM_patients = as.integer(group_table["nonLM"]),
    n_sig_padj_005 = sum(!is.na(res_df$padj) & res_df$padj < 0.05),
    stringsAsFactors = FALSE
  )
}

if (length(run_summary) > 0) {
  run_summary_df <- bind_rows(run_summary)
  write.csv(run_summary_df, file.path(outdir, "DESeq2_run_summary.csv"), row.names = FALSE)
  print(run_summary_df)
} else {
  message("No analyses were run.")
}

