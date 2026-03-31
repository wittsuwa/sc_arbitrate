library(Seurat)
library(Matrix)

amp_dir <- "/QRISdata/Q7458/Projects/ARBITRATE_singlecell/AMP Data"

export_seurat <- function(rds_path, out_dir, prefix) {
  cat("Loading:", rds_path, "\n")
  obj <- readRDS(rds_path)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Add prefix to cluster names
  obj$clusters_names <- paste0(prefix, obj$clusters_names)
  cat("Cell types:\n")
  print(table(obj$clusters_names))

  # Raw counts matrix (genes x cells sparse)
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  writeMM(counts, file = file.path(out_dir, "counts.mtx"))

  # Gene names and cell barcodes
  write.csv(data.frame(gene = rownames(counts)),
            file = file.path(out_dir, "genes.csv"), row.names = FALSE)
  write.csv(data.frame(barcode = colnames(counts)),
            file = file.path(out_dir, "barcodes.csv"), row.names = FALSE)

  # Cell metadata — at minimum clusters_names
  meta <- obj@meta.data
  write.csv(meta, file = file.path(out_dir, "metadata.csv"), row.names = TRUE)

  cat("Exported to:", out_dir, "\n\n")
}

export_seurat(
  rds_path = file.path(amp_dir, "dunlapetal_finalCD4.rds"),
  out_dir  = file.path(amp_dir, "dunlapetal_finalCD4_export"),
  prefix   = "CD4_"
)

export_seurat(
  rds_path = file.path(amp_dir, "dunlapetal_finalCD8.rds"),
  out_dir  = file.path(amp_dir, "dunlapetal_finalCD8_export"),
  prefix   = "CD8_"
)

cat("Done.\n")