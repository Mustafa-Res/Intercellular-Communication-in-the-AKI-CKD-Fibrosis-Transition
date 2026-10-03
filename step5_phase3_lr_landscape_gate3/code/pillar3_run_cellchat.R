#!/usr/bin/env Rscript
# Pillar 3 CellChat engine (one shard = one dataset x condition).
# Ingests the forR_<tag> bundle written by pillar3_prep_and_liana.py (IDENTICAL cells
# to the LIANA run), builds a Seurat object with normalized data, then FAITHFULLY calls
# the skill's run_cellchat_analysis() (sourced, not reimplemented) with db_category="all".
# Extracts PT->receiver LR tables + saves the CellChat object for the differential step.
#
# Usage: Rscript pillar3_run_cellchat.R <forR_dir> <out_dir> <tag>

.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(CellChat)
})
args <- commandArgs(trailingOnly = TRUE)
forR <- args[1]; out_dir <- args[2]; tag <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
SKILL <- "/mnt/skills/system/cell-cell-communication/01_bulk_signature/code/scripts/run_cellchat.R"
source(SKILL)  # provides run_cellchat_analysis()

RECEIVERS <- c("receiver_myeloid", "receiver_stromal", "receiver_endothelial")

# ---- load bundle ----
cat("Loading bundle", forR, "\n")
data_gc  <- as(Matrix::readMM(gzfile(file.path(forR, "matrix.mtx.gz"))), "CsparseMatrix")   # genes x cells (normalized)
counts_gc<- as(Matrix::readMM(gzfile(file.path(forR, "counts.mtx.gz"))), "CsparseMatrix")   # genes x cells (raw)
genes <- readLines(file.path(forR, "genes.tsv"))
barcodes <- readLines(file.path(forR, "barcodes.tsv"))
meta <- read.table(file.path(forR, "meta.tsv"), header = TRUE, sep = "\t",
                   row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
rownames(data_gc) <- genes; colnames(data_gc) <- barcodes
rownames(counts_gc) <- genes; colnames(counts_gc) <- barcodes
meta <- meta[barcodes, , drop = FALSE]
cat("  ", length(barcodes), "cells x", length(genes), "genes\n")
cat("  compartments:", paste(names(table(meta$compartment)), table(meta$compartment), sep="=", collapse=" "), "\n")

# ---- build Seurat object with normalized data in the 'data' slot ----
so <- CreateSeuratObject(counts = counts_gc, meta.data = meta)
so <- SetAssayData(so, layer = "data", new.data = data_gc)   # our precomputed log-norm
so$celltype <- factor(as.character(so$compartment),
                      levels = c("PT", RECEIVERS))
Idents(so) <- "celltype"

# ---- run CellChat via the skill function (all signaling categories) ----
cc <- run_cellchat_analysis(so, species = "human", group.by = "celltype",
                            db_category = "all", min.cells = 10)

# ---- save the full object (needed for mergeCellChat differential) ----
saveRDS(cc, file.path(out_dir, paste0("cellchat_", tag, ".rds")))

# ---- extract PT->receiver LR tables (source=PT) ----
present <- levels(cc@idents)
rec_present <- intersect(RECEIVERS, present)
if (!("PT" %in% present) || length(rec_present) == 0) {
  cat("WARNING: PT or receivers missing in", tag, "\n")
}
df <- tryCatch(
  subsetCommunication(cc, sources.use = "PT", targets.use = rec_present),
  error = function(e) { cat("subsetCommunication error:", conditionMessage(e), "\n"); NULL })
if (!is.null(df) && nrow(df) > 0) {
  write.table(df, file.path(out_dir, paste0("cellchat_LR_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  wrote", nrow(df), "PT->receiver LR rows\n")
} else {
  # still write an empty stub with header so downstream is robust
  write.table(data.frame(source=character(), target=character(), ligand=character(),
                         receptor=character(), prob=numeric(), pval=numeric(),
                         interaction_name=character(), pathway_name=character()),
              file.path(out_dir, paste0("cellchat_LR_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  wrote EMPTY PT->receiver table for", tag, "\n")
}

# ---- also dump pathway-level PT->receiver (netP) ----
dfp <- tryCatch(
  subsetCommunication(cc, sources.use = "PT", targets.use = rec_present, slot.name = "netP"),
  error = function(e) NULL)
if (!is.null(dfp) && nrow(dfp) > 0) {
  write.table(dfp, file.path(out_dir, paste0("cellchat_pathway_", tag, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  wrote", nrow(dfp), "PT->receiver pathway rows\n")
}
cat("CellChat DONE for", tag, "\n")
