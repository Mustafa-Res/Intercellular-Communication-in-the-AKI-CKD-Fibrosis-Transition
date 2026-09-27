#!/usr/bin/env Rscript
# Map GPL16686 (Affymetrix Human Gene 2.0 ST) probe IDs -> gene symbols via
# RefSeq accessions (GB_ACC) using org.Hs.eg.db. Writes feature->gene_symbol map.
suppressMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})
args <- commandArgs(trailingOnly = TRUE)
acc_path <- args[1]
out_path <- args[2]

acc <- read.table(acc_path, header = TRUE, sep = "\t", colClasses = "character")
acc$acc_nov <- sub("\\..*$", "", acc$GB_ACC)  # strip version suffix
m <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(acc$acc_nov),
                           keytype = "REFSEQ", columns = "SYMBOL")
m <- m[!is.na(m$SYMBOL), ]
merged <- merge(acc, m, by.x = "acc_nov", by.y = "REFSEQ")
out <- unique(data.frame(feature = merged$ID, gene_symbol = merged$SYMBOL))
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("mapped features:", nrow(out), "| unique genes:",
    length(unique(out$gene_symbol)), "\n")
