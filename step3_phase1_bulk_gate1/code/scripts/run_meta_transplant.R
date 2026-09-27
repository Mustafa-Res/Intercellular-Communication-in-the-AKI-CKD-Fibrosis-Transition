#!/usr/bin/env Rscript
# Analysis B: transplant-ATN random-effects meta-analysis (GSE30718 + GSE53769)
# Per-gene inverse-variance random-effects meta on log2FC effect sizes.
# Primary model = DerSimonian-Laird (DL); sensitivity = REML.
# Effect = logFC, variance = (moderated SE)^2 (from limma), per approved plan.
suppressMessages({library(metafor)})

args <- commandArgs(trailingOnly = TRUE)
# args: <de1.tsv> <de2.tsv> <out.tsv>
de1_path <- args[1]; de2_path <- args[2]; out_path <- args[3]

de1 <- read.delim(de1_path, stringsAsFactors = FALSE)
de2 <- read.delim(de2_path, stringsAsFactors = FALSE)

# Keep gene, logFC, variance; ensure finite
prep <- function(d) {
  d <- d[, c("gene", "logFC", "variance")]
  d <- d[is.finite(d$logFC) & is.finite(d$variance) & d$variance > 0, ]
  d[!duplicated(d$gene), ]
}
de1 <- prep(de1); de2 <- prep(de2)

shared <- intersect(de1$gene, de2$gene)
cat(sprintf("Shared genes for meta: %d\n", length(shared)))
rownames(de1) <- de1$gene; rownames(de2) <- de2$gene

res <- data.frame(gene = shared,
                  yi1 = de1[shared, "logFC"], vi1 = de1[shared, "variance"],
                  yi2 = de2[shared, "logFC"], vi2 = de2[shared, "variance"],
                  stringsAsFactors = FALSE)

# Per-gene random-effects meta
out <- data.frame(gene = shared, pooled_logFC = NA_real_, se = NA_real_,
                  ci_lb = NA_real_, ci_ub = NA_real_, zval = NA_real_,
                  pval = NA_real_, Q = NA_real_, Qp = NA_real_, I2 = NA_real_,
                  tau2 = NA_real_, pooled_logFC_reml = NA_real_, pval_reml = NA_real_,
                  stringsAsFactors = FALSE)

for (i in seq_along(shared)) {
  yi <- c(res$yi1[i], res$yi2[i]); vi <- c(res$vi1[i], res$vi2[i])
  fit <- tryCatch(rma(yi = yi, vi = vi, method = "DL"), error = function(e) NULL)
  if (!is.null(fit)) {
    out$pooled_logFC[i] <- as.numeric(fit$b)
    out$se[i] <- fit$se; out$ci_lb[i] <- fit$ci.lb; out$ci_ub[i] <- fit$ci.ub
    out$zval[i] <- fit$zval; out$pval[i] <- fit$pval
    out$Q[i] <- fit$QE; out$Qp[i] <- fit$QEp
    out$I2[i] <- fit$I2; out$tau2[i] <- fit$tau2
  }
  fit2 <- tryCatch(rma(yi = yi, vi = vi, method = "REML"), error = function(e) NULL)
  if (!is.null(fit2)) { out$pooled_logFC_reml[i] <- as.numeric(fit2$b); out$pval_reml[i] <- fit2$pval }
  if (i %% 3000 == 0) cat(sprintf("  ...%d/%d\n", i, length(shared)))
}

out$adj_pval <- p.adjust(out$pval, method = "BH")
out$adj_pval_reml <- p.adjust(out$pval_reml, method = "BH")
out <- out[order(out$adj_pval), ]
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote %d genes -> %s\n", nrow(out), out_path))
cat(sprintf("Sig (BH<0.05) DL: %d | REML: %d\n",
            sum(out$adj_pval < 0.05, na.rm = TRUE),
            sum(out$adj_pval_reml < 0.05, na.rm = TRUE)))
cat(sprintf("Median I2: %.1f%%\n", median(out$I2, na.rm = TRUE)))
