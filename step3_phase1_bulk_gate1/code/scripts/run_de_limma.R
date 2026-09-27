#!/usr/bin/env Rscript
# Per-cohort differential expression with limma (PLAN v1.3/v1.6).
#
# Input : an expression matrix TSV (features x samples) + a sample-meta TSV
#         (sample, group in {case,control}[, covariates]) + a feature->gene map TSV
#         (feature, gene_symbol).  All produced by the Python loader.
# Output: per-gene DE table with log2FC, SE (for meta-analysis variance), t, P,
#         adj.P (BH), plus AveExpr and n_case/n_control.
#
# Design: ~ group (+ optional covariate columns passed via --covars).
# Effect-size convention for downstream random-effects meta-analysis:
#   effect   = logFC (case vs control)
#   variance = SE^2   (SE = logFC / t, recovered from limma moderated t)
#
# Data are assumed log2-scale on input (Python loader logs linear data first).

suppressMessages({
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NA) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) args[i + 1] else default
}

expr_path  <- getarg("--expr")
meta_path  <- getarg("--meta")
map_path   <- getarg("--map")     # optional; if NA, index already gene symbols
out_path   <- getarg("--out")
covars     <- getarg("--covars", "")  # comma-separated meta columns
cohort_id  <- getarg("--id", "cohort")

stopifnot(!is.na(expr_path), !is.na(meta_path), !is.na(out_path))

expr <- as.matrix(read.table(expr_path, header = TRUE, sep = "\t",
                             row.names = 1, check.names = FALSE))
meta <- read.table(meta_path, header = TRUE, sep = "\t",
                   row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# align
common <- intersect(colnames(expr), rownames(meta))
expr <- expr[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]
cat(sprintf("[%s] %d features x %d samples (%s)\n", cohort_id,
            nrow(expr), ncol(expr), paste(names(table(meta$group)),
            table(meta$group), sep = "=", collapse = ", ")))

# ---- collapse features -> gene symbol (max mean-expression probe per gene) ----
if (!is.na(map_path)) {
  map <- read.table(map_path, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "", comment.char = "",
                    colClasses = "character")
  colnames(map)[1:2] <- c("feature", "gene_symbol")
  map$feature <- as.character(map$feature)
  rownames(expr) <- as.character(rownames(expr))
  map <- map[map$gene_symbol != "" & !is.na(map$gene_symbol), ]
  map <- map[map$feature %in% rownames(expr), ]
  expr <- expr[map$feature, , drop = FALSE]
  gene <- map$gene_symbol
  # collapse: for each gene keep the probe with highest average expression
  ave <- rowMeans(expr, na.rm = TRUE)
  ord <- order(ave, decreasing = TRUE)
  expr <- expr[ord, , drop = FALSE]
  gene <- gene[ord]
  keep <- !duplicated(gene)
  expr <- expr[keep, , drop = FALSE]
  rownames(expr) <- gene[keep]
  cat(sprintf("[%s] collapsed to %d genes\n", cohort_id, nrow(expr)))
}

# ---- design ----
grp <- factor(meta$group, levels = c("control", "case"))
form_terms <- "grp"
if (nzchar(covars)) {
  cv <- strsplit(covars, ",")[[1]]
  for (c in cv) meta[[c]] <- factor(meta[[c]])
  form_terms <- c("grp", cv)
}
design <- model.matrix(as.formula(paste("~", paste(form_terms, collapse = " + "))),
                       data = cbind(grp = grp, meta))

# ---- limma ----
fit <- lmFit(expr, design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
coef_name <- "grpcase"
tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")

# recover SE = logFC / t  (moderated); variance = SE^2 for meta-analysis
se <- tt$logFC / tt$t
res <- data.frame(
  gene      = rownames(tt),
  logFC     = tt$logFC,
  SE        = se,
  variance  = se^2,
  t         = tt$t,
  P.Value   = tt$P.Value,
  adj.P.Val = tt$adj.P.Val,
  AveExpr   = tt$AveExpr,
  n_case    = sum(grp == "case"),
  n_control = sum(grp == "control"),
  cohort    = cohort_id,
  stringsAsFactors = FALSE
)
res <- res[order(res$P.Value), ]
write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[%s] wrote %d genes -> %s\n", cohort_id, nrow(res), out_path))
cat(sprintf("[%s] n sig (adj.P<0.05): %d\n", cohort_id, sum(res$adj.P.Val < 0.05, na.rm=TRUE)))
