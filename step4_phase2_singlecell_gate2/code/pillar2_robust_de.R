#!/usr/bin/env Rscript
# Pillar 2 Step 5 (R part) — robustness DE rankings for PT AKI-vs-normal.
# Produces logFC rankings (feature_id + symbol) for:
#   (1) per-dataset within-atlas AKI-vs-normal (a12ccb9b, dea717d4)  -> 2-dataset rule
#   (2) leave-one-dataset-out (drop each contributing dataset, pooled)
#   (3) down-sampled control: subsample normal donors to n=AKI, B=1000, edgeR collapsed
#       -> emits a MATRIX of EV-ranking-relevant logFC per iteration is too big;
#          instead we save per-iteration EV_TIER1 median logFC + full-gene logFC needed
#          for GSEA is heavy, so we save per-iteration logFC vectors ONLY for genes in
#          the union of all frozen gene sets (keeps files small, GSEA still valid on that
#          restricted axis? -> NO. GSEA needs full axis. So instead we compute NES in
#          Python by re-ranking; here we just save the per-iteration FULL logFC vector as
#          a compact feather-like TSV.) -> We save full logFC per iteration to parquet-free
#          TSV.GZ in a long format is huge. Compromise: save per-iteration logFC as columns.
# Engine note: primary contrast used limma-voom+dupCor; for the down-sample LOOP we use
# edgeR-QLF on donor-COLLAPSED counts (fast; voom-edgeR logFC r=0.95 for PT). Documented.
suppressMessages({library(edgeR); library(limma)})
set.seed(1)

DURABLE <- "/mnt/shared-workspace/shared/pillar2"
IN <- file.path(DURABLE, "de_inputs")
OUTL <- "/workspace/pillar2_robust_out"; dir.create(OUTL, showWarnings=FALSE, recursive=TRUE)
gene_map <- read.delim(file.path(DURABLE,"gene_map.tsv"), stringsAsFactors=FALSE)
rownames(gene_map) <- gene_map$feature_id

cnt <- as.matrix(read.delim(gzfile(file.path(IN,"counts_PT.tsv.gz")), row.names=1, check.names=FALSE))
cold <- read.delim(file.path(IN,"coldata_PT.tsv"), stringsAsFactors=FALSE); rownames(cold) <- cold$sample_id
stopifnot(all(colnames(cnt)==rownames(cold)))

# collapse to donor-level (sum suspensions) — used by per-dataset, LOO, downsample
collapse_donor <- function(y, meta){
  ud <- unique(meta$donor_id)
  yc <- sapply(ud, function(d) rowSums(y[, meta$donor_id==d, drop=FALSE]))
  yc <- as.matrix(yc); colnames(yc) <- ud
  md <- data.frame(donor_id=ud, stringsAsFactors=FALSE); rownames(md) <- ud
  md$disease_code <- meta$disease_code[match(ud, meta$donor_id)]
  md$dataset_id   <- meta$dataset_id[match(ud, meta$donor_id)]
  list(y=yc, md=md)
}

# edgeR QLF AKI vs normal on collapsed counts; returns data.frame(feature_id,symbol,logFC,FDR)
edger_aki_normal <- function(y, md, with_dataset=TRUE){
  md$disease_code <- factor(md$disease_code, levels=c("normal","AKI"))
  if (with_dataset && nlevels(factor(md$dataset_id))>=2) {
    md$dataset_id <- factor(md$dataset_id)
    mm <- model.matrix(~ dataset_id + disease_code, data=md)
    if (qr(mm)$rank < ncol(mm)) mm <- model.matrix(~ disease_code, data=md)
  } else mm <- model.matrix(~ disease_code, data=md)
  dge <- DGEList(counts=y[, rownames(md), drop=FALSE])
  keepg <- filterByExpr(dge, group=md$disease_code)
  dge <- dge[keepg,,keep.lib.sizes=FALSE]; dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, mm); fit <- glmQLFit(dge, mm)
  qlf <- glmQLFTest(fit, coef="disease_codeAKI")
  tt <- topTags(qlf, n=Inf, sort.by="none")$table
  data.frame(feature_id=rownames(tt), symbol=gene_map[rownames(tt),"feature_name"],
             logFC=tt$logFC, FDR=tt$FDR, stringsAsFactors=FALSE)
}

## ---------- (1) per-dataset within-atlas AKI-vs-normal ----------
per_ds_datasets <- c("a12ccb9b-4fbe-457d-8590-ac78053259ef",
                     "dea717d4-7bc0-4e46-950f-fd7e1cc8df7d")
for (ds in per_ds_datasets) {
  m <- cold[cold$dataset_id==ds & cold$disease_code %in% c("AKI","normal"), ]
  cp <- collapse_donor(cnt[, rownames(m), drop=FALSE], m)
  res <- edger_aki_normal(cp$y, cp$md, with_dataset=FALSE)  # single dataset
  write.table(res, file.path(OUTL, paste0("perdataset_PT_AKIvsnormal_", substr(ds,1,8),".tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("per-dataset", substr(ds,1,8), ": n_AKI=", sum(cp$md$disease_code=="AKI"),
      " n_normal=", sum(cp$md$disease_code=="normal"), " genes=", nrow(res), "\n")
}

## ---------- (2) leave-one-dataset-out (pooled, collapsed) ----------
contrib <- cold[cold$disease_code %in% c("AKI","normal"), ]
ds_all <- unique(contrib$dataset_id)
for (drop_ds in c("__none__", ds_all)) {
  m <- contrib[contrib$dataset_id != drop_ds, ]
  # need both AKI and normal present after drop
  if (length(unique(m$disease_code[m$disease_code=="AKI"]))==0 ||
      sum(m$disease_code=="AKI")<3 || sum(m$disease_code=="normal")<3) {
    cat("LOO drop", substr(drop_ds,1,8), ": insufficient after drop, skip\n"); next }
  cp <- collapse_donor(cnt[, rownames(m), drop=FALSE], m)
  res <- edger_aki_normal(cp$y, cp$md, with_dataset=TRUE)
  tag <- if(drop_ds=="__none__") "ALL" else substr(drop_ds,1,8)
  write.table(res, file.path(OUTL, paste0("loo_PT_AKIvsnormal_drop_", tag,".tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
  cat("LOO drop", tag, ": n_AKI=", sum(cp$md$disease_code=="AKI"),
      " n_normal=", sum(cp$md$disease_code=="normal"), " genes=", nrow(res), "\n")
}

## ---------- (3) down-sampled control: subsample normal to n_AKI, B=1000 ----------
# Save per-iteration FULL logFC vectors (genes x B) as compact TSV.GZ for Python GSEA.
contrib <- cold[cold$disease_code %in% c("AKI","normal"), ]
cp_all <- collapse_donor(cnt[, rownames(contrib), drop=FALSE], contrib)
aki_donors <- cp_all$md$donor_id[cp_all$md$disease_code=="AKI"]
norm_donors <- cp_all$md$donor_id[cp_all$md$disease_code=="normal"]
n_aki <- length(aki_donors)
B <- 1000
cat("\ndown-sample: n_AKI=", n_aki, " normal pool=", length(norm_donors), " B=", B, "\n")
# reference full gene axis = union of genes tested in the ALL-normal collapsed run
ref <- edger_aki_normal(cp_all$y, cp_all$md, with_dataset=TRUE)
gene_axis <- ref$feature_id
mat <- matrix(NA_real_, nrow=length(gene_axis), ncol=B,
              dimnames=list(gene_axis, paste0("iter",1:B)))
t0 <- Sys.time()
for (b in 1:B) {
  sub_norm <- sample(norm_donors, n_aki, replace=FALSE)
  keep_d <- c(aki_donors, sub_norm)
  md <- cp_all$md[cp_all$md$donor_id %in% keep_d, ]
  y <- cp_all$y[, md$donor_id, drop=FALSE]
  r <- tryCatch(edger_aki_normal(y, md, with_dataset=TRUE), error=function(e) NULL)
  if (!is.null(r)) mat[r$feature_id, b] <- r$logFC
  if (b %% 200 == 0) cat("  iter", b, " elapsed", round(difftime(Sys.time(),t0,units="mins"),1),"min\n")
}
# write with symbol column for Python GSEA
outm <- data.frame(feature_id=gene_axis, symbol=gene_map[gene_axis,"feature_name"],
                   mat, check.names=FALSE)
gz <- gzfile(file.path(OUTL,"downsample_PT_AKIvsnormal_logFC.tsv.gz"),"w")
write.table(outm, gz, sep="\t", quote=FALSE, row.names=FALSE); close(gz)
cat("wrote down-sample logFC matrix:", length(gene_axis),"genes x",B,"iters\n")

# copy durable via shell cp later (R file.copy 0-byte bug); just report
cat("\nDONE robust DE ->", OUTL, "\n")
