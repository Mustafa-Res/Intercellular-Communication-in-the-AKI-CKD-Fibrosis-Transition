#!/usr/bin/env Rscript
# Pillar 2 Step 5 fix — down-sampled control ONLY (per-dataset & LOO already done).
# Bug fix: assign logFC by intersecting each iteration's genes with the reference axis
# (subsample can pass filterByExpr on genes absent from the full-run reference, or vice versa).
suppressMessages({library(edgeR); library(limma)})
set.seed(1)
DURABLE <- "/mnt/shared-workspace/shared/pillar2"; IN <- file.path(DURABLE,"de_inputs")
OUTL <- "/workspace/pillar2_robust_out"; dir.create(OUTL, showWarnings=FALSE, recursive=TRUE)
gene_map <- read.delim(file.path(DURABLE,"gene_map.tsv"), stringsAsFactors=FALSE)
rownames(gene_map) <- gene_map$feature_id
cnt <- as.matrix(read.delim(gzfile(file.path(IN,"counts_PT.tsv.gz")), row.names=1, check.names=FALSE))
cold <- read.delim(file.path(IN,"coldata_PT.tsv"), stringsAsFactors=FALSE); rownames(cold) <- cold$sample_id

collapse_donor <- function(y, meta){
  ud <- unique(meta$donor_id)
  yc <- sapply(ud, function(d) rowSums(y[, meta$donor_id==d, drop=FALSE]))
  yc <- as.matrix(yc); colnames(yc) <- ud
  md <- data.frame(donor_id=ud, stringsAsFactors=FALSE); rownames(md) <- ud
  md$disease_code <- meta$disease_code[match(ud, meta$donor_id)]
  md$dataset_id   <- meta$dataset_id[match(ud, meta$donor_id)]
  list(y=yc, md=md)
}
edger_aki_normal <- function(y, md){
  md$disease_code <- factor(md$disease_code, levels=c("normal","AKI"))
  md$dataset_id <- factor(md$dataset_id)
  mm <- if (nlevels(md$dataset_id)>=2) model.matrix(~ dataset_id + disease_code, data=md) else model.matrix(~ disease_code, data=md)
  if (qr(mm)$rank < ncol(mm)) mm <- model.matrix(~ disease_code, data=md)
  dge <- DGEList(counts=y[, rownames(md), drop=FALSE])
  keepg <- filterByExpr(dge, group=md$disease_code)
  dge <- dge[keepg,,keep.lib.sizes=FALSE]; dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, mm); fit <- glmQLFit(dge, mm)
  qlf <- glmQLFTest(fit, coef="disease_codeAKI")
  tt <- topTags(qlf, n=Inf, sort.by="none")$table
  setNames(tt$logFC, rownames(tt))
}

contrib <- cold[cold$disease_code %in% c("AKI","normal"), ]
cp_all <- collapse_donor(cnt[, rownames(contrib), drop=FALSE], contrib)
aki_donors  <- cp_all$md$donor_id[cp_all$md$disease_code=="AKI"]
norm_donors <- cp_all$md$donor_id[cp_all$md$disease_code=="normal"]
n_aki <- length(aki_donors); B <- 1000
ref_lfc <- edger_aki_normal(cp_all$y, cp_all$md)
gene_axis <- names(ref_lfc)
mat <- matrix(NA_real_, nrow=length(gene_axis), ncol=B, dimnames=list(gene_axis, paste0("iter",1:B)))
cat("down-sample: n_AKI=",n_aki," normal pool=",length(norm_donors)," B=",B," ref genes=",length(gene_axis),"\n")
t0 <- Sys.time()
for (b in 1:B) {
  keep_d <- c(aki_donors, sample(norm_donors, n_aki, replace=FALSE))
  md <- cp_all$md[cp_all$md$donor_id %in% keep_d, ]
  y <- cp_all$y[, md$donor_id, drop=FALSE]
  lfc <- tryCatch(edger_aki_normal(y, md), error=function(e) NULL)
  if (!is.null(lfc)) {
    common <- intersect(names(lfc), gene_axis)     # <-- FIX: intersect before assign
    mat[common, b] <- lfc[common]
  }
  if (b %% 200 == 0) cat("  iter",b," elapsed",round(difftime(Sys.time(),t0,units="mins"),1),"min\n")
}
outm <- data.frame(feature_id=gene_axis, symbol=gene_map[gene_axis,"feature_name"], mat, check.names=FALSE)
gz <- gzfile(file.path(OUTL,"downsample_PT_AKIvsnormal_logFC.tsv.gz"),"w")
write.table(outm, gz, sep="\t", quote=FALSE, row.names=FALSE); close(gz)
cat("wrote down-sample logFC matrix:",length(gene_axis),"genes x",B,"iters\n")
