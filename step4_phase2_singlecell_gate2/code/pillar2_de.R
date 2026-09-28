#!/usr/bin/env Rscript
# Pillar 2 Step 3 — pseudobulk DE per compartment, 3 contrasts.
# PRIMARY  : limma-voom + duplicateCorrelation (block=donor_id), full donor x suspension set.
# SENSITIVITY: edgeR-QLF on donor-collapsed counts (cell+nucleus summed -> 1 sample/donor).
# Model: ~ dataset_id + suspension_type + disease_code (drop collinear w/ logged reason).
# Primary contrast: AKI vs normal | Secondary: CKD vs normal | Exploratory: AKI vs CKD.
suppressMessages({library(edgeR); library(limma)})

DURABLE <- "/mnt/shared-workspace/shared/pillar2"
IN  <- file.path(DURABLE, "de_inputs")
OUTL <- "/workspace/pillar2_de_out"
dir.create(OUTL, showWarnings = FALSE, recursive = TRUE)
gene_map <- read.delim(file.path(DURABLE, "gene_map.tsv"), stringsAsFactors = FALSE)
rownames(gene_map) <- gene_map$feature_id

COMPARTMENTS <- c("PT","tubular_other","receiver_myeloid","receiver_stromal",
                  "receiver_endothelial","injury_context")
CONTRASTS <- list(
  AKI_vs_normal = c("AKI","normal"),   # PRIMARY
  CKD_vs_normal = c("CKD","normal"),   # SECONDARY
  AKI_vs_CKD    = c("AKI","CKD")       # EXPLORATORY
)
MIN_DONORS_PER_GROUP <- 3

log_lines <- c()
lg <- function(...) { s <- paste0(...); log_lines[[length(log_lines)+1]] <<- s; cat(s, "\n") }

# build a full-rank design from covariates + disease_code; drop aliased covars (log reason)
build_design <- function(md, covars = c("dataset_id","suspension_type")) {
  keep <- c(); dropped <- c()
  for (cv in covars) {
    x <- factor(as.character(md[[cv]]))
    if (nlevels(x) < 2) { dropped <- c(dropped, paste0(cv,"(1 level)")); next }
    md[[cv]] <- x; keep <- c(keep, cv)
  }
  md$disease_code <- factor(md$disease_code)
  terms <- c(keep, "disease_code")
  form <- as.formula(paste("~", paste(terms, collapse = " + ")))
  mm <- model.matrix(form, data = md)
  if (qr(mm)$rank < ncol(mm)) {
    for (cv in rev(keep)) {
      tt <- setdiff(terms, cv); ff <- as.formula(paste("~", paste(tt, collapse=" + ")))
      mmt <- model.matrix(ff, data = md)
      if (qr(mmt)$rank == ncol(mmt)) { dropped <- c(dropped, paste0(cv,"(aliased)"))
        terms <- tt; form <- ff; mm <- mmt; break }
    }
    if (qr(mm)$rank < ncol(mm)) { form <- ~ disease_code; mm <- model.matrix(form, data=md)
      dropped <- unique(c(dropped, "all_covars(rank-deficient)")) }
  }
  list(form = form, mm = mm, dropped = dropped, md = md)
}

run_meta <- list()

for (comp in COMPARTMENTS) {
  cfile <- file.path(IN, paste0("counts_", comp, ".tsv.gz"))
  mfile <- file.path(IN, paste0("coldata_", comp, ".tsv"))
  if (!file.exists(cfile)) { lg("SKIP ", comp); next }
  cnt <- as.matrix(read.delim(gzfile(cfile), row.names = 1, check.names = FALSE))
  cold <- read.delim(mfile, stringsAsFactors = FALSE); rownames(cold) <- cold$sample_id
  stopifnot(all(colnames(cnt) == rownames(cold)))
  lg("\n==== ", comp, " : ", ncol(cnt), " samples, ", nrow(cnt), " genes ====")

  for (cn in names(CONTRASTS)) {
    grp <- CONTRASTS[[cn]]; num <- grp[1]; den <- grp[2]
    md <- cold[cold$disease_code %in% grp, , drop = FALSE]
    md$disease_code <- factor(md$disease_code, levels = c(den, num))  # den = reference
    nd_num <- length(unique(md$donor_id[md$disease_code==num]))
    nd_den <- length(unique(md$donor_id[md$disease_code==den]))
    if (nd_num < MIN_DONORS_PER_GROUP || nd_den < MIN_DONORS_PER_GROUP) {
      lg("  [", cn, "] SKIP donors ", num,"=",nd_num," ",den,"=",nd_den); next }
    coef_name <- paste0("disease_code", num)

    ## ---------- PRIMARY: limma-voom + duplicateCorrelation ----------
    y <- cnt[, rownames(md), drop = FALSE]
    dge <- DGEList(counts = y)
    des <- build_design(md); mm <- des$mm; md2 <- des$md
    keepg <- filterByExpr(dge, group = md2$disease_code)
    dge <- dge[keepg, , keep.lib.sizes = FALSE]; dge <- calcNormFactors(dge)
    v <- voom(dge, mm)
    # within-donor correlation (block = donor)
    corfit <- duplicateCorrelation(v, mm, block = md2$donor_id)
    dcor <- corfit$consensus.correlation
    v <- voom(dge, mm, block = md2$donor_id, correlation = dcor)   # re-voom w/ block
    corfit <- duplicateCorrelation(v, mm, block = md2$donor_id)
    dcor <- corfit$consensus.correlation
    vfit <- lmFit(v, mm, block = md2$donor_id, correlation = dcor)
    vfit <- eBayes(vfit)
    vtt <- topTable(vfit, coef = coef_name, n = Inf, sort.by = "none")
    res <- data.frame(feature_id = rownames(vtt),
                      symbol = gene_map[rownames(vtt), "feature_name"],
                      logFC = vtt$logFC, AveExpr = vtt$AveExpr, t = vtt$t,
                      PValue = vtt$P.Value, FDR = vtt$adj.P.Val,
                      stringsAsFactors = FALSE)

    ## ---------- SENSITIVITY: edgeR-QLF on donor-collapsed counts ----------
    # collapse: sum counts across suspensions within donor
    donor_of <- md$donor_id
    ud <- unique(donor_of)
    yc <- sapply(ud, function(d) rowSums(y[, donor_of==d, drop=FALSE]))
    yc <- as.matrix(yc); colnames(yc) <- ud
    # collapsed coldata: disease per donor (unique), dataset per donor (mode), suspension = dominant
    cmd <- data.frame(donor_id = ud, stringsAsFactors = FALSE)
    rownames(cmd) <- ud
    cmd$disease_code <- md$disease_code[match(ud, donor_of)]
    cmd$disease_code <- factor(as.character(cmd$disease_code), levels=c(den,num))
    # dataset_id: donors are within one dataset (verified) -> take first
    cmd$dataset_id <- md$dataset_id[match(ud, donor_of)]
    # suspension: dominant modality by n_cells
    dom_susp <- sapply(ud, function(d){ sub <- md[donor_of==d,,drop=FALSE]
      sub$suspension_type[which.max(sub$n_cells)] })
    cmd$suspension_type <- dom_susp
    dge2 <- DGEList(counts = yc[, ud, drop=FALSE])
    des2 <- build_design(cmd); mm2 <- des2$mm; cmd2 <- des2$md
    keepg2 <- filterByExpr(dge2, group = cmd2$disease_code)
    dge2 <- dge2[keepg2, , keep.lib.sizes=FALSE]; dge2 <- calcNormFactors(dge2)
    dge2 <- estimateDisp(dge2, mm2)
    fit2 <- glmQLFit(dge2, mm2)
    qlf <- glmQLFTest(fit2, coef = coef_name)
    qtt <- topTags(qlf, n = Inf, sort.by = "none")$table
    res$edgeR_logFC <- qtt[res$feature_id, "logFC"]
    res$edgeR_FDR   <- qtt[res$feature_id, "FDR"]

    fout <- file.path(OUTL, paste0("de_", comp, "__", cn, ".tsv"))
    write.table(res, fout, sep = "\t", quote = FALSE, row.names = FALSE)
    nsig <- sum(res$FDR < 0.05, na.rm = TRUE)
    lfc_cor <- suppressWarnings(cor(res$logFC, res$edgeR_logFC, use="complete.obs"))
    lg("  [", cn, "] genes=", nrow(res), " voomN(", num,"/",den,")=",
       sum(md2$disease_code==num),"/",sum(md2$disease_code==den),
       " donors=",nd_num,"/",nd_den, " dupCor=", round(dcor,3),
       " sig(FDR<.05)=", nsig, " voom-edgeR logFCcor=", round(lfc_cor,3),
       if(length(des$dropped)) paste0(" | dropped:",paste(des$dropped,collapse=",")) else "")
    run_meta[[length(run_meta)+1]] <- data.frame(
      compartment=comp, contrast=cn, engine_primary="limma-voom+dupCor",
      genes_tested=nrow(res), n_samp_num=sum(md2$disease_code==num),
      n_samp_den=sum(md2$disease_code==den), donors_num=nd_num, donors_den=nd_den,
      dup_correlation=round(dcor,4), n_sig_fdr05=nsig,
      voom_edgeR_logFC_cor=round(lfc_cor,4),
      model=paste(as.character(des$form),collapse=" "),
      dropped_covars=paste(des$dropped,collapse=";"), stringsAsFactors=FALSE)
  }
}

rm_df <- do.call(rbind, run_meta)
write.table(rm_df, file.path(OUTL, "de_run_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
writeLines(unlist(log_lines), file.path(OUTL, "de_run_log.txt"))
cat("\n==== DE RUN SUMMARY ====\n"); print(rm_df)

# copy durable
dir.create(file.path(DURABLE,"de_out"), showWarnings=FALSE)
for (f in list.files(OUTL)) file.copy(file.path(OUTL,f), file.path(DURABLE,"de_out",f), overwrite=TRUE)
cat("\nDONE ->", file.path(DURABLE,"de_out"), "\n")
