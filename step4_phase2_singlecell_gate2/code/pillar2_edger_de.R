#!/usr/bin/env Rscript
# Pillar 2 Step 3 — edgeR QLF pseudobulk DE per compartment, 3 contrasts.
# Model: ~ dataset_id + suspension_type + disease_code  (drop collinear cols w/ logged reason)
# Primary: AKI vs normal | Secondary: CKD vs normal | Exploratory: AKI vs CKD
# edgeR QLF primary + limma-voom sensitivity. filterByExpr low-count filter.
suppressMessages({library(edgeR); library(limma)})

DURABLE <- "/mnt/shared-workspace/shared/pillar2"
IN  <- file.path(DURABLE, "de_inputs")
OUTL <- "/workspace/pillar2_de_out"          # local write
dir.create(OUTL, showWarnings = FALSE, recursive = TRUE)
gene_map <- read.delim(file.path(DURABLE, "gene_map.tsv"), stringsAsFactors = FALSE)
rownames(gene_map) <- gene_map$feature_id

COMPARTMENTS <- c("PT","tubular_other","receiver_myeloid","receiver_stromal",
                  "receiver_endothelial","injury_context")
# contrast defs: name -> c(numerator, denominator) in disease_code
CONTRASTS <- list(
  AKI_vs_normal = c("AKI","normal"),   # PRIMARY
  CKD_vs_normal = c("CKD","normal"),   # SECONDARY
  AKI_vs_CKD    = c("AKI","CKD")       # EXPLORATORY
)
MIN_DONORS_PER_GROUP <- 3   # need >=3 per side to attempt a contrast (design safety)

log_lines <- c()
lg <- function(...) { s <- paste0(...); log_lines[[length(log_lines)+1]] <<- s; cat(s, "\n") }

# helper: drop covariate columns that are collinear/aliased/zero-variance in a design
# returns list(formula_terms=..., dropped=...) operating on a data.frame `md`
build_design <- function(md) {
  # candidate covariates in priority order (disease_code always last/kept)
  covars <- c("dataset_id","suspension_type")
  keep <- c()
  dropped <- c()
  for (cv in covars) {
    x <- factor(as.character(md[[cv]]))
    if (nlevels(x) < 2) { dropped <- c(dropped, paste0(cv,"(1 level)")); next }
    # confounded with disease? cross-tab: if each disease group maps to disjoint covar levels
    tab <- table(md$disease_code, x)
    # aliasing check done jointly below via model matrix rank; here just keep multi-level
    md[[cv]] <- x
    keep <- c(keep, cv)
  }
  md$disease_code <- factor(md$disease_code)  # 2 levels after subsetting
  terms <- c(keep, "disease_code")
  form <- as.formula(paste("~", paste(terms, collapse = " + ")))
  mm <- model.matrix(form, data = md)
  # detect rank deficiency -> drop covariates until full rank (protect disease_code)
  qrr <- qr(mm)
  if (qrr$rank < ncol(mm)) {
    # iteratively drop the lowest-priority covariate (suspension_type first, then dataset_id)
    for (cv in rev(keep)) {
      terms_try <- setdiff(terms, cv)
      form_try <- as.formula(paste("~", paste(terms_try, collapse = " + ")))
      mm_try <- model.matrix(form_try, data = md)
      if (qr(mm_try)$rank == ncol(mm_try)) {
        dropped <- c(dropped, paste0(cv,"(aliased/rank-deficient)"))
        terms <- terms_try; form <- form_try; mm <- mm_try
        break
      }
    }
    # if still deficient, drop both covariates
    if (qr(mm)$rank < ncol(mm)) {
      form <- ~ disease_code
      mm <- model.matrix(form, data = md)
      dropped <- unique(c(dropped, "all_covars(rank-deficient)"))
    }
  }
  list(form = form, mm = mm, dropped = dropped, md = md)
}

all_de <- list()
run_meta <- list()

for (comp in COMPARTMENTS) {
  cfile <- file.path(IN, paste0("counts_", comp, ".tsv.gz"))
  mfile <- file.path(IN, paste0("coldata_", comp, ".tsv"))
  if (!file.exists(cfile)) { lg("SKIP ", comp, " (no counts)"); next }
  cnt <- read.delim(gzfile(cfile), row.names = 1, check.names = FALSE)
  cnt <- as.matrix(cnt)
  cold <- read.delim(mfile, stringsAsFactors = FALSE)
  rownames(cold) <- cold$sample_id
  stopifnot(all(colnames(cnt) == rownames(cold)))
  lg("\n==== ", comp, " : ", ncol(cnt), " samples, ", nrow(cnt), " genes ====")

  for (cn in names(CONTRASTS)) {
    grp <- CONTRASTS[[cn]]; num <- grp[1]; den <- grp[2]
    keep_s <- cold$disease_code %in% grp
    md <- cold[keep_s, , drop = FALSE]
    md$disease_code <- factor(md$disease_code, levels = c(den, num))  # den=ref
    n_num <- sum(md$disease_code == num); n_den <- sum(md$disease_code == den)
    nd_num <- length(unique(md$donor_id[md$disease_code == num]))
    nd_den <- length(unique(md$donor_id[md$disease_code == den]))
    if (nd_num < MIN_DONORS_PER_GROUP || nd_den < MIN_DONORS_PER_GROUP) {
      lg("  [", cn, "] SKIP: donors ", num, "=", nd_num, " ", den, "=", nd_den); next
    }
    y <- cnt[, rownames(md), drop = FALSE]
    dge <- DGEList(counts = y)
    des <- build_design(md)
    mm <- des$mm; md2 <- des$md
    if (length(des$dropped) > 0)
      lg("  [", cn, "] dropped covars: ", paste(des$dropped, collapse=", "))
    # filterByExpr using disease group
    keepg <- filterByExpr(dge, group = md2$disease_code)
    dge <- dge[keepg, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, mm)
    # ---- edgeR QLF (primary) ----
    fit <- glmQLFit(dge, mm)
    coef_name <- paste0("disease_code", num)  # since ref=den
    qlf <- glmQLFTest(fit, coef = coef_name)
    tt <- topTags(qlf, n = Inf, sort.by = "none")$table
    res <- data.frame(feature_id = rownames(tt),
                      symbol = gene_map[rownames(tt), "feature_name"],
                      logFC = tt$logFC, logCPM = tt$logCPM,
                      F = tt$F, PValue = tt$PValue, FDR = tt$FDR,
                      stringsAsFactors = FALSE)
    # ---- limma-voom (sensitivity) ----
    v <- voom(dge, mm)
    vfit <- lmFit(v, mm)
    vfit <- eBayes(vfit)
    vtt <- topTable(vfit, coef = coef_name, n = Inf, sort.by = "none")
    res$voom_logFC <- vtt[res$feature_id, "logFC"]
    res$voom_P     <- vtt[res$feature_id, "P.Value"]
    res$voom_FDR   <- vtt[res$feature_id, "adj.P.Val"]

    fout <- file.path(OUTL, paste0("de_", comp, "__", cn, ".tsv"))
    write.table(res, fout, sep = "\t", quote = FALSE, row.names = FALSE)
    nsig <- sum(res$FDR < 0.05, na.rm = TRUE)
    lg("  [", cn, "] genes_tested=", nrow(res), " n_", num, "=", n_num,
       " n_", den, "=", n_den, " donors(", num, "/", den, ")=", nd_num, "/", nd_den,
       " sig(FDR<.05)=", nsig,
       " | edgeR-voom logFC cor=", round(cor(res$logFC, res$voom_logFC, use="complete.obs"),3))
    run_meta[[length(run_meta)+1]] <- data.frame(
      compartment=comp, contrast=cn, genes_tested=nrow(res),
      n_num=n_num, n_den=n_den, donors_num=nd_num, donors_den=nd_den,
      n_sig_fdr05=nsig, model=paste(as.character(des$form), collapse=" "),
      dropped_covars=paste(des$dropped, collapse=";"), stringsAsFactors=FALSE)
  }
}

rm_df <- do.call(rbind, run_meta)
write.table(rm_df, file.path(OUTL, "de_run_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
writeLines(unlist(log_lines), file.path(OUTL, "de_run_log.txt"))
cat("\n==== DE run summary ====\n")
print(rm_df)
cat("\nDONE. Outputs in", OUTL, "\n")
