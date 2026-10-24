#!/usr/bin/env Rscript
# 40e_revision_multivariable.R -- R3 major 2: gene sets defined for multivariable prediction, evaluated as multivariable
# classifiers rather than as single genes.
# For each set (LASSO 22, SVM-RFE 11, bootstrap-stable genes with selection frequency >= 0.5, core 4 for reference),
# a ridge logistic model is fitted on the pooled discovery cohorts (z-scored within cohort, the same input as the
# selection step) and applied unchanged to TCGA-UCEC and GSE146889 (z-scored within cohort).
# Null: the identical pipeline (refit with internal cross-validation) on random gene sets of equal size drawn from
# (a) the recurrent consensus (selection-matched) and (b) non-consensus genes of the discovery universe.
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(glmnet)})
set.seed(42); NNULL <- 500

fast_auc <- function(x, isT) { r <- rank(x); nT <- sum(isT); nN <- sum(!isT); (sum(r[isT]) - nT * (nT + 1) / 2) / (nT * nN) }

## ---- discovery matrix for every gene of the shared universe (same construction as R/09)
MAT <- c(GSE17025 = "data_processed/expression/GSE17025_gene_matrix_bestprobe.tsv",
         GSE63678 = "data_processed/expression/GSE63678_gene_matrix_bestprobe.tsv",
         GSE106191 = "data_processed/expression/GSE106191_all97_gene_matrix_bestprobe.tsv")
man <- read.delim("metadata/sample_manifest.tsv", stringsAsFactors = FALSE)
m6  <- read.delim("metadata/sample_manifest_GSE106191.tsv", stringsAsFactors = FALSE, quote = "")
V <- names(MAT)
L <- lapply(V, function(x) { t <- read.delim(sprintf("results/DEG/%s_all_genes.tsv", x), stringsAsFactors = FALSE); t <- t[!duplicated(t$gene), ]; t$gene })
uni <- Reduce(intersect, L)
blocks <- lapply(V, function(g) {
  x <- as.matrix(read.delim(MAT[g], row.names = 1, check.names = FALSE))
  mm <- if (g == "GSE106191") m6 else subset(man, dataset == g & include_primary == "yes")
  x <- x[intersect(uni, rownames(x)), mm$sample_id]
  list(z = t(scale(t(x))), y = as.integer(mm$normal_or_tumor == "tumor"))
})
gD <- Reduce(intersect, lapply(blocks, function(b) rownames(b$z)))
XD <- t(do.call(cbind, lapply(blocks, function(b) b$z[gD, ]))); yD <- unlist(lapply(blocks, `[[`, "y"))
cat(sprintf("discovery: %d samples x %d genes | tumour %d / control %d\n", nrow(XD), ncol(XD), sum(yD), sum(!yD)))
ml <- read.delim("results/ML/ML_input_matrix.tsv", check.names = FALSE)
chk <- max(abs(XD[ml$sample, "NUMA1"] - ml$NUMA1)); cat(sprintf("matches R/09 input (max |diff| NUMA1) = %.2e\n", chk)); stopifnot(chk < 1e-8)

## ---- validation cohorts, z-scored within cohort
tc <- function() {
  se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
  cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop = FALSE]); code <- substr(cdl$barcode, 14, 15)
  EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
  SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
  keep <- (code == "01" & cdl$primary_diagnosis %in% c(EN, SR)) | code == "11"
  cnt <- assay(se, "unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
  o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]; cnt <- cnt[i, ]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt) != "", ]
  grp <- factor(ifelse(code[keep] == "11", "N", "T"))
  d <- DGEList(cnt); d <- d[filterByExpr(d, group = grp), , keep.lib.sizes = FALSE]; d <- calcNormFactors(d)
  list(z = t(scale(t(cpm(d, log = TRUE)))), isT = grp == "T")
}
g146 <- function() {
  man <- subset(read.delim("metadata/sample_manifest_GSE146889.tsv", stringsAsFactors = FALSE), tissue == "Endometrial")
  raw <- read.delim("data_raw/GEO/GSE146889/GSE146889_GeneCount.tsv.gz", check.names = FALSE)
  gc <- raw[, paste0(man$column, "_count")]; colnames(gc) <- man$column
  sym <- raw$GeneName; ord <- order(sym, -rowMeans(gc)); idx <- ord[!duplicated(sym[ord])]
  m <- as.matrix(gc[idx, ]); rownames(m) <- sym[idx]; m <- m[!is.na(rownames(m)) & rownames(m) != "", ]
  gv <- factor(man$group); d <- DGEList(m); d <- d[filterByExpr(d, group = gv), , keep.lib.sizes = FALSE]; d <- calcNormFactors(d)
  list(z = t(scale(t(cpm(d, log = TRUE)))), isT = gv == "tumor")
}
VAL <- list(`TCGA-UCEC` = tc(), GSE146889 = g146())

## ---- the multivariable pipeline: ridge logistic, lambda by 5-fold CV on discovery only
measurable <- Reduce(intersect, c(list(colnames(XD)), lapply(VAL, function(v) rownames(v$z))))
fit_eval <- function(genes) {                     # one model per gene set, fitted on genes measurable in every cohort
  gg <- intersect(genes, measurable)
  cv <- cv.glmnet(XD[, gg, drop = FALSE], yD, family = "binomial", alpha = 0, nfolds = 5, type.measure = "deviance")
  b  <- as.numeric(coef(cv, s = "lambda.min"))[-1]
  vapply(names(VAL), function(nm) fast_auc(as.numeric(t(VAL[[nm]]$z[gg, , drop = FALSE]) %*% b), VAL[[nm]]$isT), 0)
}

rec  <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors = FALSE)$gene
las  <- read.delim("results/ML/lasso_selected_genes.tsv")$gene
svm  <- subset(read.delim("results/ML/svm_rfe_selected_genes.tsv"), selected %in% c("True", TRUE))$gene
stab <- subset(read.delim("results/ML/bootstrap_stability_recurrent.tsv"), freq_intersection >= 0.5)$gene
SETS <- list(`LASSO-selected` = las, `SVM-RFE-selected` = svm, `bootstrap-stable (frequency >= 0.5)` = stab,
             `core four genes (reference)` = c("HAND2-AS1","OSR2","PEG3","RUNX1T1"))
pool_sel <- intersect(rec, measurable); pool_non <- intersect(setdiff(uni, rec), measurable)

res <- do.call(rbind, lapply(names(SETS), function(s) {
  gs <- SETS[[s]]; k <- length(intersect(gs, measurable)); obs <- fit_eval(gs)
  ns <- replicate(NNULL, fit_eval(sample(pool_sel, k))); nn <- replicate(NNULL, fit_eval(sample(pool_non, k)))
  do.call(rbind, lapply(names(VAL), function(v) data.frame(
    gene_set = s, n_genes = k, cohort = v, AUC = round(obs[v], 4),
    selmatched_null_median = round(median(ns[v, ], na.rm = TRUE), 4),
    selmatched_P = signif((sum(ns[v, ] >= obs[v], na.rm = TRUE) + 1) / (sum(!is.na(ns[v, ])) + 1), 3),
    noncons_null_median = round(median(nn[v, ], na.rm = TRUE), 4),
    noncons_P = signif((sum(nn[v, ] >= obs[v], na.rm = TRUE) + 1) / (sum(!is.na(nn[v, ])) + 1), 3))))
}))
print(res, row.names = FALSE)
write.table(res, "results/revision/multivariable_gene_sets.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("===== 40e DONE =====\n")
