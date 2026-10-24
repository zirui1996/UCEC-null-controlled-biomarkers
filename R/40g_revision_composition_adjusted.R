#!/usr/bin/env Rscript
# 40g_revision_composition_adjusted.R (v2) -- R3 major 5: tumour-versus-normal effects of the four core genes with and without a
# fibroblast marker score as a covariate, fitted in the same genome-wide voom-limma framework as the primary analyses:
#   TCGA-UCEC  voom on all filtered genes, design ~ group (as R/28b), +/- fibroblast score
#   GSE146889  voom on the 35 complete pairs, design ~ 0 + group + patient (as R/02d), +/- fibroblast score
# The fibroblast score is the mean z-score (log-CPM) of COL1A1, COL1A2 and DCN, the canonical fibroblast markers used for lineage
# assignment in the single-nucleus analysis. The output is the attenuation of the tumour coefficient after adjustment for this proxy;
# it does not estimate the share of the difference caused by composition.
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(limma)})
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1"); FIB <- c("COL1A1","COL1A2","DCN")
fast_auc <- function(x, isT) { r <- rank(x); nT <- sum(isT); nN <- sum(!isT); (sum(r[isT]) - nT * (nT + 1) / 2) / (nT * nN) }

tumour_coef <- function(v_fun, design, coef_name) {
  v <- v_fun(design); fit <- eBayes(lmFit(v, design)); tt <- topTable(fit, coef = coef_name, number = Inf)
  tt[CORE, c("logFC", "P.Value", "adj.P.Val")]
}

## ---------------------------------------------------------------- TCGA-UCEC
se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop = FALSE]); code <- substr(cdl$barcode, 14, 15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code == "01" & cdl$primary_diagnosis %in% c(EN, SR)) | code == "11"
cnt <- assay(se, "unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]; cnt <- cnt[i, ]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt) != "", ]
grp <- factor(ifelse(code[keep] == "11", "Normal", "Tumor"), levels = c("Normal","Tumor"))
dT <- DGEList(cnt); dT <- dT[filterByExpr(dT, group = grp), , keep.lib.sizes = FALSE]; dT <- calcNormFactors(dT)
lcT <- cpm(dT, log = TRUE); fibT <- colMeans(t(scale(t(lcT[FIB, ]))))
u <- tumour_coef(function(D) voom(dT, D), model.matrix(~ grp), "grpTumor")
a <- tumour_coef(function(D) voom(dT, D), model.matrix(~ grp + fibT), "grpTumor")
ref <- read.delim("results/DEG/TCGA_UCEC_voom_all_genes.tsv"); rownames(ref) <- ref$gene
stopifnot(max(abs(u$logFC - ref[CORE, "logFC"])) < 1e-6)          # reproduces the primary TCGA voom analysis
zT <- t(scale(t(lcT))); pT <- -colMeans(zT[CORE, ]); isT <- grp == "Tumor"
gT <- data.frame(cohort = "TCGA-UCEC", gene = CORE, log2FC_unadjusted = round(u$logFC, 2), adjP_unadjusted = signif(u$adj.P.Val, 3),
                 log2FC_adjusted = round(a$logFC, 2), adjP_adjusted = signif(a$adj.P.Val, 3),
                 attenuation_pct = round(100 * (1 - a$logFC / u$logFC), 1))
sT <- data.frame(cohort = "TCGA-UCEC", n_tumour = sum(isT), n_normal = sum(!isT), design = "voom, ~ group (+ fibroblast score)",
                 fibroblast_score_AUC_normal_higher = round(1 - fast_auc(fibT, isT), 3), rho_panel_fibroblast = round(cor(pT, fibT, method = "spearman"), 3),
                 panel_AUC = round(fast_auc(pT, isT), 4), panel_AUC_residualized = round(fast_auc(resid(lm(pT ~ fibT)), isT), 4))

## ---------------------------------------------------------------- GSE146889 (paired)
man <- subset(read.delim("metadata/sample_manifest_GSE146889.tsv", stringsAsFactors = FALSE), tissue == "Endometrial")
raw <- read.delim("data_raw/GEO/GSE146889/GSE146889_GeneCount.tsv.gz", check.names = FALSE)
gc <- raw[, paste0(man$column, "_count")]; colnames(gc) <- man$column
sym <- raw$GeneName; ord <- order(sym, -rowMeans(gc)); idx <- ord[!duplicated(sym[ord])]
m <- as.matrix(gc[idx, ]); rownames(m) <- sym[idx]; m <- m[!is.na(rownames(m)) & rownames(m) != "", ]
dV <- DGEList(m, group = man$group); dV <- dV[filterByExpr(dV, group = man$group), , keep.lib.sizes = FALSE]; dV <- calcNormFactors(dV)
paired <- man$pair_id %in% names(which(table(man$pair_id) == 2))
dp <- dV[, paired]; f <- factor(man$group[paired], levels = c("normal","tumor")); pp <- factor(man$pair_id[paired])
lcV <- cpm(dp, log = TRUE); fibV <- colMeans(t(scale(t(lcV[FIB, ]))))
D0 <- model.matrix(~ 0 + f + pp); colnames(D0)[1:2] <- c("N","T")
D1 <- model.matrix(~ 0 + f + pp + fibV); colnames(D1)[1:2] <- c("N","T")
fitc <- function(D) { v <- voom(dp, D); tt <- topTable(eBayes(contrasts.fit(lmFit(v, D), makeContrasts(T - N, levels = D))), number = Inf); tt[CORE, ] }
u2 <- fitc(D0); a2 <- fitc(D1)
ref2 <- read.delim("results/DEG/GSE146889_all_genes.tsv"); ref2 <- ref2[!duplicated(ref2$gene), ]; rownames(ref2) <- ref2$gene
stopifnot(max(abs(u2$logFC - ref2[CORE, "logFC"])) < 1e-6)        # reproduces Table 2
zV <- t(scale(t(lcV))); pV <- -colMeans(zV[CORE, ]); isTV <- f == "tumor"
gV <- data.frame(cohort = "GSE146889", gene = CORE, log2FC_unadjusted = round(u2$logFC, 2), adjP_unadjusted = signif(u2$adj.P.Val, 3),
                 log2FC_adjusted = round(a2$logFC, 2), adjP_adjusted = signif(a2$adj.P.Val, 3),
                 attenuation_pct = round(100 * (1 - a2$logFC / u2$logFC), 1))
sV <- data.frame(cohort = "GSE146889", n_tumour = sum(isTV), n_normal = sum(!isTV), design = "voom, 35 pairs, ~ 0 + group + patient (+ fibroblast score)",
                 fibroblast_score_AUC_normal_higher = round(1 - fast_auc(fibV, isTV), 3), rho_panel_fibroblast = round(cor(pV, fibV, method = "spearman"), 3),
                 panel_AUC = round(fast_auc(pV, isTV), 4), panel_AUC_residualized = round(fast_auc(resid(lm(pV ~ fibV)), isTV), 4))

G <- rbind(gT, gV); S <- rbind(sT, sV)
print(G, row.names = FALSE); print(S, row.names = FALSE)
write.table(G, "results/revision/composition_adjusted_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(S, "results/revision/composition_adjusted_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("===== 40g v2 DONE =====\n")
