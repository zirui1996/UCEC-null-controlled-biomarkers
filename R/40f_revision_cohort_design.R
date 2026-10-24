#!/usr/bin/env Rscript
# 40f_revision_cohort_design.R
#   R1 minor 10  full results of the rejected four-cohort (>=3 of 4, GSE146889 voting) design, without the withdrawn veto
#   R3 major 4   how unequal cohort sizes shape consensus discovery: GSE17025 (91 EC / 12 controls) is downsampled to the
#                size of the other cohorts and the consensus is rebuilt, holding the other two cohorts' calls fixed
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(limma)})
set.seed(42); FDR <- 0.05; LFC <- 1; NREP <- 200; NNULL <- 1000

load1 <- function(g) { t <- read.delim(sprintf("results/DEG/%s_all_genes.tsv", g), stringsAsFactors = FALSE); t <- t[!duplicated(t$gene), ]; rownames(t) <- t$gene; t }
rec  <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors = FALSE)$gene
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")
consensus <- function(tabs, k) {
  u <- Reduce(intersect, lapply(tabs, rownames))
  sig <- sapply(tabs, function(t) t[u, "adj.P.Val"] < FDR & abs(t[u, "logFC"]) >= LFC)
  lfc <- sapply(tabs, function(t) t[u, "logFC"]); rownames(sig) <- rownames(lfc) <- u
  ok <- vapply(seq_along(u), function(i) { s <- which(sig[i, ]); length(s) >= k && length(unique(sign(lfc[i, s]))) == 1 }, TRUE)
  list(genes = u[ok], universe = u, sig = sig, lfc = lfc)
}

## ---- R1 minor 10: four-cohort design, no veto
T3 <- lapply(c("GSE17025","GSE63678","GSE106191"), load1)
T4 <- c(T3, list(load1("GSE146889")))
c3 <- consensus(T3, 2); c4 <- consensus(T4, 3); c4core <- consensus(T4, 4)
stopifnot(setequal(c3$genes, rec))
cmb <- table(apply(c4$sig[c4$genes, , drop = FALSE], 1, function(x) paste(c("17025","63678","106191","146889")[x], collapse = "+")))
cat(sprintf("three-cohort recurrent %d | four-cohort >=3/4: %d (universe %d) | 4/4: %d (%s)\n", length(c3$genes), length(c4$genes),
            length(c4$universe), length(c4core$genes), paste(c4core$genes, collapse = ",")))
cat(sprintf("four-cohort set within the three-cohort set: %d kept, %d new\n", length(intersect(c4$genes, c3$genes)), length(setdiff(c4$genes, c3$genes))))
print(sort(cmb, decreasing = TRUE))

# discrimination of the four-cohort signature in TCGA (GSE146889 helped define it and cannot evaluate it)
se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop = FALSE]); code <- substr(cdl$barcode, 14, 15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code == "01" & cdl$primary_diagnosis %in% c(EN, SR)) | code == "11"
cnt <- assay(se, "unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]; cnt <- cnt[i, ]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt) != "", ]
grp <- factor(ifelse(code[keep] == "11", "N", "T")); isT <- grp == "T"
d <- DGEList(cnt); d <- d[filterByExpr(d, group = grp), , keep.lib.sizes = FALSE]; d <- calcNormFactors(d); z <- t(scale(t(cpm(d, log = TRUE))))
mlfc <- rowMeans(c3$lfc); fast_auc <- function(x) { r <- rank(x); nT <- sum(isT); nN <- sum(!isT); (sum(r[isT]) - nT * (nT + 1) / 2) / (nT * nN) }
sc <- function(gg) { gg <- intersect(gg, rownames(z)); up <- gg[mlfc[gg] > 0]; dn <- gg[mlfc[gg] < 0]
  (if (length(up)) colMeans(z[up, , drop = FALSE]) else 0) - (if (length(dn)) colMeans(z[dn, , drop = FALSE]) else 0) }
g4 <- intersect(c4$genes, rownames(z)); pool <- intersect(setdiff(c3$universe, rec), rownames(z))
a4 <- fast_auc(sc(g4)); n4 <- replicate(NNULL, fast_auc(sc(sample(pool, length(g4)))))
four <- data.frame(design = c("three-cohort recurrent (>=2/3, reported)", "four-cohort recurrent (>=3/4, rejected)"),
                   n_genes = c(length(c3$genes), length(c4$genes)), universe = c(length(c3$universe), length(c4$universe)),
                   n_TCGA = c(length(intersect(rec, rownames(z))), length(g4)),
                   TCGA_AUC = round(c(fast_auc(sc(intersect(rec, rownames(z)))), a4), 4),
                   null_median = round(c(NA, median(n4)), 4), empirical_P = c(0.193, (sum(n4 >= a4) + 1) / (NNULL + 1)),
                   median_abs_log2FC = round(c(median(abs(mlfc[c3$genes])), median(abs(rowMeans(c4$lfc[c4$genes, ])))), 2))
print(four, row.names = FALSE)

## ---- R3 major 4: downsample GSE17025 to 12 EC + 12 controls (GSE63678 has 7 + 5, GSE106191 64 + 33)
x <- read.delim("data_processed/expression/GSE17025_gene_matrix_bestprobe.tsv", row.names = 1, check.names = FALSE)
man <- subset(read.delim("metadata/sample_manifest.tsv", stringsAsFactors = FALSE), dataset == "GSE17025" & include_primary == "yes")
x <- as.matrix(x[, man$sample_id]); tum <- man$sample_id[man$normal_or_tumor == "tumor"]; nor <- man$sample_id[man$normal_or_tumor != "tumor"]
deg <- function(cols) { f <- factor(ifelse(cols %in% tum, "T", "N"), levels = c("N","T")); dsg <- model.matrix(~ f)
  tt <- topTable(eBayes(lmFit(x[, cols], dsg)), coef = 2, number = Inf); tt[, c("logFC","adj.P.Val")] }
down <- function(nt) t(replicate(NREP, {
  tt <- deg(c(sample(tum, nt), nor)); cs <- consensus(list(tt, T3[[2]], T3[[3]]), 2)
  c(n_DE_GSE17025 = sum(tt$adj.P.Val < FDR & abs(tt$logFC) >= LFC), n_recurrent = length(cs$genes),
    retained_of_182 = length(intersect(cs$genes, rec)), core4_retained = sum(CORE %in% consensus(list(tt, T3[[2]], T3[[3]]), 3)$genes))
}))
ds <- lapply(c(12, 30), down); names(ds) <- c("12 EC + 12 controls", "30 EC + 12 controls")
ds_tab <- do.call(rbind, lapply(names(ds), function(n) data.frame(GSE17025_size = n,
  t(apply(ds[[n]], 2, function(v) sprintf("%.0f [%.0f-%.0f]", median(v), quantile(v, .025), quantile(v, .975)))), check.names = FALSE)))
full <- data.frame(GSE17025_size = "full (91 EC + 12 controls)", n_DE_GSE17025 = sum(T3[[1]]$adj.P.Val < FDR & abs(T3[[1]]$logFC) >= LFC),
                   n_recurrent = length(rec), retained_of_182 = length(rec), core4_retained = 4, check.names = FALSE)
ds_tab <- rbind(full, ds_tab); print(ds_tab, row.names = FALSE)
for (n in names(ds)) cat(sprintf("%s: all four core genes retained in %d/%d downsamples; per gene: %s\n", n,
  sum(ds[[n]][, "core4_retained"] == 4), NREP, paste(table(factor(ds[[n]][, "core4_retained"], levels = 0:4)), collapse = "/")))
write.table(four, "results/revision/four_cohort_design.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(as.data.frame(cmb), "results/revision/four_cohort_combinations.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ds_tab, "results/revision/GSE17025_downsampling.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("===== 40f DONE =====\n")
