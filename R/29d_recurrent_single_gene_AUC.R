#!/usr/bin/env Rscript
# 29d_recurrent_single_gene_AUC.R -- single-gene TCGA AUC for every recurrent gene with empirical P from the exhaustive null (source data for Figure 4B)
# preprocessing and AUC definition follow R/29b_TCGA_null_controls.R exactly; the null is the 11,800 non-recurrent gene AUCs saved by 29b
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR)})
se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop=FALSE])
code <- substr(cdl$barcode,14,15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code=="01" & cdl$primary_diagnosis %in% c(EN,SR)) | code=="11"
cnt <- assay(se,"unstranded")[,keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]
cnt <- cnt[i,]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt)!="",]
grp <- factor(ifelse(code[keep]=="11","Normal","Tumor"), levels=c("Normal","Tumor"))
d <- DGEList(cnt); d <- d[filterByExpr(d,group=grp),,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
lc <- cpm(d, log=TRUE)
isT <- grp=="Tumor"; nT <- sum(isT); nN <- sum(!isT)
auc1 <- function(x){ r <- rank(x); a <- (sum(r[isT]) - nT*(nT+1)/2)/(nT*nN); max(a, 1-a) }
nd <- readRDS("results/validation/TCGA_null_distributions.rds")$null_single
rec <- read.delim("results/consensus/recurrent_consensus_annotated.tsv", stringsAsFactors=FALSE)
bs  <- read.delim("results/ML/bootstrap_stability_recurrent.tsv", stringsAsFactors=FALSE)
g <- intersect(rec$gene, rownames(lc))
res <- data.frame(gene=g, AUC=sapply(g, function(x) auc1(lc[x,])))
res$empirical_P <- sapply(res$AUC, function(a) mean(nd >= a))
res <- merge(res, rec[, c("gene","direction","tier","cell_cycle")], by="gene")
res <- merge(res, bs[, c("gene","freq_intersection")], by="gene", all.x=TRUE)
res <- res[order(-res$freq_intersection), ]
write.table(res, "results/validation/TCGA_recurrent_single_gene_AUC.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("recurrent genes %d, measurable in TCGA %d | null over %d genes (95th percentile %.3f)\n", nrow(rec), nrow(res), length(nd), quantile(nd,.95)))
cat(sprintf("recurrent genes with empirical P<0.05: %d / %d\n", sum(res$empirical_P < 0.05), nrow(res)))
cat(sprintf("Spearman rho between bootstrap frequency and AUC = %.3f (P = %.2g)\n",
    cor(res$freq_intersection, res$AUC, method="spearman"), cor.test(res$freq_intersection, res$AUC, method="spearman")$p.value))
print(head(res[, c("gene","freq_intersection","AUC","empirical_P","tier")], 8), row.names=FALSE)
print(res[res$tier=="core_3of3", c("gene","freq_intersection","AUC","empirical_P")], row.names=FALSE)
cat("===== 29d DONE =====\n")
