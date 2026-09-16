#!/usr/bin/env Rscript
# 28c_signature_null_save.R -- save the random null distribution of the recurrent signature, which R/28 only printed, as a result table (same seed as R/28)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(pROC)})
set.seed(42)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(pROC)})
set.seed(42)

se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop=FALSE])
code <- substr(cdl$barcode,14,15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code=="01" & cdl$primary_diagnosis %in% c(EN,SR)) | code=="11"
cnt <- assay(se,"unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]
cnt <- cnt[i,]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt)!="", ]
grp <- factor(ifelse(code[keep]=="11","Normal","Tumor"), levels=c("Normal","Tumor"))
d <- DGEList(cnt); d <- d[filterByExpr(d, group=grp),,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
lc <- cpm(d, log=TRUE); z <- t(scale(t(lc)))
cat(sprintf("TCGA-UCEC: %d samples (Tumour %d / Normal %d) x %d genes\n\n",
    ncol(z), sum(grp=="Tumor"), sum(grp=="Normal"), nrow(z)))

rec  <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors=FALSE)
core <- read.delim("results/consensus/core_consensus_3of3.tsv", stringsAsFactors=FALSE)
ncc  <- read.delim("results/consensus/recurrent_noCellCycle.tsv", stringsAsFactors=FALSE)
lf <- setNames(rowMeans(rec[,c("GSE17025_logFC","GSE63678_logFC","GSE106191_logFC")]), rec$gene)

score <- function(g) {
  g <- intersect(g, rownames(z))
  up <- g[lf[g] > 0]; dn <- g[lf[g] < 0]
  (if (length(up)) colMeans(z[up,,drop=FALSE]) else 0) -
  (if (length(dn)) colMeans(z[dn,,drop=FALSE]) else 0) }

SETS <- list(recurrent=rec$gene, `core 4`=core$gene, noCellCycle=ncc$gene,
             recurrent_up=rec$gene[rec$direction=="up"],
             recurrent_down=rec$gene[rec$direction=="down"])
## ---- required by v2.4 section 0.7.4: null-distribution control ----
NSIG <- length(intersect(rec$gene, rownames(z)))
cat(sprintf("\n=== null-distribution control (random sets of %d genes, matching the recurrent genes measurable in TCGA; direction taken from the sign of the discovery logFC) ===\n", NSIG))
V <- c("GSE17025","GSE63678","GSE106191")
L <- lapply(V, function(g){t<-read.delim(sprintf("results/DEG/%s_all_genes.tsv",g),stringsAsFactors=FALSE)
  t<-t[!duplicated(t$gene),]; setNames(t$logFC,t$gene)})
uni <- Reduce(intersect, lapply(L,names))
mlfc <- rowMeans(sapply(L, function(x) x[uni])); names(mlfc) <- uni
pool <- intersect(setdiff(uni, rec$gene), rownames(z))
sc2 <- function(g){ up<-g[mlfc[g]>0]; dn<-g[mlfc[g]<0]
  (if(length(up)) colMeans(z[up,,drop=FALSE]) else 0) - (if(length(dn)) colMeans(z[dn,,drop=FALSE]) else 0) }
real <- as.numeric(roc(grp, score(rec$gene), quiet=TRUE)$auc)
null <- replicate(1000, as.numeric(roc(grp, sc2(sample(pool, NSIG)), quiet=TRUE)$auc))
cat(sprintf("observed recurrent signature (%d genes) AUC = %.4f\n", NSIG, real))
cat(sprintf("random signature AUC: median %.4f | 95%% interval %.4f-%.4f | maximum %.4f\n",
    median(null), quantile(null,.025), quantile(null,.975), max(null)))
cat(sprintf("empirical P = %.4f  ->  %s\n", mean(null>=real),
    ifelse(mean(null>=real)<0.05, "significantly better than random; usable as valid validation evidence",
           "indistinguishable from random; MUST NOT be used as validation evidence (as for GSE146889)")))
res <- data.frame(set="recurrent signature", n_genes=NSIG, AUC=round(real,4),
                  null_median=round(median(null),4), null_q025=round(quantile(null,.025),4),
                  null_q975=round(quantile(null,.975),4), null_max=round(max(null),4),
                  empirical_P=round(mean(null>=real),4), n_draws=1000)
write.table(res, "results/validation/TCGA_signature_null.tsv", sep="\t", row.names=FALSE, quote=FALSE)
print(res, row.names=FALSE)
cat("===== 28c DONE =====\n")
