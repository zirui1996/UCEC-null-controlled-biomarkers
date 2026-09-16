#!/usr/bin/env Rscript
# 29_core4_panel.R -- full evaluation of the four-gene core consensus panel
# the core 4 definition (significant in 3/3 cohorts, same direction, not vetoed by GSE56087) was fixed in v2.3 section 0.6.1, before any validation
# every discrimination metric is reported alongside a null distribution, as required by v2.4 section 0.7.4
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(pROC); library(survival)})
set.seed(42)
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")

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
lc <- cpm(d, log=TRUE); z <- t(scale(t(lc)))
cat(sprintf("TCGA: Tumor %d / Normal %d\n\n", sum(grp=="Tumor"), sum(grp=="Normal")))

## ---- 1. panel score (four downregulated genes, negative mean) + null distribution ----
g <- intersect(CORE, rownames(z))
s <- -colMeans(z[g,,drop=FALSE])
r <- roc(grp, s, quiet=TRUE); ci <- as.numeric(ci.auc(r))
cat(sprintf("=== core 4 panel（TCGA）===\nAUC = %.4f (95%%CI %.4f-%.4f) | Wilcoxon P = %.3g\n",
    as.numeric(r$auc), ci[1], ci[3], wilcox.test(s~grp)$p.value))
V <- c("GSE17025","GSE63678","GSE106191")
L <- lapply(V, function(x){t<-read.delim(sprintf("results/DEG/%s_all_genes.tsv",x),stringsAsFactors=FALSE)
  t<-t[!duplicated(t$gene),]; setNames(t$logFC,t$gene)})
uni <- Reduce(intersect, lapply(L,names)); mlfc <- rowMeans(sapply(L,function(x)x[uni])); names(mlfc)<-uni
rec <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors=FALSE)$gene
pool <- intersect(setdiff(uni, rec), rownames(z))
sc <- function(gg){ up<-gg[mlfc[gg]>0]; dn<-gg[mlfc[gg]<0]
  (if(length(up)) colMeans(z[up,,drop=FALSE]) else 0) - (if(length(dn)) colMeans(z[dn,,drop=FALSE]) else 0) }
null <- replicate(2000, as.numeric(roc(grp, sc(sample(pool,4)), quiet=TRUE)$auc))
null <- pmax(null, 1-null)
cat(sprintf("null distribution (2000 random four-gene panels): median %.3f | 95%% interval %.3f-%.3f | maximum %.3f\n",
    median(null), quantile(null,.025), quantile(null,.975), max(null)))
cat(sprintf("empirical P = %.4f  ->  %s\n\n", mean(null>=as.numeric(r$auc)),
    ifelse(mean(null>=as.numeric(r$auc))<0.05,"SIGNIFICANTLY BETTER THAN RANDOM; usable as validation evidence","indistinguishable from random; must not be used as evidence")))

## ---- 2. per gene (TCGA + matched pairs) ----
cat("=== per gene (TCGA) ===\n")
pid <- substr(cdl$barcode[keep],1,12)
pr <- intersect(pid[grp=="Tumor"], pid[grp=="Normal"])
for (x in CORE) {
  if (!x %in% rownames(lc)) next
  a <- roc(grp, lc[x,], quiet=TRUE); ax <- as.numeric(a$auc)
  pt <- sapply(pr, function(p) mean(lc[x, pid==p & grp=="Tumor"]))
  pn <- sapply(pr, function(p) mean(lc[x, pid==p & grp=="Normal"]))
  cat(sprintf("  %-11s Normal %6.2f -> Tumor %6.2f | Wilcoxon %.2e | AUC %.3f | paired (n=%d) P=%.2e\n",
      x, median(lc[x,grp=="Normal"]), median(lc[x,grp=="Tumor"]),
      wilcox.test(lc[x,]~grp)$p.value, ax, length(pr), wilcox.test(pt,pn,paired=TRUE)$p.value))
}

## ---- 3. independent validation in GSE146889 ----
cat("\n=== GSE146889 (independent, never used for selection) ===\n")
v <- read.delim("results/DEG/GSE146889_all_genes.tsv", stringsAsFactors=FALSE); rownames(v)<-v$gene
for (x in CORE) {
  if (x %in% rownames(v)) {
    cat(sprintf("  %-11s logFC %+.2f | adj.P %.2e | direction %s\n", x, v[x,"logFC"], v[x,"adj.P.Val"],
        ifelse(v[x,"logFC"]<0,"concordant (down)","CONFLICT")))
  } else cat(sprintf("  %-11s NOT MEASURABLE\n", x))
}

## ---- 4. prognosis (rerun for the new primary gene set, and labelled as such) ----
cat("\n=== prognosis (overall survival, TCGA; n and events below) ===\n")
sv <- readRDS("results/survival/survival_objects.rds")$d
tz <- z[, grp=="Tumor"]; colnames(tz) <- substr(cdl$barcode[keep][grp=="Tumor"],1,12)
sv$pid <- substr(rownames(sv), 1, 12)
if (is.null(rownames(sv)) || all(is.na(sv$pid))) sv$pid <- substr(cdl$barcode[keep][grp=="Tumor"],1,12)[seq_len(nrow(sv))]
m <- match(sv$pid, colnames(tz))
sv$panel <- -colMeans(tz[g, m, drop=FALSE])
for (x in c(CORE,"panel")) {
  vv <- if (x=="panel") sv$panel else as.numeric(tz[x, m])
  if (all(is.na(vv))) next
  sv$tmp <- as.numeric(scale(vv))
  u <- coxph(Surv(OS_time,OS_event) ~ tmp, data=sv); ss <- summary(u)
  cat(sprintf("  %-11s HR = %.3f (%.3f-%.3f), P = %.3f | C-index %.3f\n", x,
      ss$coef[,"exp(coef)"], ss$conf.int[,3], ss$conf.int[,4], ss$coef[,"Pr(>|z|)"], ss$concordance[1]))
}
cat(sprintf("  (n = %d, events = %d)\n", sum(!is.na(sv$panel)), sum(sv$OS_event[!is.na(sv$panel)])))
cat("\n===== 29_core4_panel DONE =====\n")
