#!/usr/bin/env Rscript
# 28_signature_AUC.R -- discrimination of the recurrent consensus signature in TCGA-UCEC (v2.6: 182 genes)
# validity: the recurrent genes were selected entirely from the three GEO cohorts; TCGA played no part in selection (independence statement of section 15)
# per the rule of v2.4 section 0.7.4: a null distribution from random gene sets of the same size must be reported alongside
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
cat("=== signature discrimination (tumour vs solid tissue normal) ===\n")
res <- data.frame()
for (nm in names(SETS)) {
  g <- intersect(SETS[[nm]], rownames(z))
  s <- score(SETS[[nm]]); r <- roc(grp, s, quiet=TRUE); ci <- as.numeric(ci.auc(r))
  w <- wilcox.test(s ~ grp)
  cat(sprintf("  %-18s n=%3d genes  AUC = %.4f (95%%CI %.4f-%.4f)  Wilcoxon P = %.3g\n",
      nm, length(g), as.numeric(r$auc), ci[1], ci[3], w$p.value))
  res <- rbind(res, data.frame(set=nm, n_genes=length(g), AUC=round(as.numeric(r$auc),4),
    CI_low=round(ci[1],4), CI_high=round(ci[3],4), P=signif(w$p.value,3))) }

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
cat(sprintf("\nreference: fraction of genes at FDR<0.05 genome-wide in TCGA = %.1f%% (70.7%% in GSE146889)\n",
    { des <- model.matrix(~grp); v <- limma::voom(d, des)
      f <- limma::eBayes(limma::lmFit(v, des)); 100*mean(limma::topTable(f, coef=2, number=Inf)$adj.P.Val<0.05) }))

pdf("results/validation/figures/TCGA_signature_ROC.pdf", width=5.5, height=5.5)
r <- roc(grp, score(rec$gene), quiet=TRUE)
plot(r, col="#B2182B", lwd=2.5, main=sprintf("%d-gene recurrent consensus signature\nTCGA-UCEC  AUC = %.3f", nrow(rec), as.numeric(r$auc)))
r2 <- roc(grp, -colMeans(z[intersect(core$gene, rownames(z)),,drop=FALSE]), quiet=TRUE)
plot(r2, col="#2166AC", lwd=2, add=TRUE)
r3 <- roc(grp, z["NUMA1",], quiet=TRUE)
plot(r3, col="#4D4D4D", lwd=2, lty=2, add=TRUE)
legend("bottomright", bty="n", lwd=2, col=c("#B2182B","#2166AC","#4D4D4D"), lty=c(1,1,2),
  legend=sprintf(c(paste0(nrow(rec), "-gene signature %.3f"),"core 4 genes %.3f","NUMA1 alone %.3f"),
                 c(as.numeric(r$auc), as.numeric(r2$auc), as.numeric(r3$auc))))
dev.off()
write.table(res, "results/validation/TCGA_signature_AUC.tsv", sep="\t", row.names=FALSE, quote=FALSE)
writeLines(capture.output(sessionInfo()), "logs/28_signature_sessionInfo.txt")
cat("\n===== 28_signature_AUC DONE =====\n")
