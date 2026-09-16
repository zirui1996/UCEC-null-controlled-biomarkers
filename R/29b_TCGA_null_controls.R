#!/usr/bin/env Rscript
# 29b_TCGA_null_controls.R -- null distributions for discrimination in TCGA-UCEC (written to disk in v2.6)
#   single gene: AUC computed for every gene in the random pool (exhaustive, no Monte Carlo error); AUC taken as max(AUC, 1-AUC), the observed value on the same basis
#   core 4 panel: 10,000 random four-gene panels (2,000 in v2.5; the P value varied between 0.0025 and 0.0070 as the random pool was adjusted)
#   random pool = genes shared by the three discovery cohorts, outside the recurrent consensus, and measurable in TCGA; panel direction is the sign of the mean discovery logFC
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(pROC)})
set.seed(42); CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1"); NP <- 10000
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
V <- c("GSE17025","GSE63678","GSE106191")
L <- lapply(V, function(x){t<-read.delim(sprintf("results/DEG/%s_all_genes.tsv",x),stringsAsFactors=FALSE)
  t<-t[!duplicated(t$gene),]; setNames(t$logFC,t$gene)})
uni <- Reduce(intersect, lapply(L,names)); mlfc <- rowMeans(sapply(L,function(x)x[uni])); names(mlfc) <- uni
rec <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors=FALSE)$gene
pool <- intersect(setdiff(uni, rec), rownames(z))
cat(sprintf("TCGA: Tumour %d / Normal %d | random pool %d genes\n", sum(grp=="Tumor"), sum(grp=="Normal"), length(pool)))

## fast AUC (Mann-Whitney rank formula, identical to pROC)
nT <- sum(grp=="Tumor"); nN <- sum(grp=="Normal"); isT <- grp=="Tumor"
auc1 <- function(x){ r <- rank(x); a <- (sum(r[isT]) - nT*(nT+1)/2)/(nT*nN); max(a, 1-a) }
null1 <- apply(lc[pool,], 1, auc1)
cat(sprintf("\nsingle-gene null distribution (exhaustive over %d genes): median %.3f | Q3 %.3f | 95th percentile %.3f\n",
    length(null1), median(null1), quantile(null1,.75), quantile(null1,.95)))
chk <- as.numeric(roc(grp, lc["NUMA1",], quiet=TRUE)$auc); stopifnot(abs(chk - auc1(lc["NUMA1",])) < 1e-8)

GENES <- c(CORE, "NUMA1", "UCK2", "NAALAD2")
res1 <- do.call(rbind, lapply(GENES, function(g) { a <- auc1(lc[g,])
  data.frame(feature=g, n_genes=1, AUC=a, null_median=median(null1), null_q95=quantile(null1,.95),
             empirical_P=mean(null1 >= a), null="all pool genes (exhaustive)") }))

sc <- function(gg){ up<-gg[mlfc[gg]>0]; dn<-gg[mlfc[gg]<0]
  (if(length(up)) colMeans(z[up,,drop=FALSE]) else 0) - (if(length(dn)) colMeans(z[dn,,drop=FALSE]) else 0) }
obsP <- auc1(-colMeans(z[CORE,,drop=FALSE]))
nullP <- replicate(NP, auc1(sc(sample(pool, 4))))
resP <- data.frame(feature="core 4 panel", n_genes=4, AUC=obsP, null_median=median(nullP), null_q95=quantile(nullP,.95),
                   empirical_P=(sum(nullP >= obsP) + 1)/(NP + 1), null=sprintf("random 4-gene panels (%s)", format(NP, big.mark=",")))
res <- rbind(res1, resP); rownames(res) <- NULL
cat(sprintf("panel null distribution (%d draws): median %.3f | 95%% interval %.3f-%.3f | %d draws at or above the observed value\n",
    NP, median(nullP), quantile(nullP,.025), quantile(nullP,.975), sum(nullP >= obsP)))
res$AUC <- round(res$AUC,4); res$null_median <- round(res$null_median,3); res$null_q95 <- round(res$null_q95,3)
res$empirical_P <- signif(res$empirical_P, 3)
print(res, row.names=FALSE)
write.table(res, "results/validation/TCGA_null_controls.tsv", sep="\t", row.names=FALSE, quote=FALSE)
saveRDS(list(null_single=null1, null_panel=nullP, obs_panel=obsP), "results/validation/TCGA_null_distributions.rds")
cat("===== 29b DONE =====\n")
