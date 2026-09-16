#!/usr/bin/env Rscript
# 04b_GSE106191_criteria.R -- evaluation of criteria 1 and 2 plus the S-a / S-b sensitivity analyses
# the three acceptance criteria of protocol section 4.2; criterion 3 was settled by the scan-date analysis
proj <- Sys.getenv("UCEC_PROJ", unset = getwd()); setwd(proj)
suppressPackageStartupMessages(library(limma))
set.seed(42)

man <- read.delim("metadata/sample_manifest_GSE106191.tsv", stringsAsFactors=FALSE, quote="")

fit_deg <- function(x, grp, cov=NULL) {
  f <- factor(grp, levels=c("normal","tumor"))
  if (is.null(cov) || length(unique(cov))<2) {
    d <- model.matrix(~0+f); colnames(d) <- c("N","T")
  } else {
    cv <- factor(make.names(cov))
    d <- model.matrix(~0+f+cv)
    colnames(d) <- make.names(colnames(d))
    colnames(d)[1:2] <- c("N","T")
  }
  if (qr(d)$rank < ncol(d)) {                       # fall back to a model without the covariate when it is collinear with the grouping
    d <- model.matrix(~0+f); colnames(d) <- c("N","T")
  }
  fit <- eBayes(contrasts.fit(lmFit(x, d), makeContrasts(T-N, levels=d)))
  topTable(fit, number=Inf, adjust.method="BH")
}

evaluate <- function(x, sub, cov_col, label, nperm=1000) {
  x <- x[, sub$sample_id, drop=FALSE]
  grp <- sub$normal_or_tumor
  cov <- if (!is.null(cov_col)) sub[[cov_col]] else NULL
  tt  <- fit_deg(x, grp, cov)
  obs_p001 <- sum(tt$P.Value < 0.001)
  perm <- replicate(nperm, sum(fit_deg(x, sample(grp), cov)$P.Value < 0.001))
  empP <- mean(perm >= obs_p001)
  cat(sprintf("\n===== %s =====\n", label))
  cat(sprintf("  samples: tumour %d / normal %d   covariate: %s\n",
              sum(grp=="tumor"), sum(grp=="normal"), ifelse(is.null(cov_col),"none",cov_col)))
  cat(sprintf("  criterion 1  DEGs at FDR<0.05 = %d  (|lfc|>=1: %d, |lfc|>=2: %d)  smallest adj.P = %s\n",
              sum(tt$adj.P.Val<0.05),
              sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>=1),
              sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>=2),
              signif(min(tt$adj.P.Val),3)))
  cat(sprintf("  criterion 2  observed rawP<0.001 = %d | permuted median %.0f, 95%% percentile %.0f | empirical P = %.4f\n",
              obs_p001, median(perm), quantile(perm,.95), empP))
  cat(sprintf("  verdict   criterion 1 %s   criterion 2 %s\n",
              ifelse(sum(tt$adj.P.Val<0.05)>0,"pass","fail"),
              ifelse(empP<0.05,"pass","fail")))
  tt$gene <- rownames(tt)
  list(tt=tt, empP=empP, obs=obs_p001, nsig=sum(tt$adj.P.Val<0.05))
}

x <- as.matrix(read.delim("data_processed/expression/GSE106191_gene_matrix_bestprobe.tsv",
                          check.names=FALSE, row.names=1))
p16 <- subset(man, include_primary=="yes")
res_main <- evaluate(x, p16, "scan_day", "main analysis: 2016 batch, model ~ scan_day + group")

mix <- subset(man, mixed_day=="yes")
res_sa <- evaluate(x, mix, "scan_day", "S-a: the 3 mixed scan days only (cleanest batch structure)")

xall <- as.matrix(read.delim("data_processed/expression/GSE106191_all97_gene_matrix_bestprobe.tsv",
                             check.names=FALSE, row.names=1))
res_sb <- evaluate(xall, man, "era", "S-b: all 97 samples with an era covariate (for comparison only, not the main analysis)")

dir.create("results/DEG", showWarnings=FALSE, recursive=TRUE)
for (nm in c("main","Sa","Sb")) {
  r <- get(paste0("res_", ifelse(nm=="main","main",nm)))
  tt <- r$tt[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val")]
  write.table(tt, sprintf("results/DEG/GSE106191_%s_all_genes.tsv", nm),
              sep="\t", row.names=FALSE, quote=FALSE)
  write.table(subset(tt, adj.P.Val<0.05 & abs(logFC)>=1),
              sprintf("results/DEG/GSE106191_%s_DEG_primary.tsv", nm),
              sep="\t", row.names=FALSE, quote=FALSE)
}

cat("\n\n===== agreement between the main analysis and S-a =====\n")
a <- subset(res_main$tt, adj.P.Val<0.05 & abs(logFC)>=1)$gene
b <- subset(res_sa$tt,   adj.P.Val<0.05 & abs(logFC)>=1)$gene
cat(sprintf("  main-analysis DEGs %d | S-a DEGs %d | intersection %d\n", length(a), length(b), length(intersect(a,b))))
if (length(b)) cat(sprintf("  %.1f%% of the S-a DEGs are reproduced by the main analysis\n", 100*length(intersect(a,b))/length(b)))
cat("\n===== CRITERIA DONE =====\n")
