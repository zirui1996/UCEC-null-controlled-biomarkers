#!/usr/bin/env Rscript
# 08_validate_GSE146889.R -- GSE146889 as an independent validation cohort
# validity: GSE146889 played no part in the consensus or in any feature selection, so it is genuinely held out for the recurrent consensus genes.
# four things are checked: direction concordance, effect-size correlation, independent replication of significance, and sample-level discrimination (AUC)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(limma); library(edgeR); library(pROC)})
set.seed(42)
dir.create("results/validation/figures", showWarnings=FALSE, recursive=TRUE)

rec  <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors=FALSE)
core <- read.delim("results/consensus/core_consensus_3of3.tsv", stringsAsFactors=FALSE)
v    <- read.delim("results/DEG/GSE146889_all_genes.tsv", stringsAsFactors=FALSE)   # paired analysis
rownames(v) <- v$gene
disc_lfc <- rowMeans(rec[, c("GSE17025_logFC","GSE63678_logFC","GSE106191_logFC")])
names(disc_lfc) <- rec$gene

ok <- intersect(rec$gene, rownames(v))
cat(sprintf("of the %d recurrent genes, %d are measurable in GSE146889\n\n", nrow(rec), length(ok)))

## ---- 1. direction concordance ----
same <- sign(disc_lfc[ok]) == sign(v[ok,"logFC"])
bt <- binom.test(sum(same), length(ok), 0.5)
cat(sprintf("[1] direction concordance %d/%d = %.1f%%  (binomial P = %.3g)\n",
            sum(same), length(ok), 100*mean(same), bt$p.value))
cat(sprintf("     concordant among upregulated genes %d/%d | among downregulated genes %d/%d\n",
    sum(same & disc_lfc[ok]>0), sum(disc_lfc[ok]>0),
    sum(same & disc_lfc[ok]<0), sum(disc_lfc[ok]<0)))

## ---- 2. effect-size correlation ----
ct <- cor.test(disc_lfc[ok], v[ok,"logFC"], method="spearman")
cp <- cor.test(disc_lfc[ok], v[ok,"logFC"], method="pearson")
cat(sprintf("\n[2] effect-size correlation  Spearman rho = %.3f (P = %.3g) | Pearson r = %.3f\n",
            ct$estimate, ct$p.value, cp$estimate))

## ---- 3. independent replication of significance ----
for (thr in c(0.05, 0.01)) {
  n <- sum(v[ok,"adj.P.Val"] < thr)
  nd <- sum(v[ok,"adj.P.Val"] < thr & same)
  cat(sprintf("\n[3] replication at FDR<%.2f: %d/%d = %.1f%% (of which %d also share the direction, %.1f%%)\n",
      thr, n, length(ok), 100*n/length(ok), nd, 100*nd/length(ok)))
}
# background control: draw the same number of genes at random and measure the replication rate
bgset <- setdiff(rownames(v), rec$gene)
nullrep <- replicate(1000, mean(v[sample(bgset, length(ok)), "adj.P.Val"] < 0.05))
cat(sprintf("     fraction at FDR<0.05 for random gene sets of the same size: median %.1f%% (95%% interval %.1f-%.1f%%)\n",
    100*median(nullrep), 100*quantile(nullrep,.025), 100*quantile(nullrep,.975)))

## ---- 4. sample-level discrimination ----
man <- subset(read.delim("metadata/sample_manifest_GSE146889.tsv", stringsAsFactors=FALSE), tissue=="Endometrial")
raw <- read.delim("data_raw/GEO/GSE146889/GSE146889_GeneCount.tsv.gz", check.names=FALSE)
cnt <- raw[, paste0(man$column,"_count")]; colnames(cnt) <- man$column
sym <- raw$GeneName
ord <- order(sym, -rowMeans(cnt)); idx <- ord[!duplicated(sym[ord])]
m <- as.matrix(cnt[idx,]); rownames(m) <- sym[idx]; m <- m[!is.na(rownames(m)) & rownames(m)!="",]
d <- DGEList(m); d <- d[filterByExpr(d, group=man$group),,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
lc <- cpm(d, log=TRUE)
z  <- t(scale(t(lc)))

score <- function(g) {
  up <- intersect(g[disc_lfc[g]>0], rownames(z)); dn <- intersect(g[disc_lfc[g]<0], rownames(z))
  colMeans(z[up,,drop=FALSE]) - colMeans(z[dn,,drop=FALSE])
}
lab <- factor(man$group, levels=c("normal","tumor"))
for (nm in c("recurrent","core_4")) {
  g <- if (nm=="recurrent") rec$gene else core$gene
  g <- intersect(g, rownames(z))
  s <- if (nm=="core_4") -colMeans(z[g,,drop=FALSE]) else score(g)   # the core genes are all downregulated, so the sign is flipped
  r <- roc(lab, s, quiet=TRUE)
  w <- wilcox.test(s ~ lab)
  cat(sprintf("\n[4] %s (%d genes)  AUC = %.3f (95%%CI %.3f-%.3f) | Wilcoxon P = %.3g\n",
      nm, length(g), as.numeric(r$auc), as.numeric(ci.auc(r))[1], as.numeric(ci.auc(r))[3], w$p.value))
  pdf(sprintf("results/validation/figures/ROC_GSE146889_%s.pdf", nm), width=5, height=5)
  plot(r, main=sprintf("GSE146889 (independent)\n%s  AUC=%.3f", nm, as.numeric(r$auc)), col="#B2182B", lwd=2)
  dev.off()
  # within-pair comparison
  pr <- man$pair_id %in% names(which(table(man$pair_id)==2))
  dif <- tapply(s[pr], man$pair_id[pr], function(x) diff(x[order(man$group[pr][man$pair_id[pr]==names(which(table(man$pair_id)==2))[1]])]))
  cat(sprintf("     pairs with score(tumor)>score(normal) among the 35 matched pairs: %d/35\n",
      sum(sapply(split(data.frame(s=s[pr], g=man$group[pr]), man$pair_id[pr]),
                 function(x) x$s[x$g=="tumor"] > x$s[x$g=="normal"]))))
}

## ---- 5. per-gene validation table ----
out <- data.frame(gene=ok, discovery_meanLFC=round(disc_lfc[ok],3),
                  GSE146889_logFC=round(v[ok,"logFC"],3), GSE146889_adjP=signif(v[ok,"adj.P.Val"],3),
                  direction_match=same, replicated_FDR05=v[ok,"adj.P.Val"]<0.05,
                  tier=ifelse(ok %in% core$gene,"core_3of3","recurrent_2of3"))
out <- out[order(out$GSE146889_adjP),]
write.table(out, "results/validation/GSE146889_validation_recurrent.tsv", sep="\t", row.names=FALSE, quote=FALSE)
writeLines(capture.output(sessionInfo()), "logs/08_validate_sessionInfo.txt")
cat("\n===== 08_validate_GSE146889 DONE =====\n")
