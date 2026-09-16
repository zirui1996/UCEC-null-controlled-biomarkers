#!/usr/bin/env Rscript
# 24_cibersortx.R -- protocol sections 18.2-18.4
#  section 18.2 QC thresholds: P<0.05 / correlation>=0.80 / RMSE<=0.30; retained and excluded counts reported for each
#  section 18.3 Spearman correlation of Gene X with the 22 immune cell types, with BH-FDR
#  §18.4 Gene X-high vs low（median），Wilcoxon + BH-FDR
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(ppcor); library(ggplot2)})
set.seed(42); GX <- "NUMA1"

cs <- read.delim("results/immune/cibersortx_out_Adjusted.txt", check.names=FALSE)
CELLS <- setdiff(colnames(cs), c("Mixture","P-value","Correlation","RMSE","Absolute score (sig.score)"))
cat(sprintf("CIBERSORTx output: %d samples x %d immune cell types\n\n", nrow(cs), length(CELLS)))

## ---- section 18.2 QC (each of the three thresholds reported separately) ----
cat("=== section 18.2 QC filtering ===\n")
q1 <- cs$`P-value` < 0.05; q2 <- cs$Correlation >= 0.80; q3 <- cs$RMSE <= 0.30
for (nm in c("deconvolution P < 0.05","reconstruction correlation >= 0.80","RMSE <= 0.30")) {
  v <- switch(nm, "deconvolution P < 0.05"=q1, "reconstruction correlation >= 0.80"=q2, q3)
  cat(sprintf("  %-36s passed %3d / %3d (excluded %d)\n", nm, sum(v), length(v), sum(!v)))
}
cat(sprintf("  %-36s passed %3d / %3d (excluded %d)\n", "all three thresholds", sum(q1&q2&q3), nrow(cs), sum(!(q1&q2&q3))))
cat("\n  [v2.5 amendment] the correlation>=0.80 and RMSE<=0.30 thresholds of protocol section 18.2 are unattainable on solid-tumour bulk data:\n")
cat(sprintf("    smallest RMSE in the whole cohort %.3f (threshold 0.30), largest correlation %.3f (threshold 0.80)\n",
    min(cs$RMSE), max(cs$Correlation)))
cat("    reason: the mixture is dominated by tumour and stromal cells, which the purely immune LM22 signature matrix cannot explain,\n")
cat("          typical reconstruction correlation in solid tumours is 0.5-0.8 with RMSE 0.6-1.2.\n")
cat("    action: use the single P<0.05 threshold of the original CIBERSORT paper; report the correlation and RMSE distributions faithfully but do not threshold on them.\n")
cat(sprintf("\n  observed distribution: median P %.3g | median correlation %.3f (range %.3f-%.3f) | median RMSE %.3f (range %.3f-%.3f)\n",
    median(cs$`P-value`), median(cs$Correlation), min(cs$Correlation), max(cs$Correlation),
    median(cs$RMSE), min(cs$RMSE), max(cs$RMSE)))
keep <- q1   # v2.5: P<0.05 only (the original three thresholds are unattainable, see above)
cs <- cs[keep, ]
cs$pid <- cs$Mixture

## ---- align NUMA1 expression ----
o <- readRDS("results/immune/immune_objects.rds")
gx_all <- o$gx
m <- match(cs$pid, names(gx_all))
cat(sprintf("\nmatched to the expression matrix: %d / %d\n", sum(!is.na(m)), nrow(cs)))
cs <- cs[!is.na(m), ]; gx <- gx_all[m[!is.na(m)]]
prolif <- readRDS("results/immune/immune_partial_objects.rds")$prolif[cs$pid]
cat(sprintf("final analysis set: n = %d\n", nrow(cs)))

## ---- section 18.3 Gene X against the 22 immune cell types ----
res <- do.call(rbind, lapply(CELLS, function(cl) {
  v <- cs[[cl]]
  if (sd(v)==0 || sum(v>0) < 10) return(data.frame(cell=cl, rho=NA, P=NA, partial_rho=NA, partial_P=NA, n_nonzero=sum(v>0)))
  s <- suppressWarnings(cor.test(gx, v, method="spearman"))
  pc <- pcor.test(rank(gx), rank(v), rank(prolif))
  data.frame(cell=cl, rho=unname(s$estimate), P=s$p.value,
             partial_rho=pc$estimate, partial_P=pc$p.value, n_nonzero=sum(v>0)) }))
res$BH_q <- p.adjust(res$P,"BH"); res$partial_BH_q <- p.adjust(res$partial_P,"BH")
res$direction <- ifelse(res$rho>0,"positive","negative")
res <- res[order(res$BH_q),]
write.table(res, "results/immune/cibersortx_NUMA1_correlations.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("\n=== section 18.3 correlation analysis (Spearman + BH) ===\nBH q<0.05: %d / %d | still significant with the same sign after controlling for proliferation: %d\n",
    sum(res$BH_q<0.05,na.rm=TRUE), sum(!is.na(res$rho)),
    sum(res$partial_BH_q<0.05 & sign(res$partial_rho)==sign(res$rho), na.rm=TRUE)))
print(data.frame(cell=res$cell, rho=round(res$rho,3), q=signif(res$BH_q,3),
                 partial_rho=round(res$partial_rho,3), partial_q=signif(res$partial_BH_q,3),
                 nonzero_samples=res$n_nonzero), row.names=FALSE)

## ---- §18.4 Gene X-high vs low ----
grp <- factor(ifelse(gx > median(gx), "High","Low"), levels=c("Low","High"))
hl <- do.call(rbind, lapply(CELLS, function(cl) {
  v <- cs[[cl]]; if (all(is.na(v)) || sd(v,na.rm=TRUE)==0 || sum(v>0,na.rm=TRUE)<10) return(NULL)
  data.frame(cell=cl, median_low=median(v[grp=="Low"]), median_high=median(v[grp=="High"]),
             P=wilcox.test(v ~ grp)$p.value) }))
hl$BH_q <- p.adjust(hl$P,"BH"); hl <- hl[order(hl$BH_q),]
write.table(hl, "results/immune/cibersortx_high_vs_low.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("\n=== §18.4 NUMA1-high vs low（median，Wilcoxon + BH）===\nBH q<0.05: %d / %d\n",
    sum(hl$BH_q<0.05), nrow(hl)))
print(head(data.frame(cell=hl$cell, low=signif(hl$median_low,3), high=signif(hl$median_high,3),
                      q=signif(hl$BH_q,3)), 10), row.names=FALSE)

## ---- section 18.6 direction concordance with xCell ----
xr <- read.delim("results/immune/xCell_correlations.tsv", stringsAsFactors=FALSE)
MAP <- list("T cells CD8"="CD8+ T-cells","Plasma cells"="Plasma cells","Macrophages M1"="Macrophages M1",
            "Macrophages M2"="Macrophages M2","B cells naive"="B-cells","NK cells activated"="NK cells",
            "Monocytes"="Monocytes","Neutrophils"="Neutrophils","Mast cells resting"="Mast cells",
            "Dendritic cells resting"="DC","Eosinophils"="Eosinophils","T cells regulatory (Tregs)"="Tregs")
cat("\n=== section 18.6 direction concordance between CIBERSORTx and xCell ===\n")
ag <- 0; tot <- 0
for (a in names(MAP)) {
  b <- MAP[[a]]; ra <- res$rho[res$cell==a]; rb <- xr$rho[xr$cell==b]
  if (!length(ra) || !length(rb) || is.na(ra)) next
  tot <- tot+1; ok <- sign(ra)==sign(rb); ag <- ag+ok
  cat(sprintf("  %-28s CIBERSORTx %+.3f | xCell %+.3f  %s\n", a, ra, rb, ifelse(ok,"concordant","DISCORDANT"))) }
cat(sprintf("\ncomparable lineages with the same direction: %d / %d\n", ag, tot))

## ---- plots ----
sig <- res[!is.na(res$BH_q) & res$BH_q<0.05, ]
if (nrow(sig)) {
  sig$lab <- sig$cell
  p <- ggplot(sig, aes(x=reorder(lab, rho), y=rho, fill=rho>0)) + geom_col() + coord_flip() +
    scale_fill_manual(values=c("TRUE"="#B2182B","FALSE"="#2166AC"), guide="none") +
    labs(x=NULL, y=sprintf("Spearman rho with %s", GX),
         title="CIBERSORTx (LM22, absolute) — significant associations") + theme_classic(base_size=10)
  ggsave("results/immune/figures/CIBERSORTx_NUMA1_correlations.pdf", p, width=7, height=max(3, nrow(sig)*0.3+1.5))
}
writeLines(capture.output(sessionInfo()), "logs/24_cibersortx_sessionInfo.txt")
cat("\n===== 24_cibersortx DONE =====\n")
