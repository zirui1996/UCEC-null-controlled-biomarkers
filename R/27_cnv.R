#!/usr/bin/env Rscript
# 27_cnv.R -- protocol section 22.2: dosage effect between CNV and expression
# GISTIC thresholds: -2 deep deletion / -1 deletion / 0 diploid / +1 gain / +2 amplification
# primary test: Expression_GeneX ~ CNV_GeneX (Kruskal-Wallis or Spearman)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
dir.create("results/cnv/figures", showWarnings=FALSE, recursive=TRUE)
set.seed(42); GX <- "NUMA1"

cn <- read.delim(gzfile("data_raw/Xena/TCGA_UCEC_GISTIC2_thresholded.gz"), row.names=1, check.names=FALSE)
cat(sprintf("GISTIC2: %d genes x %d samples\n", nrow(cn), ncol(cn)))
if (!GX %in% rownames(cn)) stop(GX, " is not in the GISTIC matrix")

o <- readRDS("results/immune/immune_objects.rds"); gx <- o$gx
names(gx) <- substr(names(gx),1,15)
colnames(cn) <- substr(colnames(cn),1,15)
common <- intersect(names(gx), colnames(cn))
cat("samples shared with the expression matrix:", length(common), "\n\n")

d <- data.frame(sample=common, expr=as.numeric(gx[common]), cnv=as.numeric(cn[GX, common]))
lab <- c("-2 deep deletion","-1 deletion","0 diploid","+1 gain","+2 amplification")
d$cnv_f <- factor(d$cnv, levels=-2:2, labels=lab)
cat("=== CNV distribution for NUMA1 ===\n"); print(table(d$cnv_f))

cat("\n=== section 22.2 primary test, expression ~ CNV ===\n")
kw <- kruskal.test(expr ~ cnv_f, data=d)
sp <- suppressWarnings(cor.test(d$cnv, d$expr, method="spearman"))
cat(sprintf("Kruskal-Wallis P = %.3g\n", kw$p.value))
cat(sprintf("Spearman (CNV as ordinal) rho = %.3f, P = %.3g\n", sp$estimate, sp$p.value))
st <- aggregate(expr ~ cnv_f, d, function(x) c(n=length(x), median=median(x), IQR=IQR(x)))
st <- data.frame(CNV=st$cnv_f, n=st$expr[,"n"], median=round(st$expr[,"median"],3), IQR=round(st$expr[,"IQR"],3))
print(st, row.names=FALSE)

# pairwise comparisons (only between categories with >=10 samples)
keep <- names(which(table(d$cnv_f) >= 10))
if (length(keep) >= 2) {
  pw <- pairwise.wilcox.test(d$expr[d$cnv_f %in% keep], droplevels(d$cnv_f[d$cnv_f %in% keep]), p.adjust.method="BH")
  cat("\npairwise comparisons (BH-adjusted, only categories with n>=10):\n"); print(signif(pw$p.value,3))
}

pdf("results/cnv/figures/NUMA1_CNV_dosage.pdf", width=6, height=5)
boxplot(expr ~ cnv_f, d, col="grey85", outline=FALSE, xlab="GISTIC2 copy-number call",
        ylab=sprintf("%s expression (log2 CPM)", GX), las=1,
        main=sprintf("%s dosage effect\nKruskal-Wallis P = %.3g", GX, kw$p.value))
points(jitter(as.numeric(d$cnv_f),.6), d$expr, pch=16, cex=.4, col=adjustcolor("black",.35))
dev.off()

write.table(d, "results/cnv/NUMA1_CNV_expression.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("\nverdict: %s\n", ifelse(kw$p.value < 0.05,
  "dosage effect present -> report as supporting mechanism", "no dosage effect -> report the negative result faithfully in the Supplementary (§22.2)")))
writeLines(capture.output(sessionInfo()), "logs/27_cnv_sessionInfo.txt")
cat("\n===== 27_cnv DONE =====\n")
