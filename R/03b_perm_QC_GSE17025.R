#!/usr/bin/env Rscript
# 03b_perm_QC_GSE17025.R -- completes the section 8 admission permutation test (GSE17025)
# the statistics follow exactly the same definition as GSE63678 / GSE106191 / GSE56087 / GSE115810:
#   number of genes with raw P<0.001 after limma; null = 1000 label permutations with group sizes preserved
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages(library(limma))
set.seed(42)

x <- read.delim("data_processed/expression/GSE17025_gene_matrix_bestprobe.tsv",
                row.names=1, check.names=FALSE)
man <- read.delim("metadata/sample_manifest.tsv", stringsAsFactors=FALSE)
man <- subset(man, dataset=="GSE17025" & include_primary=="yes")
stopifnot(setequal(man$sample_id, colnames(x)))
x <- x[, man$sample_id]
grp <- ifelse(man$normal_or_tumor=="tumor","tumor","normal")
cat("samples:", table(grp)["tumor"], "tumour /", table(grp)["normal"], "normal | genes", nrow(x), "\n")

deg <- function(xx, g) {
  f <- factor(g, levels=c("normal","tumor")); d <- model.matrix(~0+f); colnames(d) <- c("N","T")
  topTable(eBayes(contrasts.fit(lmFit(xx,d), makeContrasts(T-N, levels=d))), number=Inf)
}
tt <- deg(x, grp)
cat(sprintf("criterion 1  FDR<0.05: %d | smallest adj.P = %s\n",
            sum(tt$adj.P.Val<0.05), signif(min(tt$adj.P.Val),3)))

obs  <- sum(tt$P.Value < 0.001)
perm <- replicate(1000, sum(deg(x, sample(grp))$P.Value < 0.001))
empP <- mean(perm >= obs)
cat(sprintf("criterion 2  observed rawP<0.001 = %d | permuted median %.0f, 95%% percentile %.0f | empirical P = %.4f\n",
            obs, median(perm), quantile(perm,.95), empP))
cat(sprintf("verdict   criterion 1 %s   criterion 2 %s\n",
            ifelse(sum(tt$adj.P.Val<0.05)>0,"pass","fail"),
            ifelse(empP<0.05,"pass","fail")))
writeLines(capture.output(sessionInfo()), "logs/03b_perm_GSE17025_sessionInfo.txt")
cat("===== DONE =====\n")
