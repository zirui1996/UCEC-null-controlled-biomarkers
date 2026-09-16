#!/usr/bin/env Rscript
# 02c_GSE56087.R -- Discovery 4: read FPKM -> limma -> admission QC (permutation test)
# 9 tumours vs 9 adjacent normals, unpaired (GEO provides no recoverable pairing information)
proj <- Sys.getenv("UCEC_PROJ", unset = getwd()); setwd(proj)
suppressPackageStartupMessages(library(limma)); set.seed(42)

d <- "data_raw/GEO/GSE56087/fpkm"
fs <- sort(list.files(d, pattern="fpkm.txt.gz$", full.names=TRUE))
gsm <- sub("_.*$","", basename(fs))

grp_map <- c(GSM1354841="normal",GSM1354842="tumor",GSM1354843="normal",GSM1354844="tumor",
             GSM1354845="normal",GSM1354846="tumor",GSM1354847="normal",GSM1354848="tumor",
             GSM1354849="normal",GSM1354850="tumor",GSM1354851="tumor",GSM1354852="normal",
             GSM1354853="normal",GSM1354854="tumor",GSM1354855="normal",GSM1354856="tumor",
             GSM1354857="normal",GSM1354858="tumor")

lst <- lapply(fs, function(f) {
  t <- read.delim(gzfile(f), stringsAsFactors=FALSE)
  t <- t[t$FPKM_status=="OK" & !is.na(t$FPKM), ]
  # for duplicate symbols keep the highest FPKM (consistent with the microarray best-probe rule)
  ag <- aggregate(FPKM ~ gene_short_name, data=t, FUN=max)
  setNames(ag$FPKM, ag$gene_short_name)
})
genes <- Reduce(intersect, lapply(lst, names))
m <- sapply(lst, function(v) v[genes]); rownames(m) <- genes; colnames(m) <- gsm
cat("genes before filtering:", nrow(m), " samples:", ncol(m), "\n")

grp <- unname(grp_map[colnames(m)])
cat("groups:", paste(names(table(grp)), table(grp), collapse="  "), "\n")

# low-expression filter: FPKM>=1 in at least as many samples as the smallest group (9)
keep <- rowSums(m >= 1) >= 9
m <- m[keep, ]; cat("genes after the low-expression filter:", nrow(m), "\n")
x <- log2(m + 1)

write.table(data.frame(gene=rownames(x), x, check.names=FALSE),
            "data_processed/expression/GSE56087_gene_matrix_bestprobe.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

deg <- function(xx, g) {
  f <- factor(g, levels=c("normal","tumor")); d <- model.matrix(~0+f); colnames(d) <- c("N","T")
  topTable(eBayes(contrasts.fit(lmFit(xx,d), makeContrasts(T-N, levels=d)), trend=TRUE), number=Inf)
}
tt <- deg(x, grp); tt$gene <- rownames(tt)
cat(sprintf("\ncriterion 1  FDR<0.05: %d | +|lfc|>=1: %d | +>=2: %d | smallest adj.P = %s\n",
            sum(tt$adj.P.Val<0.05), sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>=1),
            sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>=2), signif(min(tt$adj.P.Val),3)))

obs <- sum(tt$P.Value<0.001)
perm <- replicate(1000, sum(deg(x, sample(grp))$P.Value < 0.001))
empP <- mean(perm >= obs)
cat(sprintf("criterion 2  observed rawP<0.001 = %d | permuted median %.0f, 95%% percentile %.0f | empirical P = %.4f\n",
            obs, median(perm), quantile(perm,.95), empP))
cat(sprintf("verdict   criterion 1 %s   criterion 2 %s\n",
            ifelse(sum(tt$adj.P.Val<0.05)>0,"pass","fail"),
            ifelse(empP<0.05,"pass","fail")))

o <- tt[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val")]
write.table(o, "results/DEG/GSE56087_all_genes.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(subset(o, adj.P.Val<0.05 & abs(logFC)>=1),
            "results/DEG/GSE56087_DEG_primary.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("\n===== GSE56087 DONE =====\n")
