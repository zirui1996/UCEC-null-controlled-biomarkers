#!/usr/bin/env Rscript
# 03_QC.R -- Phase 1, Step 2
# per-cohort quality control. Protocol section 8
# usage: Rscript R/03_QC.R GSE17025
# note: QC plots are not used to drop samples on the basis of how favourable the result looks (protocol section 8)

args <- commandArgs(trailingOnly = TRUE)
gse  <- args[1]
proj <- Sys.getenv("UCEC_PROJ", unset = getwd()); setwd(proj)
suppressPackageStartupMessages({ library(ggplot2) })

man <- read.delim("metadata/sample_manifest.tsv", stringsAsFactors = FALSE, quote = "")
man <- man[man$dataset == gse & man$include_primary == "yes", ]

x <- read.delim(sprintf("data_processed/expression/%s_gene_matrix_bestprobe.tsv", gse),
                check.names = FALSE, row.names = 1)
x <- as.matrix(x)
man <- man[match(colnames(x), man$sample_id), ]
grp <- factor(man$normal_or_tumor, levels = c("normal", "tumor"))

outdir <- "results/QC"; dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(outdir, sprintf("%s_QC.pdf", gse)), width = 9, height = 7)

# 1) density
d <- do.call(rbind, lapply(seq_len(ncol(x)), function(i)
  data.frame(value = x[, i], sample = colnames(x)[i], group = grp[i])))
print(ggplot(d, aes(value, group = sample, colour = group)) + geom_density(alpha = .3) +
      labs(title = paste(gse, "- expression density"), x = "log2 RMA") + theme_classic())

# 2) boxplot
par(mar = c(8, 4, 3, 1))
boxplot(x, las = 2, outline = FALSE, cex.axis = .5,
        col = c("grey70", "indianred")[as.integer(grp)],
        main = paste(gse, "- per-sample boxplot"), ylab = "log2 RMA")
legend("topright", legend = levels(grp), fill = c("grey70", "indianred"), bty = "n")

# 3) PCA
pca <- prcomp(t(x), scale. = FALSE)
ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pd  <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group = grp, id = colnames(x),
                  title = man$title)
p <- ggplot(pd, aes(PC1, PC2, shape = group, colour = group)) +
  geom_point(size = 3) +
  labs(title = paste(gse, "- PCA"), x = sprintf("PC1 (%.1f%%)", ve[1]),
       y = sprintf("PC2 (%.1f%%)", ve[2])) + theme_classic()

# GSE115810: flag the 4 samples whose labels conflict with GSE36389 (protocol sections 9.4 / 8)
if (gse == "GSE115810") {
  disputed <- c("endometrium13", "endometrium16", "endometrium17", "endometrium24")
  pd$disputed <- sub("_.*$", "", pd$title) %in% disputed
  p <- p + geom_point(data = subset(pd, disputed), aes(PC1, PC2),
                      shape = 21, size = 6, colour = "blue", fill = NA, stroke = 1.2) +
    ggrepel::geom_text_repel(data = subset(pd, disputed), aes(label = title), size = 3) +
    labs(subtitle = "blue circles = the 4 samples labelled control in GSE36389 (protocol §9.4)")
}
print(p)

# 4) sample correlation heatmap
cm <- cor(x, method = "spearman")
heatmap(cm, symm = TRUE, margins = c(8, 8), cexRow = .5, cexCol = .5,
        main = paste(gse, "- Spearman correlation"),
        ColSideColors = c("grey70", "indianred")[as.integer(grp)])

# 5) hierarchical clustering
hc <- hclust(as.dist(1 - cm), method = "average")
plot(hc, labels = paste0(substr(man$normal_or_tumor, 1, 1), "_", man$sample_id),
     cex = .5, main = paste(gse, "- hierarchical clustering (1 - Spearman)"), xlab = "")
dev.off()

# QC summary
qs <- data.frame(sample_id = colnames(x), group = as.character(grp),
                 median = apply(x, 2, median), iqr = apply(x, 2, IQR),
                 cor_to_median = round(apply(cm, 2, function(z) median(z[z < 1])), 3),
                 PC1 = round(pca$x[,1], 2), PC2 = round(pca$x[,2], 2))
write.table(qs, file.path(outdir, sprintf("%s_QC_summary.tsv", gse)),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("[", gse, "] QC done ->", file.path(outdir, sprintf("%s_QC.pdf", gse)), "\n")
cat("  samples:", ncol(x), " genes:", nrow(x), "\n")
cat("  variance explained by PC1/PC2:", ve[1], "/", ve[2], "%\n")
if (gse == "GSE115810") {
  cat("\n  ---- position of the 4 disputed samples in the PCA (recorded as required by protocol §9.4) ----\n")
  print(pd[pd$disputed, c("title", "group", "PC1", "PC2")], row.names = FALSE)
  cat("  ---- 3 prespecified controls ----\n")
  print(pd[pd$group == "normal", c("title", "PC1", "PC2")], row.names = FALSE)
}
