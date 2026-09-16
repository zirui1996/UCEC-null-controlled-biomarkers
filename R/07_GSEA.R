#!/usr/bin/env Rscript
# 07_GSEA.R -- Phase 3 / protocol section 11.4: pre-ranked Hallmark GSEA
# the interface is frozen: fgseaMultilevel (not fgseaSimple(nperm=)), minSize 15 / maxSize 500, FDR<0.05
# each voting cohort is ranked by its own limma t statistic; cohorts are not pooled (section 7.3 forbids concatenating across platforms)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(fgsea); library(ggplot2)})
dir.create("results/enrichment/figures", showWarnings=FALSE, recursive=TRUE)
set.seed(42)

GMT <- "data_raw/genesets/h.all.v2025.1.Hs.symbols.gmt"
pathways <- gmtPathways(GMT)
cat("Hallmark gene sets:", length(pathways), " source:", GMT, "\n")
cat("MD5:", system(sprintf("md5 -q %s", GMT), intern=TRUE), "\n")

VOTE <- c("GSE17025","GSE63678","GSE106191")
res <- list()
for (g in VOTE) {
  t <- read.delim(sprintf("results/DEG/%s_all_genes.tsv", g), stringsAsFactors=FALSE)
  t <- t[!duplicated(t$gene) & !is.na(t$t), ]
  ranks <- setNames(t$t, t$gene); ranks <- sort(ranks, decreasing=TRUE)
  fg <- fgseaMultilevel(pathways=pathways, stats=ranks, minSize=15, maxSize=500)
  fg <- fg[order(fg$padj), ]
  res[[g]] <- fg
  cat(sprintf("%-11s testable pathways %d | FDR<0.05 %d (NES>0 %d / NES<0 %d)\n", g, nrow(fg),
      sum(fg$padj<0.05), sum(fg$padj<0.05 & fg$NES>0), sum(fg$padj<0.05 & fg$NES<0)))
  o <- as.data.frame(fg); o$leadingEdge <- sapply(fg$leadingEdge, paste, collapse=";")
  write.table(o, sprintf("results/enrichment/GSEA_Hallmark_%s.tsv", g), sep="\t", row.names=FALSE, quote=FALSE)
}

## cross-cohort summary: in how many cohorts the set reaches FDR<0.05 with the same NES sign
allp <- unique(unlist(lapply(res, function(x) x$pathway)))
mat <- sapply(VOTE, function(g){ x <- res[[g]]; setNames(x$NES, x$pathway)[allp] })
sig <- sapply(VOTE, function(g){ x <- res[[g]]; setNames(x$padj, x$pathway)[allp] < 0.05 })
rownames(mat) <- rownames(sig) <- allp
sig[is.na(sig)] <- FALSE
nsig <- rowSums(sig)
dir_ok <- apply(cbind(sig, mat), 1, function(r){
  s <- which(as.logical(r[1:3])); if(!length(s)) return(NA)
  length(unique(sign(as.numeric(r[4:6])[s])))==1 })

cat("\n=== cross-cohort Hallmark summary (number of cohorts at FDR<0.05) ===\n"); print(table(n_cohorts=nsig))
sum3 <- data.frame(pathway=allp, n_sig=nsig, dir_consistent=dir_ok,
                   NES_GSE17025=mat[,1], NES_GSE63678=mat[,2], NES_GSE106191=mat[,3])
sum3 <- sum3[order(-sum3$n_sig, sum3$pathway), ]
write.table(sum3, "results/enrichment/GSEA_Hallmark_crosscohort.tsv", sep="\t", row.names=FALSE, quote=FALSE)

cat("\n--- significant in 3/3 cohorts with the same direction ---\n")
x <- subset(sum3, n_sig==3 & dir_consistent %in% TRUE)
print(data.frame(pathway=sub("^HALLMARK_","",x$pathway), round(x[,4:6],2)), row.names=FALSE)

## heatmap data used for Figure 3E
pdf("results/enrichment/figures/GSEA_Hallmark_heatmap.pdf", width=8, height=10)
m <- mat[nsig>=2, , drop=FALSE]; m <- m[order(-rowMeans(m, na.rm=TRUE)), , drop=FALSE]
par(mar=c(6,18,3,2))
image(t(m[nrow(m):1,]), col=colorRampPalette(c("#2166AC","white","#B2182B"))(51),
      zlim=c(-max(abs(m),na.rm=TRUE), max(abs(m),na.rm=TRUE)), axes=FALSE,
      main="Hallmark GSEA NES (FDR<0.05 in >=2 cohorts)")
axis(1, at=seq(0,1,length.out=3), labels=VOTE, las=2, cex.axis=0.8)
axis(2, at=seq(0,1,length.out=nrow(m)), labels=sub("^HALLMARK_","",rownames(m)[nrow(m):1]), las=2, cex.axis=0.6)
dev.off()

saveRDS(res, "results/enrichment/GSEA_objects.rds")
writeLines(capture.output(sessionInfo()), "logs/07_GSEA_sessionInfo.txt")
cat("\n===== 07_GSEA DONE =====\n")
