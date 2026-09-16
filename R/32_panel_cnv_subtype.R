#!/usr/bin/env Rscript
# 32_panel_cnv_subtype.R -- CNV dosage effect (section 22.2) and molecular subtype differences (section 16.3 supplementary figure) for the core 4 panel
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR)})
set.seed(42)
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")

se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop=FALSE])
code <- substr(cdl$barcode,14,15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
tk <- code=="01" & cdl$primary_diagnosis %in% c(EN,SR)
cnt <- assay(se,"unstranded")[,tk]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]
cnt <- cnt[i,]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt)!="",]
d <- DGEList(cnt); d <- d[rowSums(cpm(cnt)>0.5)>=50,,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
lc <- cpm(d, log=TRUE); z <- t(scale(t(lc)))
colnames(z) <- colnames(lc) <- substr(cdl$barcode[tk],1,15)
g <- intersect(CORE, rownames(z))
panel <- -colMeans(z[g,,drop=FALSE])
cat("panel genes:", paste(g, collapse=", "), sprintf("| n=%d tumours\n\n", length(panel)))

## ---------- A. section 22.2 CNV dosage effect (per gene and for the panel) ----------
cn <- read.delim(gzfile("data_raw/Xena/TCGA_UCEC_GISTIC2_thresholded.gz"), row.names=1, check.names=FALSE)
colnames(cn) <- substr(colnames(cn),1,15)
cm <- intersect(names(panel), colnames(cn))
cat(sprintf("=== section 22.2 CNV dosage effect (%d shared samples) ===\n", length(cm)))
lab <- c("-2","-1","0","+1","+2")
cnvres <- data.frame()
for (x in c(g,"panel")) {
  if (x=="panel") {
    cat("  panel       NOT COMPUTED: the four genes lie on chr4 / chr8 / chr8 / chr19,\n")
    cat("              averaging GISTIC copy number has no biological meaning; dosage must be tested per gene\n")
    next }
  if (!x %in% rownames(cn)) { cat(sprintf("  %-11s not in the GISTIC matrix\n", x)); next }
  cv <- as.numeric(cn[x, cm]); ev <- as.numeric(lc[x, cm])
  f <- factor(cv, levels=-2:2, labels=lab)
  kw <- kruskal.test(ev ~ f); s <- suppressWarnings(cor.test(cv, ev, method="spearman"))
  md <- tapply(ev, f, median); nn <- table(f)
  cat(sprintf("  %-11s KW P = %-9.3g Spearman rho = %+.3f (P=%.3g) | n: %s | median: %s\n", x, kw$p.value,
      s$estimate, s$p.value, paste(sprintf("%s=%d",lab,as.integer(nn)),collapse=" "),
      paste(sprintf("%s=%.2f",lab,md),collapse=" ")))
  cnvres <- rbind(cnvres, data.frame(gene=x, test="Kruskal", stat=kw$statistic, P=kw$p.value, n=length(cm)))
}
cnvres$BH_q <- p.adjust(cnvres$P,"BH")
cat("\nafter BH adjustment:\n"); print(data.frame(gene=cnvres$gene, P=signif(cnvres$P,3), BH_q=signif(cnvres$BH_q,3)), row.names=FALSE)
write.table(cnvres, "results/cnv/panel_cnv.tsv", sep="\t", row.names=FALSE, quote=FALSE)

pdf("results/cnv/figures/core4_CNV_dosage.pdf", width=9, height=8)
par(mfrow=c(2,2), mar=c(4,4,3,1))
for (x in g) { if (!x %in% rownames(cn)) next
  cv <- as.numeric(cn[x,cm]); ev <- as.numeric(lc[x,cm]); f <- factor(cv, levels=-2:2, labels=lab)
  boxplot(ev ~ f, col="grey85", outline=FALSE, xlab="GISTIC2 call", ylab="log2 CPM",
          main=sprintf("%s  KW P = %.2g", x, kruskal.test(ev~f)$p.value))
  points(jitter(as.numeric(f),.5), ev, pch=16, cex=.3, col=adjustcolor("black",.3)) }
dev.off()

## ---------- B. molecular subtypes ----------
cd2 <- cdl[tk,]; cd2$pid15 <- substr(cd2$barcode,1,15)
st <- as.character(cd2$paper_IntegrativeCluster)[match(names(panel), cd2$pid15)]   # as.character is required: ifelse returns integer codes for a factor
st <- factor(ifelse(st %in% c("CN high","CN low","MSI","POLE"), st, NA), levels=c("CN low","CN high","MSI","POLE"))
cat(sprintf("\n=== molecular subtypes (available for %d / %d) ===\n", sum(!is.na(st)), length(panel)))
print(table(st))
kw <- kruskal.test(panel[!is.na(st)] ~ droplevels(st[!is.na(st)]))
cat(sprintf("\nKruskal-Wallis P for the panel score across subtypes = %.3g\n", kw$p.value))
print(round(tapply(panel[!is.na(st)], droplevels(st[!is.na(st)]), median),3))
pw <- pairwise.wilcox.test(panel[!is.na(st)], droplevels(st[!is.na(st)]), p.adjust.method="BH")
cat("\npairwise comparisons (BH-adjusted):\n"); print(signif(pw$p.value,3))

cat("\nper gene across subtypes (BH adjustment over 4 tests):\n")
pg <- sapply(g, function(x) kruskal.test(as.numeric(z[x,!is.na(st)]) ~ droplevels(st[!is.na(st)]))$p.value)
print(data.frame(gene=names(pg), KW_P=signif(pg,3), BH_q=signif(p.adjust(pg,"BH"),3)), row.names=FALSE)

pdf("results/survival/figures/core4_panel_by_subtype.pdf", width=5.5, height=5)
boxplot(panel[!is.na(st)] ~ droplevels(st[!is.na(st)]), col="grey85", outline=FALSE,
        ylab="core-4 panel score", xlab="TCGA molecular subtype",
        main=sprintf("Kruskal-Wallis P = %.3g", kw$p.value))
points(jitter(as.numeric(droplevels(st[!is.na(st)])),.6), panel[!is.na(st)], pch=16, cex=.5, col=adjustcolor("black",.4))
dev.off()
cat("\n===== 32_panel_cnv_subtype DONE =====\n")
