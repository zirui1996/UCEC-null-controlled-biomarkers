#!/usr/bin/env Rscript
# 23_geneX_pathways.R -- protocol section 20: pathways associated with Gene X
# split TCGA-UCEC into NUMA1-high/low at the median -> genome-wide limma moderated t -> Hallmark GSEA
# no direction is specified in advance (section 20)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(limma); library(edgeR); library(fgsea); library(ggplot2)})
dir.create("results/pathways/figures", showWarnings=FALSE, recursive=TRUE)
set.seed(42); GX <- "NUMA1"

se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop=FALSE])
code <- substr(cdl$barcode,14,15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- code=="01" & cdl$primary_diagnosis %in% c(EN, SR)
cnt <- assay(se,"unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]
cnt <- cnt[i,]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt)!="", ]

d <- DGEList(cnt); d <- d[filterByExpr(d),,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
lc <- cpm(d, log=TRUE)
gx <- lc[GX, ]
grp <- factor(ifelse(gx > median(gx), "High","Low"), levels=c("Low","High"))
cat(sprintf("TCGA-UCEC tumours: %d; NUMA1-high %d / low %d; genes tested %d\n\n",
    ncol(lc), sum(grp=="High"), sum(grp=="Low"), nrow(lc)))

## ---- limma: genome-wide moderated t ----
des <- model.matrix(~0+grp); colnames(des) <- c("Low","High")
v <- voom(d, des)
fit <- eBayes(contrasts.fit(lmFit(v, des), makeContrasts(High-Low, levels=des)))
tt <- topTable(fit, number=Inf); tt$gene <- rownames(tt)
cat(sprintf("NUMA1-high vs low: %d genes at FDR<0.05 (up %d / down %d)\n",
    sum(tt$adj.P.Val<0.05), sum(tt$adj.P.Val<0.05 & tt$logFC>0), sum(tt$adj.P.Val<0.05 & tt$logFC<0)))
write.table(tt[,c("gene","logFC","AveExpr","t","P.Value","adj.P.Val")],
            "results/pathways/NUMA1_high_vs_low_limma.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## ---- Hallmark GSEA (interface frozen to fgseaMultilevel, section 11.4) ----
hm <- readLines("data_raw/genesets/h.all.v2025.1.Hs.symbols.gmt")
pw <- setNames(lapply(strsplit(hm,"\t"), function(x) x[-(1:2)]), sapply(strsplit(hm,"\t"), `[`, 1))
ranks <- sort(setNames(tt$t, tt$gene), decreasing=TRUE)
fg <- fgseaMultilevel(pathways=pw, stats=ranks, minSize=15, maxSize=500)
fg <- fg[order(fg$padj), ]
cat(sprintf("\nHallmark GSEA: FDR<0.05 %d / %d（NES>0 %d / NES<0 %d）\n",
    sum(fg$padj<0.05), nrow(fg), sum(fg$padj<0.05 & fg$NES>0), sum(fg$padj<0.05 & fg$NES<0)))
o2 <- as.data.frame(fg); o2$leadingEdge <- sapply(fg$leadingEdge, paste, collapse=";")
write.table(o2, "results/pathways/NUMA1_Hallmark_GSEA.tsv", sep="\t", row.names=FALSE, quote=FALSE)

sig <- subset(fg, padj<0.05)
cat("\n--- top 12 sets enriched in NUMA1-high (NES>0) ---\n")
u <- head(sig[sig$NES>0, ][order(-sig[sig$NES>0,]$NES), ], 12)
print(data.frame(pathway=sub("^HALLMARK_","",u$pathway), NES=round(u$NES,2), padj=signif(u$padj,3)), row.names=FALSE)
cat("\n--- top 12 sets enriched in NUMA1-low (NES<0) ---\n")
dn <- head(sig[sig$NES<0, ][order(sig[sig$NES<0,]$NES), ], 12)
print(data.frame(pathway=sub("^HALLMARK_","",dn$pathway), NES=round(dn$NES,2), padj=signif(dn$padj,3)), row.names=FALSE)

## ---- required by section 20: summarise one tumour theme and one immune theme ----
IMMUNE <- c("INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","INFLAMMATORY_RESPONSE",
            "ALLOGRAFT_REJECTION","IL6_JAK_STAT3_SIGNALING","COMPLEMENT","IL2_STAT5_SIGNALING","TNFA_SIGNALING_VIA_NFKB")
im <- sig[sub("^HALLMARK_","",sig$pathway) %in% IMMUNE, ]
cat(sprintf("\n=== immune theme: %d Hallmark immune sets significant, of which %d have NES<0 ===\n", nrow(im), sum(im$NES<0)))
if (nrow(im)) print(data.frame(pathway=sub("^HALLMARK_","",im$pathway), NES=round(im$NES,2), padj=signif(im$padj,3)), row.names=FALSE)

## ---- plots ----
top <- rbind(head(sig[sig$NES>0,][order(-sig[sig$NES>0,]$NES),],15), head(sig[sig$NES<0,][order(sig[sig$NES<0,]$NES),],15))
top$lab <- sub("^HALLMARK_","",top$pathway)
p <- ggplot(top, aes(x=reorder(lab, NES), y=NES, fill=NES>0)) + geom_col() + coord_flip() +
  scale_fill_manual(values=c("TRUE"="#B2182B","FALSE"="#2166AC"), labels=c("NUMA1-low","NUMA1-high"), name="") +
  labs(x=NULL, y="Normalized enrichment score", title="Hallmark GSEA: NUMA1-high vs low (TCGA-UCEC)") +
  theme_classic(base_size=10)
ggsave("results/pathways/figures/NUMA1_Hallmark_GSEA.pdf", p, width=8, height=8)

## ---- co-expression network (descriptive) ----
cr <- apply(lc, 1, function(r) cor(gx, r, method="spearman"))
cr <- sort(cr, decreasing=TRUE)
cat("\n=== top 15 genes positively co-expressed with NUMA1 ===\n"); print(round(head(cr[names(cr)!=GX],15),3))
cat("\n=== top 15 negatively correlated ===\n"); print(round(tail(cr,15),3))
write.table(data.frame(gene=names(cr), spearman_rho=cr), "results/pathways/NUMA1_coexpression.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

writeLines(capture.output(sessionInfo()), "logs/23_pathways_sessionInfo.txt")
cat("\n===== 23_geneX_pathways DONE =====\n")
