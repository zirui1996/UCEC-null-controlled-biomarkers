#!/usr/bin/env Rscript
# 34_panel_mut_pathway.R -- mutations (section 22.1) and pathways (section 20) for the core 4 panel
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(limma); library(fgsea); library(dplyr)})
set.seed(42); CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")

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
g <- intersect(CORE, rownames(z)); panel <- -colMeans(z[g,,drop=FALSE])
pid <- substr(cdl$barcode[tk],1,12); names(panel) <- pid
cat(sprintf("panel: %d tumours, genes %s\n\n", length(panel), paste(g, collapse=", ")))

## ---------- A. section 22.1 mutations ----------
maf <- as.data.frame(readRDS("data_raw/TCGA/mut_cnv/TCGA_UCEC_maf.rds"))
maf$pid <- substr(maf$Tumor_Sample_Barcode,1,12)
mut <- unique(maf[,c("pid","Hugo_Symbol")])
cd <- read.delim("metadata/TCGA_UCEC_colData.tsv", stringsAsFactors=FALSE); cd$pid <- substr(cd$barcode,1,12)
sub <- as.character(cd$paper_IntegrativeCluster)[match(names(panel), cd$pid)]
sub <- ifelse(sub %in% c("CN high","CN low","MSI","POLE"), sub, NA)
tmb <- table(maf$pid); com <- intersect(names(panel), names(tmb))
dd <- data.frame(pid=com, panel=panel[com], tmb=as.numeric(tmb[com]), subtype=sub[match(com,names(panel))])
dd <- dd[!duplicated(dd$pid),]
ct <- suppressWarnings(cor.test(dd$panel, log10(dd$tmb+1), method="spearman"))
cat(sprintf("=== section 22.1 mutations ===\npanel vs log10(TMB): rho = %.3f, P = %.3g\n", ct$estimate, ct$p.value))
dd$grp <- ifelse(dd$panel > median(dd$panel), "High","Low")
GEN <- names(which(table(mut$Hugo_Symbol) >= 0.05*nrow(dd)))
fis <- function(sd, tag) { if (nrow(sd)<20) return(NULL)
  r <- do.call(rbind, lapply(GEN, function(x){
    has <- sd$pid %in% mut$pid[mut$Hugo_Symbol==x]
    tb <- table(factor(sd$grp,c("Low","High")), factor(has,c(FALSE,TRUE)))
    if (any(dim(tb)!=2) || sum(tb[,2])<3) return(NULL)
    f <- fisher.test(tb); data.frame(stratum=tag, gene=x, n=nrow(sd),
      pct_low=round(100*tb[1,2]/sum(tb[1,]),1), pct_high=round(100*tb[2,2]/sum(tb[2,]),1),
      OR=unname(f$estimate), P=f$p.value)}))
  if(is.null(r)) return(NULL); r$BH_q <- p.adjust(r$P,"BH"); r[order(r$BH_q),] }
ra <- fis(dd,"all")
cat(sprintf("whole cohort, %d genes (mutated in >=5%% of patients), BH q<0.05: %d\n", nrow(ra), sum(ra$BH_q<0.05)))
print(head(data.frame(gene=ra$gene, low=ra$pct_low, high=ra$pct_high, OR=round(ra$OR,2), q=signif(ra$BH_q,3)),6), row.names=FALSE)
cat("\nstratified by molecular subtype (required by section 22.1):\n")
for (st in c("CN low","CN high","MSI","POLE")) { r <- fis(dd[which(dd$subtype==st),], st)
  if (is.null(r)) { cat(sprintf("  %-9s too few samples\n", st)); next }
  cat(sprintf("  %-9s n=%d, %d genes tested, BH q<0.05: %d\n", st, r$n[1], nrow(r), sum(r$BH_q<0.05))) }
cat("\ndirected test: SWI/SNF\n")
for (x in c("ARID1A","ARID1B","ARID2","SMARCA4","SMARCB1","SMARCC1","SMARCC2","PBRM1")) {
  has <- dd$pid %in% mut$pid[mut$Hugo_Symbol==x]; if (sum(has)<5) next
  cat(sprintf("  %-9s mutated %3d (%.1f%%) | median panel score, mutated %.3f vs wild-type %.3f | P = %.3g\n",
      x, sum(has), 100*mean(has), median(dd$panel[has]), median(dd$panel[!has]),
      wilcox.test(dd$panel[has], dd$panel[!has])$p.value)) }
write.table(ra, "results/mutation/panel_fisher_all.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## ---------- B. section 20 pathways ----------
grp <- factor(ifelse(panel > median(panel), "High","Low"), levels=c("Low","High"))
des <- model.matrix(~0+grp); colnames(des) <- c("Low","High")
v <- voom(d, des); fit <- eBayes(contrasts.fit(lmFit(v,des), makeContrasts(High-Low, levels=des)))
tt <- topTable(fit, number=Inf); tt$gene <- rownames(tt)
cat(sprintf("\n=== section 20 pathways: panel-high vs low, %d genes at FDR<0.05 ===\n", sum(tt$adj.P.Val<0.05)))
hm <- readLines("data_raw/genesets/h.all.v2025.1.Hs.symbols.gmt")
pw <- setNames(lapply(strsplit(hm,"\t"), function(x) x[-(1:2)]), sapply(strsplit(hm,"\t"),`[`,1))
fg <- fgseaMultilevel(pw, sort(setNames(tt$t, tt$gene), decreasing=TRUE), minSize=15, maxSize=500)
fg <- fg[order(fg$padj),]; sig <- subset(fg, padj<0.05)
cat(sprintf("Hallmark GSEA: %d/50 significant (NES>0 %d / NES<0 %d)\n", nrow(sig), sum(sig$NES>0), sum(sig$NES<0)))
cat("\ntop 10 sets enriched in panel-high:\n")
u <- head(sig[sig$NES>0,][order(-sig[sig$NES>0,]$NES),],10)
print(data.frame(pathway=sub("^HALLMARK_","",u$pathway), NES=round(u$NES,2), padj=signif(u$padj,3)), row.names=FALSE)
cat("\ntop 10 sets enriched in panel-low:\n")
dn <- head(sig[sig$NES<0,][order(sig[sig$NES<0,]$NES),],10)
print(data.frame(pathway=sub("^HALLMARK_","",dn$pathway), NES=round(dn$NES,2), padj=signif(dn$padj,3)), row.names=FALSE)
IM <- c("INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE","INFLAMMATORY_RESPONSE","ALLOGRAFT_REJECTION",
        "IL6_JAK_STAT3_SIGNALING","COMPLEMENT","IL2_STAT5_SIGNALING","TNFA_SIGNALING_VIA_NFKB")
im <- sig[sub("^HALLMARK_","",sig$pathway) %in% IM,]
cat(sprintf("\nimmune theme: %d significant sets, of which %d have NES<0\n", nrow(im), sum(im$NES<0)))
if (nrow(im)) print(data.frame(pathway=sub("^HALLMARK_","",im$pathway), NES=round(im$NES,2), padj=signif(im$padj,3)), row.names=FALSE)
o2 <- as.data.frame(fg); o2$leadingEdge <- sapply(fg$leadingEdge, paste, collapse=";")
write.table(o2, "results/pathways/panel_Hallmark_GSEA.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(tt[,c("gene","logFC","t","P.Value","adj.P.Val")], "results/pathways/panel_high_vs_low_limma.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
cat("\n===== 34_panel_mut_pathway DONE =====\n")
