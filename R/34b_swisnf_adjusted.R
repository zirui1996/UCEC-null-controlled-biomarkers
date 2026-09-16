#!/usr/bin/env Rscript
# 34b_swisnf_adjusted.R -- targeted SWI/SNF tests: BH adjustment, TMB adjustment, and exclusion of hypermutated subtypes
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

## SWI/SNF: raw Wilcoxon + BH(8) | logistic mutation ~ panel(z) + log10(TMB) | CN low / CN high only
SW <- c("ARID1A","ARID1B","ARID2","SMARCA4","SMARCB1","SMARCC1","SMARCC2","PBRM1")
dd$pz <- as.numeric(scale(dd$panel)); dd$ltmb <- log10(dd$tmb+1)
res <- do.call(rbind, lapply(SW, function(x) {
  dd$has <- as.integer(dd$pid %in% mut$pid[mut$Hugo_Symbol==x])
  p_raw <- wilcox.test(dd$panel[dd$has==1], dd$panel[dd$has==0])$p.value
  f <- summary(glm(has ~ pz + ltmb, data=dd, family=binomial))$coefficients
  nm <- dd[which(dd$subtype %in% c("CN low","CN high")),]
  p_nm <- if (sum(nm$has)>=5 && sum(nm$has)<nrow(nm)-5) wilcox.test(nm$panel[nm$has==1], nm$panel[nm$has==0])$p.value else NA
  data.frame(gene=x, n_mut=sum(dd$has), P_raw=p_raw,
             OR_panel_adjTMB=exp(f["pz","Estimate"]), P_adjTMB=f["pz","Pr(>|z|)"],
             n_mut_nonhyper=sum(nm$has), P_nonhyper=p_nm)
}))
res$q_raw <- p.adjust(res$P_raw,"BH"); res$q_adjTMB <- p.adjust(res$P_adjTMB,"BH")
res$q_nonhyper <- p.adjust(res$P_nonhyper,"BH")
cat(sprintf("median TMB by subtype: %s\n", paste(tapply(dd$tmb, dd$subtype, median), names(tapply(dd$tmb, dd$subtype, median)), sep=" ", collapse=" | ")))
print(transform(res, P_raw=signif(P_raw,3), q_raw=signif(q_raw,3), OR_panel_adjTMB=round(OR_panel_adjTMB,2),
      P_adjTMB=signif(P_adjTMB,3), q_adjTMB=signif(q_adjTMB,3), P_nonhyper=signif(P_nonhyper,3), q_nonhyper=signif(q_nonhyper,3)), row.names=FALSE)
write.table(res, "results/mutation/panel_SWISNF_adjusted.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("\nBH q<0.05  unadjusted: %d/8 | TMB-adjusted: %d/8 | non-hypermutated subtypes: %d/%d\n",
    sum(res$q_raw<0.05), sum(res$q_adjTMB<0.05), sum(res$q_nonhyper<0.05, na.rm=TRUE), sum(!is.na(res$q_nonhyper))))
cat("===== 34b DONE =====\n")
