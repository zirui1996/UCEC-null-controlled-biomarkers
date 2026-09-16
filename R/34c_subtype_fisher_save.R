#!/usr/bin/env Rscript
# 34c_subtype_fisher_save.R -- save the subtype-stratified Fisher results for the panel (lines 1-45 of R/34 reused verbatim; R/34 only printed them)
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
out <- list()
for (st in c("CN low","CN high","MSI","POLE")) { r <- fis(dd[which(dd$subtype==st),], st)
  if (is.null(r)) { cat(sprintf("  %-9s too few samples\n", st)); next }
  cat(sprintf("  %-9s n=%d, %d genes tested, BH q<0.05: %d\n", st, r$n[1], nrow(r), sum(r$BH_q<0.05))); out[[st]] <- r }
write.table(do.call(rbind, out), "results/mutation/panel_fisher_by_subtype.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("===== 34c DONE =====\n")
