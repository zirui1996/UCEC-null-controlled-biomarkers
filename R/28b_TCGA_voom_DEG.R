#!/usr/bin/env Rscript
# 28b_TCGA_voom_DEG.R -- write the genome-wide voom-limma results for TCGA-UCEC tumour vs adjacent normal
# identical to the ad hoc computation at line 67 of 28_signature_AUC.R (~grp, filterByExpr, TMM, voom); serves as source data for Figure 1c
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(limma)})
se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop=FALSE])
code <- substr(cdl$barcode,14,15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code=="01" & cdl$primary_diagnosis %in% c(EN,SR)) | code=="11"
cnt <- assay(se,"unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]
cnt <- cnt[i,]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt)!="",]
grp <- factor(ifelse(code[keep]=="11","Normal","Tumor"), levels=c("Normal","Tumor"))
d <- DGEList(cnt); d <- d[filterByExpr(d, group=grp),,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
des <- model.matrix(~grp); v <- voom(d, des); f <- eBayes(lmFit(v, des))
tt <- topTable(f, coef=2, number=Inf); tt$gene <- rownames(tt)
write.table(tt[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val")], "results/DEG/TCGA_UCEC_voom_all_genes.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("TCGA-UCEC: Tumor %d / Normal %d | genes %d | FDR<0.05 %d (%.1f%%) | +|log2FC|>=1 %d\n",
    sum(grp=="Tumor"), sum(grp=="Normal"), nrow(tt), sum(tt$adj.P.Val<0.05), 100*mean(tt$adj.P.Val<0.05),
    sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>=1)))
