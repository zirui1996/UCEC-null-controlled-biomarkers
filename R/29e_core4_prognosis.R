#!/usr/bin/env Rscript
# 29e_core4_prognosis.R -- overall-survival analysis of the core 4 genes and the panel (source data for Table 2 / Figure S3)
# v2.6: results/survival/core4_prognosis.tsv previously had no script behind it; this script recomputes and checks it on the basis of section 4 of R/29_core4_panel.R
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(survival)})
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")
se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop=FALSE])
code <- substr(cdl$barcode,14,15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code=="01" & cdl$primary_diagnosis %in% c(EN,SR)) | code=="11"
cnt <- assay(se,"unstranded")[,keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]
cnt <- cnt[i,]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt)!="",]
grp <- factor(ifelse(code[keep]=="11","Normal","Tumor"), levels=c("Normal","Tumor"))
d <- DGEList(cnt); d <- d[filterByExpr(d,group=grp),,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
z <- t(scale(t(cpm(d, log=TRUE))))
sv <- readRDS("results/survival/survival_objects.rds")$d
# the row names of survival_objects are the original row numbers of the tumour samples in R/16, not patient IDs; matching on row names as R/29 did fails for every sample.
# barcodes are reconstructed here from the sample order used in R/16, and each match is verified against the stored NUMA1_z values.
go <- readRDS("results/validation/TCGA_geneX_objects.rds")
tb <- go$cd$barcode[go$grp=="Tumor"]
sv$barcode <- tb[as.integer(rownames(sv))]
chkz <- as.numeric(scale(go$lc["NUMA1", go$grp=="Tumor"]))[as.integer(rownames(sv))]
cat(sprintf("row-index reconstruction check: largest NUMA1_z difference %.2e (n=%d)\n", max(abs(chkz - sv$NUMA1_z)), nrow(sv)))
stopifnot(max(abs(chkz - sv$NUMA1_z)) < 1e-8)
tz <- z[, grp=="Tumor"]; colnames(tz) <- cdl$barcode[keep][grp=="Tumor"]
m <- match(sv$barcode, colnames(tz))
cat(sprintf("barcodes matched %d / %d\n", sum(!is.na(m)), nrow(sv)))
sv$pid <- substr(sv$barcode, 1, 12)
for (g in CORE) sv[[g]] <- as.numeric(tz[g, m])
sv$panel <- -colMeans(tz[CORE, m, drop=FALSE])
sv <- sv[!is.na(sv$panel), ]
cat(sprintf("n = %d | events = %d | median follow-up %d days (reverse KM: %d)\n", nrow(sv), sum(sv$OS_event),
    as.integer(median(sv$OS_time)), as.integer(summary(survfit(Surv(OS_time, 1-OS_event) ~ 1, data=sv))$table["median"])))
res <- do.call(rbind, lapply(c(CORE, "panel"), function(x) {
  sv$tmp <- as.numeric(scale(sv[[x]])); s <- summary(coxph(Surv(OS_time,OS_event) ~ tmp, data=sv))
  sv$hi <- factor(ifelse(sv[[x]] > median(sv[[x]]), "High", "Low"), levels=c("Low","High"))
  lr <- survdiff(Surv(OS_time,OS_event) ~ hi, data=sv)
  data.frame(gene=x, HR=s$coef[,"exp(coef)"], lo=s$conf.int[,3], hi=s$conf.int[,4], P=s$coef[,"Pr(>|z|)"],
             C=unname(s$concordance[1]), KM_P=1-pchisq(lr$chisq, 1), n=nrow(sv), events=sum(sv$OS_event))
}))
res$BH_q <- p.adjust(res$P, "BH")
print(res, row.names=FALSE, digits=3)
old <- read.delim("archive/unscripted_results/core4_prognosis_n500_adhoc.tsv")
cmp <- merge(res, old, by="gene", suffixes=c("", "_old"))
for (k in c("HR","lo","hi","P","C","KM_P","BH_q")) cat(sprintf("  %-5s largest difference from the original table %.2e\n", k, max(abs(cmp[[k]] - cmp[[paste0(k,"_old")]]))))
write.table(res, "results/survival/core4_prognosis.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(sv[, c("pid","OS_time","OS_event","age","stage","grade","histo","subtype",CORE,"panel")],
            "results/supp/survival_core4_per_patient.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("===== 29e DONE =====\n")
