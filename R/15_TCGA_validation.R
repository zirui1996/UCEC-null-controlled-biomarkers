#!/usr/bin/env Rscript
# 15_TCGA_validation.R -- protocol section 15
# formal name: independent cross-platform expression and discriminatory validation
# (this is not "clinical diagnostic validation": the controls are adjacent tissue, not a healthy screening population)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(pROC)})
dir.create("results/validation/figures", showWarnings=FALSE, recursive=TRUE)
set.seed(42)

GENES <- c("NUMA1","UCK2","HAND2-AS1","OSR2","PEG3","RUNX1T1","NAALAD2","ISG15","SPP1","ST3GAL5")
GX <- "NUMA1"

se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cd <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop=FALSE])
cd$participant <- substr(cd$barcode,1,12); cd$code <- substr(cd$barcode,14,15)
cd$group <- ifelse(cd$code=="11","Normal", ifelse(cd$code=="01","Tumor", NA))

## ---- sample selection (config: tcga_histology_primary = endometrioid_and_serous) ----
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SE_ <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
cd$histo <- ifelse(cd$primary_diagnosis %in% EN, "endometrioid",
            ifelse(cd$primary_diagnosis %in% SE_, "serous", "other"))
# the histology restriction applies to tumours only; adjacent normal endometrium is not classified by tumour histology, so all 35 samples are kept
keep <- (cd$group=="Tumor" & cd$histo %in% c("endometrioid","serous")) | (cd$group=="Normal")
keep[is.na(keep)] <- FALSE
cat("=== inclusion summary ===\n")
cat(sprintf("553 tumours -> %d eligible by histology (endometrioid %d / serous %d); %d excluded (other/not reported)\n",
    sum(keep & cd$group=="Tumor"), sum(keep & cd$histo=="endometrioid" & cd$group=="Tumor"),
    sum(keep & cd$histo=="serous" & cd$group=="Tumor"),
    sum(cd$group=="Tumor", na.rm=TRUE)-sum(keep & cd$group=="Tumor", na.rm=TRUE)))
cat(sprintf("all %d adjacent normal samples retained (the histology restriction does not apply to adjacent tissue)\n", sum(keep & cd$group=="Normal")))
cat(sprintf("recurrent tumours (code 02) excluded: %d\n", sum(cd$code=="02")))

se2 <- se[, keep]; cd2 <- cd[keep, ]
grp <- factor(cd2$group, levels=c("Normal","Tumor"))

## ---- expression: section 7.4 uses STAR counts as the primary data (TMM-logCPM); TPM serves as a sensitivity analysis ----
cnt <- assay(se2, "unstranded"); rownames(cnt) <- rowData(se2)$gene_name
tpm <- assay(se2, "tpm_unstrand"); rownames(tpm) <- rowData(se2)$gene_name
d <- DGEList(cnt); d <- d[filterByExpr(d, group=grp), , keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
lc <- cpm(d, log=TRUE)
lt <- log2(tpm + 1)
cat(sprintf("genes after filtering %d; samples %d (Tumour %d / Normal %d)\n\n", nrow(lc), ncol(lc),
    sum(grp=="Tumor"), sum(grp=="Normal")))

pairs <- intersect(cd2$participant[grp=="Tumor"], cd2$participant[grp=="Normal"])
cat("matched patients (observed after the histology filter):", length(pairs), "\n\n")

res <- data.frame()
for (g in GENES) {
  if (!g %in% rownames(lc)) { cat(sprintf("%-11s NOT IN THE FILTERED MATRIX\n", g)); next }
  x <- lc[g, ]
  w  <- wilcox.test(x ~ grp)
  r  <- roc(grp, x, quiet=TRUE, direction="auto"); ci <- as.numeric(ci.auc(r))
  # paired sensitivity analysis (§15.2.1): paired signed-rank; no ROC is reported on the paired subset
  pt <- sapply(pairs, function(p) mean(x[cd2$participant==p & grp=="Tumor"]))
  pn <- sapply(pairs, function(p) mean(x[cd2$participant==p & grp=="Normal"]))
  wp <- wilcox.test(pt, pn, paired=TRUE)
  # TPM sensitivity analysis
  xt <- lt[g, ]; wt <- wilcox.test(xt ~ grp)
  res <- rbind(res, data.frame(gene=g,
    normal_median=round(median(x[grp=="Normal"]),3), normal_IQR=sprintf("%.2f–%.2f", quantile(x[grp=="Normal"],.25), quantile(x[grp=="Normal"],.75)),
    tumor_median=round(median(x[grp=="Tumor"]),3),  tumor_IQR=sprintf("%.2f–%.2f", quantile(x[grp=="Tumor"],.25), quantile(x[grp=="Tumor"],.75)),
    direction=ifelse(median(x[grp=="Tumor"])>median(x[grp=="Normal"]),"up","down"),
    P_wilcox=signif(w$p.value,3), AUC=round(as.numeric(r$auc),3),
    AUC_CI=sprintf("%.3f–%.3f", ci[1], ci[3]),
    P_paired=signif(wp$p.value,3), n_pairs=length(pairs),
    P_TPM=signif(wt$p.value,3)))
}
res$FDR_wilcox <- signif(p.adjust(res$P_wilcox, "BH"), 3)
print(res[, c("gene","direction","normal_median","tumor_median","P_wilcox","FDR_wilcox","AUC","AUC_CI","P_paired","P_TPM")], row.names=FALSE)
write.table(res, "results/validation/TCGA_geneX_validation.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## ---- direction concordance check against the discovery cohorts ----
rec <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors=FALSE)
cat("\n=== direction concordance with the discovery cohorts ===\n")
for (i in seq_len(nrow(res))) {
  g <- res$gene[i]; dd <- rec$direction[rec$gene==g]
  if (length(dd)) cat(sprintf("  %-11s discovery %-4s | TCGA %-4s -> %s\n", g, dd, res$direction[i],
      ifelse(dd==res$direction[i],"concordant","CONFLICT")))
}

## ---- plots for Gene X ----
x <- lc[GX, ]
pdf(sprintf("results/validation/figures/TCGA_%s_boxplot.pdf", GX), width=4.2, height=5)
boxplot(x ~ grp, col=c("#4393C3","#D6604D"), outline=FALSE, ylab="log2 CPM (TMM)", xlab="",
        main=sprintf("%s in TCGA-UCEC\nP = %.2e", GX, wilcox.test(x~grp)$p.value))
points(jitter(as.numeric(grp), .6), x, pch=16, cex=.4, col=adjustcolor("black",.35)); dev.off()
r <- roc(grp, x, quiet=TRUE)
pdf(sprintf("results/validation/figures/TCGA_%s_ROC.pdf", GX), width=5, height=5)
plot(r, col="#B2182B", lwd=2, main=sprintf("%s — TCGA-UCEC\nAUC = %.3f (95%%CI %.3f–%.3f)",
     GX, as.numeric(r$auc), as.numeric(ci.auc(r))[1], as.numeric(ci.auc(r))[3])); dev.off()

saveRDS(list(lc=lc[GENES[GENES %in% rownames(lc)],], grp=grp, cd=cd2, pairs=pairs),
        "results/validation/TCGA_geneX_objects.rds")
writeLines(capture.output(sessionInfo()), "logs/15_TCGA_validation_sessionInfo.txt")
cat("\n===== 15_TCGA_validation DONE =====\n")
