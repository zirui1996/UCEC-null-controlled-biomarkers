#!/usr/bin/env Rscript
# 18_ssGSEA_checkpoints.R -- protocol section 19: immune checkpoints + ssGSEA
# section 18.7: all immune analyses are done on TCGA-UCEC bulk data; no deconvolution is run on GSE17025 (LCM epithelium)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(GSVA); library(org.Hs.eg.db)})
dir.create("results/immune/figures", showWarnings=FALSE, recursive=TRUE)
set.seed(42)
GX <- "NUMA1"

## ---- TPM for the TCGA-UCEC tumours ----
se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop=FALSE])
code <- substr(cdl$barcode,14,15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- code=="01" & cdl$primary_diagnosis %in% c(EN, SR)
tpm <- assay(se, "tpm_unstrand")[, keep]; rownames(tpm) <- rowData(se)$gene_name
tpm <- tpm[!duplicated(rownames(tpm)) & !is.na(rownames(tpm)) & rownames(tpm)!="", ]
tpm <- tpm[rowSums(tpm > 0) >= 0.2*ncol(tpm), ]
lt <- log2(tpm + 1)
cat(sprintf("TCGA-UCEC tumours: %d samples x %d genes (TPM, as defined in section 18.1)\n\n", ncol(lt), nrow(lt)))
gx <- lt[GX, ]

## ---- section 19.1 fixed checkpoint panel (locked before analysis, all reported) ----
CP <- c("PDCD1","CD274","CTLA4","LAG3","TIGIT","HAVCR2")
cat("=== section 19.1 immune checkpoints (fixed panel, Spearman + BH) ===\n")
r <- t(sapply(CP, function(g) {
  if (!g %in% rownames(lt)) return(c(rho=NA, P=NA))
  ct <- suppressWarnings(cor.test(gx, lt[g,], method="spearman"))
  c(rho=unname(ct$estimate), P=ct$p.value) }))
r <- as.data.frame(r); r$BH_q <- p.adjust(r$P, "BH")
r$direction <- ifelse(r$rho>0, "positive", "negative")
print(data.frame(gene=rownames(r), rho=round(r$rho,3), P=signif(r$P,3),
                 BH_q=signif(r$BH_q,3), direction=r$direction), row.names=FALSE)
cat(sprintf("\ncheckpoints at BH q<0.05: %d / %d (all six are reported, not only the positive ones, section 19.1)\n",
            sum(r$BH_q<0.05, na.rm=TRUE), nrow(r)))
write.table(data.frame(gene=rownames(r), r), "results/immune/checkpoint_correlations.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

## ---- §19.2 ssGSEA ----
hm <- readLines("data_raw/genesets/h.all.v2025.1.Hs.symbols.gmt")
hm <- setNames(lapply(strsplit(hm,"\t"), function(x) x[-(1:2)]), sapply(strsplit(hm,"\t"), `[`, 1))
# additional immune signatures (GO definition, from the local org.Hs.eg.db, reproducible)
GO <- c(antigen_presentation="GO:0019882", T_cell_activation="GO:0042110",
        cytotoxicity="GO:0001913", IFNg_response="GO:0034341")
imm <- lapply(GO, function(g) unique(suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys=g, keytype="GOALL", columns="SYMBOL"))$SYMBOL))
gs <- c(hm, imm)
gs <- gs[sapply(gs, function(x) sum(x %in% rownames(lt))) >= 10]
cat(sprintf("\n=== section 19.2 ssGSEA: %d gene sets (Hallmark %d + immune signatures %d) ===\n",
            length(gs), sum(grepl("^HALLMARK", names(gs))), sum(!grepl("^HALLMARK", names(gs)))))

par <- ssgseaParam(exprData = as.matrix(lt), geneSets = gs, minSize = 10, maxSize = 500, normalize = TRUE)
ss <- gsva(par, verbose = FALSE)
cat("ssGSEA matrix:", nrow(ss), "x", ncol(ss), "\n")
saveRDS(ss, "results/immune/ssGSEA_scores.rds")

## Gene X high vs low (split at the median, Wilcoxon + BH)
grp <- factor(ifelse(gx > median(gx), "High","Low"), levels=c("Low","High"))
res <- data.frame(set=rownames(ss),
  median_low  = apply(ss[,grp=="Low"],1,median), median_high = apply(ss[,grp=="High"],1,median),
  P = apply(ss,1,function(v) wilcox.test(v ~ grp)$p.value),
  rho = apply(ss,1,function(v) suppressWarnings(cor(gx, v, method="spearman"))))
res$BH_q <- p.adjust(res$P,"BH"); res <- res[order(res$BH_q),]
write.table(res, "results/immune/ssGSEA_NUMA1_high_vs_low.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("\ngene sets at BH q<0.05: %d / %d\n", sum(res$BH_q<0.05), nrow(res)))
cat("\n--- top 15 by BH q ---\n")
print(data.frame(set=sub("^HALLMARK_","",head(res$set,15)), rho=round(head(res$rho,15),3),
                 q=signif(head(res$BH_q,15),3)), row.names=FALSE)
cat("\n--- the four additional immune signatures (all reported, significant or not) ---\n")
print(data.frame(set=res$set[res$set %in% names(imm)], rho=round(res$rho[res$set %in% names(imm)],3),
                 P=signif(res$P[res$set %in% names(imm)],3), q=signif(res$BH_q[res$set %in% names(imm)],3)),
      row.names=FALSE)

saveRDS(list(lt=lt, gx=gx, grp=grp, ss=ss, cp=r), "results/immune/immune_objects.rds")
writeLines(capture.output(sessionInfo()), "logs/18_ssGSEA_sessionInfo.txt")
cat("\n===== 18_ssGSEA_checkpoints DONE =====\n")
