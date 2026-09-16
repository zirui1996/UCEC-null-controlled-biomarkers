#!/usr/bin/env Rscript
# 19_prep_CIBERSORTx.R -- prepare the CIBERSORTx mixture file per section 18.1
# the input must be non-log TPM (not counts, not VST, not log2(TPM+1))
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages(library(SummarizedExperiment))
dir.create("results/immune/cibersortx_input", showWarnings=FALSE, recursive=TRUE)

se <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop=FALSE])
code <- substr(cdl$barcode,14,15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- code=="01" & cdl$primary_diagnosis %in% c(EN, SR)

tpm <- assay(se,"tpm_unstrand")[, keep]
rownames(tpm) <- rowData(se)$gene_name
colnames(tpm) <- cdl$barcode[keep]
# for duplicate symbols keep the one with the highest mean expression; drop all-zero rows
o <- order(rownames(tpm), -rowMeans(tpm)); i <- o[!duplicated(rownames(tpm)[o])]
tpm <- tpm[i, ]; tpm <- tpm[!is.na(rownames(tpm)) & rownames(tpm)!="" & rowSums(tpm)>0, ]
cat(sprintf("mixture: %d genes x %d samples\n", nrow(tpm), ncol(tpm)))
cat(sprintf("value range check (input must be non-log TPM): min=%.3f max=%.1f median column sum=%.0f\n",
    min(tpm), max(tpm), median(colSums(tpm))))
stopifnot(max(tpm) > 1000)   # log-transformed data could never reach values this large

out <- data.frame(GeneSymbol=rownames(tpm), tpm, check.names=FALSE)
f <- "results/immune/cibersortx_input/TCGA_UCEC_TPM_mixture.txt"
write.table(out, f, sep="\t", row.names=FALSE, quote=FALSE)
cat("wrote:", f, "\nsize:", round(file.size(f)/1024/1024,1), "MB\n")
