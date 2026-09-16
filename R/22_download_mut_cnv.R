#!/usr/bin/env Rscript
# 22_download_mut_cnv.R -- protocol section 22: TCGA-UCEC mutations (MAF) and CNV (GISTIC thresholded)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages(library(TCGAbiolinks))
dir.create("data_raw/TCGA/mut_cnv", showWarnings=FALSE, recursive=TRUE)

## ---- MAF (simple somatic mutations) ----
qm <- GDCquery(project="TCGA-UCEC", data.category="Simple Nucleotide Variation",
               data.type="Masked Somatic Mutation",
               workflow.type="Aliquot Ensemble Somatic Variant Merging and Masking")
cat("MAF files:", nrow(getResults(qm)), "\n")
GDCdownload(qm, directory="data_raw/TCGA/mut_cnv/GDCdata", files.per.chunk=50)
maf <- GDCprepare(qm, directory="data_raw/TCGA/mut_cnv/GDCdata")
cat("MAF rows:", nrow(maf), "| patients:", length(unique(substr(maf$Tumor_Sample_Barcode,1,12))), "\n")
saveRDS(maf, "data_raw/TCGA/mut_cnv/TCGA_UCEC_maf.rds")

## ---- CNV (gene level, GISTIC thresholded) ----
qc <- GDCquery(project="TCGA-UCEC", data.category="Copy Number Variation",
               data.type="Gene Level Copy Number Scores")
r <- getResults(qc); cat("\nCNV files:", nrow(r), "\n")
if (nrow(r) > 0) {
  GDCdownload(qc, directory="data_raw/TCGA/mut_cnv/GDCdata")
  cnv <- GDCprepare(qc, directory="data_raw/TCGA/mut_cnv/GDCdata")
  cat("CNV dimensions:", paste(dim(cnv), collapse=" x "), "\n")
  saveRDS(cnv, "data_raw/TCGA/mut_cnv/TCGA_UCEC_cnv_gistic.rds")
} else cat("no results for this data.type; a different CNV product is required\n")
cat("\n===== 22_download DONE =====\n")
