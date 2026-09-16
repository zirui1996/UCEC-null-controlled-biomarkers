#!/usr/bin/env Rscript
# 14_download_TCGA.R -- protocol section 7.4: download TCGA-UCEC harmonized data (GDC)
# expression: STAR counts are the primary data for statistics; the same object also carries TPM/FPKM for display and immune deconvolution
# normal samples are identified by characters 14-15 of the barcode being "11" (Solid Tissue Normal)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(TCGAbiolinks); library(SummarizedExperiment)})

q <- GDCquery(project = "TCGA-UCEC",
              data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification",
              workflow.type = "STAR - Counts")
res <- getResults(q)
cat("files returned by the query:", nrow(res), "\n")
cat("sample_type distribution:\n"); print(table(res$sample_type))

GDCdownload(q, directory = "data_raw/TCGA/GDCdata", files.per.chunk = 40)
se <- GDCprepare(q, directory = "data_raw/TCGA/GDCdata", summarizedExperiment = TRUE)
cat("\nSE dimensions:", dim(se)[1], "genes x", dim(se)[2], "samples\n")
cat("assays:", paste(assayNames(se), collapse=", "), "\n")

cd <- as.data.frame(colData(se))
cd$participant <- substr(cd$barcode, 1, 12)
cd$sample_code <- substr(cd$barcode, 14, 15)
cd$group <- ifelse(cd$sample_code == "11", "Normal",
            ifelse(cd$sample_code == "01", "Tumor", NA))
cat("\nclassified by characters 14-15 of the barcode:\n"); print(table(cd$sample_code, cd$group, useNA="ifany"))
np <- length(intersect(cd$participant[cd$group=="Tumor"], cd$participant[cd$group=="Normal"]))
cat("patients with both a tumour and a normal sample:", np, "\n")

saveRDS(se, "data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
write.table(cd, "metadata/TCGA_UCEC_colData.tsv", sep="\t", row.names=FALSE, quote=FALSE)
writeLines(c(paste("download date:", Sys.Date()),
             paste("TCGAbiolinks:", as.character(packageVersion("TCGAbiolinks"))),
             paste("workflow: STAR - Counts"),
             paste("files:", nrow(res)), paste("SE:", paste(dim(se), collapse=" x ")),
             paste("matched pairs:", np)), "logs/TCGA_download_provenance.txt")
cat("\n===== 14_download_TCGA DONE =====\n")
