#!/usr/bin/env Rscript
# 12_export_for_LOCO.R -- export per-cohort all-gene z-score matrices and labels for the nested LOCO analysis of section 14
# section 14.3: gene-wise z-scoring is done independently within each dataset (transductive; the wording is hedged accordingly in the protocol)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
dir.create("results/LOCO", showWarnings=FALSE, recursive=TRUE)
MAT <- c(GSE17025 ="data_processed/expression/GSE17025_gene_matrix_bestprobe.tsv",
         GSE63678 ="data_processed/expression/GSE63678_gene_matrix_bestprobe.tsv",
         GSE106191="data_processed/expression/GSE106191_all97_gene_matrix_bestprobe.tsv")
man <- read.delim("metadata/sample_manifest.tsv", stringsAsFactors=FALSE)
for (g in names(MAT)) {
  x <- as.matrix(read.delim(MAT[g], row.names=1, check.names=FALSE))
  z <- t(scale(t(x)))                                  # z-score within the dataset
  z <- z[complete.cases(z), , drop=FALSE]              # drop zero-variance genes
  if (g=="GSE106191") {
    m <- read.delim("metadata/sample_manifest_GSE106191.tsv", stringsAsFactors=FALSE)
  } else {
    m <- subset(man, dataset==g & include_primary=="yes")
  }
  grp <- setNames(ifelse(m$normal_or_tumor=="tumor","Tumor","Normal"), m$sample_id)[colnames(z)]
  stopifnot(!any(is.na(grp)))
  write.table(data.frame(gene=rownames(z), z, check.names=FALSE),
              sprintf("results/LOCO/z_%s.tsv", g), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(data.frame(sample=colnames(z), group=grp),
              sprintf("results/LOCO/lab_%s.tsv", g), sep="\t", row.names=FALSE, quote=FALSE)
  cat(sprintf("%-11s %5d genes x %3d samples (Tumour %d / Normal %d)\n", g, nrow(z), ncol(z),
      sum(grp=="Tumor"), sum(grp=="Normal")))
}
cat("DONE\n")
