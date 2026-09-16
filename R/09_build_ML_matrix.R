#!/usr/bin/env Rscript
# 09_build_ML_matrix.R -- build the machine-learning input matrix of section 12
# X = recurrent consensus genes (v2.6: 182); the three voting cohorts are pooled after z-scoring within each dataset (same rule as section 14.2)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
dir.create("results/ML", showWarnings=FALSE, recursive=TRUE)

rec <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors=FALSE)$gene
MAT <- c(GSE17025 ="data_processed/expression/GSE17025_gene_matrix_bestprobe.tsv",
         GSE63678 ="data_processed/expression/GSE63678_gene_matrix_bestprobe.tsv",
         GSE106191="data_processed/expression/GSE106191_all97_gene_matrix_bestprobe.tsv")
man <- read.delim("metadata/sample_manifest.tsv", stringsAsFactors=FALSE)

blocks <- list()
for (g in names(MAT)) {
  x <- read.delim(MAT[g], row.names=1, check.names=FALSE)
  x <- x[intersect(rec, rownames(x)), , drop=FALSE]
  z <- t(scale(t(as.matrix(x))))                     # z-score within the dataset
  if (g=="GSE106191") {
    m <- read.delim("metadata/sample_manifest_GSE106191.tsv", stringsAsFactors=FALSE)
    grp <- setNames(ifelse(m$normal_or_tumor=="tumor","Tumor","Normal"), m$sample_id)
  } else {
    m <- subset(man, dataset==g & include_primary=="yes")
    grp <- setNames(ifelse(m$normal_or_tumor=="tumor","Tumor","Normal"), m$sample_id)
  }
  grp <- grp[colnames(z)]
  stopifnot(!any(is.na(grp)))
  cat(sprintf("%-11s genes %3d / %d | samples %3d (Tumour %d / Normal %d)\n",
      g, nrow(z), length(rec), ncol(z), sum(grp=="Tumor"), sum(grp=="Normal")))
  blocks[[g]] <- list(z=z, grp=grp)
}

genes <- Reduce(intersect, lapply(blocks, function(b) rownames(b$z)))
cat("\ngenes measurable in all three cohorts:", length(genes), "/", length(rec), "\n")
X <- do.call(cbind, lapply(blocks, function(b) b$z[genes,,drop=FALSE]))
y <- unlist(lapply(blocks, function(b) b$grp), use.names=FALSE)
ds <- rep(names(blocks), sapply(blocks, function(b) ncol(b$z)))
cat("after merging:", ncol(X), "samples x", nrow(X), "genes | Tumour", sum(y=="Tumor"), "/ Normal", sum(y=="Normal"), "\n")

out <- data.frame(sample=colnames(X), dataset=ds, group=y, t(X), check.names=FALSE)
write.table(out, "results/ML/ML_input_matrix.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("wrote results/ML/ML_input_matrix.tsv\n")
