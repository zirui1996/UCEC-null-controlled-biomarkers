#!/usr/bin/env Rscript
# 02b_preprocess_GSE106191.R -- preprocessing of the Discovery 3 replacement cohort
# main analysis set = the 2016 scan batch (35 carcinomas / 19 hyperplasias); a second run over all 97 samples serves as the S-b comparison
proj <- Sys.getenv("UCEC_PROJ", unset = getwd()); setwd(proj)
suppressPackageStartupMessages({ library(affy); library(hgu133plus2.db) })

man <- read.delim("metadata/sample_manifest_GSE106191.tsv", stringsAsFactors=FALSE, quote="")
d   <- "data_raw/GEO/GSE106191/CEL"
stopifnot(all(file.exists(file.path(d, man$source_file))))

collapse <- function(mat, genes, how) {
  sp <- split(seq_len(nrow(mat)), genes)
  out <- t(vapply(sp, function(ix) {
    if (length(ix)==1) return(mat[ix,])
    if (how=="best") mat[ix[which.max(rowMeans(mat[ix,,drop=FALSE]))],] else colMeans(mat[ix,,drop=FALSE])
  }, numeric(ncol(mat))))
  colnames(out) <- colnames(mat); out
}

run <- function(sub, tag) {
  cat("\n==== ", tag, ": ", nrow(sub), " arrays (tumor ", sum(sub$normal_or_tumor=="tumor"),
      " / normal ", sum(sub$normal_or_tumor=="normal"), ") ====\n", sep="")
  t0 <- Sys.time()
  eset <- justRMA(filenames=sub$source_file, celfile.path=d, verbose=FALSE)
  cat("  justRMA:", round(difftime(Sys.time(),t0,units="mins"),2), "min\n")
  e <- exprs(eset)
  colnames(e) <- sub$sample_id[match(colnames(e), sub$source_file)]
  stopifnot(!any(is.na(colnames(e))))
  write.table(data.frame(probe=rownames(e), e, check.names=FALSE),
              sprintf("data_processed/expression/%s_RMA_probe_matrix.tsv", tag),
              sep="\t", row.names=FALSE, quote=FALSE)
  sym <- AnnotationDbi::mapIds(hgu133plus2.db, keys=rownames(e), column="SYMBOL",
                               keytype="PROBEID", multiVals="first")
  k <- !is.na(sym)
  cat("  mapped probes:", sum(k), "/", length(sym), " unique symbols:", length(unique(sym[k])), "\n")
  for (how in c("best","mean")) {
    gm <- collapse(e[k,,drop=FALSE], sym[k], how)
    write.table(data.frame(gene=rownames(gm), gm, check.names=FALSE),
                sprintf("data_processed/expression/%s_gene_matrix_%sprobe.tsv", tag, how),
                sep="\t", row.names=FALSE, quote=FALSE)
    cat("  ", how, "-probe:", nrow(gm), "genes\n")
  }
}

run(subset(man, include_primary=="yes"), "GSE106191")       # main analysis: the 2016 batch
run(man,                                  "GSE106191_all97") # S-b comparison: all samples
writeLines(capture.output(sessionInfo()), "logs/02b_preprocess_GSE106191_sessionInfo.txt")
cat("\n===== GSE106191 PREPROCESS DONE =====\n")
