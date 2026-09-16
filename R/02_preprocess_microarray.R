#!/usr/bin/env Rscript
# 02_preprocess_microarray.R -- Phase 1, Step 1/3
# per-cohort justRMA followed by probe -> gene symbol collapse
# protocol section 7.1 (affy::justRMA; not oligo, not vendor CHP) and section 9.1 (best-probe as the main rule, mean as a sensitivity analysis)
# usage: Rscript R/02_preprocess_microarray.R GSE17025
#   process one cohort at a time to limit peak memory use on an 8 GB machine

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript 02_preprocess_microarray.R <GSE_ID>")
gse  <- args[1]
proj <- Sys.getenv("UCEC_PROJ", unset = getwd()); setwd(proj)

suppressPackageStartupMessages({ library(affy) })

ANN <- c(GSE17025 = "hgu133plus2.db", GSE63678 = "hgu133a2.db", GSE115810 = "hgu133a.db")
stopifnot(gse %in% names(ANN))
annpkg <- ANN[[gse]]
suppressPackageStartupMessages(library(annpkg, character.only = TRUE))

man <- read.delim("metadata/sample_manifest.tsv", stringsAsFactors = FALSE, quote = "")
man <- man[man$dataset == gse & man$include_primary == "yes", ]
stopifnot(nrow(man) > 0)

celdir <- file.path("data_raw/GEO", gse, "CEL")
allcel <- list.files(celdir, pattern = "\\.CEL(\\.gz)?$", ignore.case = TRUE)

# match exactly on source_file from the manifest; inferring groups from file order is forbidden
man$cel <- man$source_file
miss <- setdiff(man$cel, allcel)
if (length(miss)) stop("CEL files listed in the manifest are missing from disk: ", paste(head(miss), collapse = ", "))
extra <- setdiff(allcel, man$cel)
if (length(extra)) cat("[note]", length(extra), "CEL files are not in the main analysis manifest and were ignored\n")

cat("[", gse, "] justRMA on", nrow(man), "arrays ...\n")
t0 <- Sys.time()
eset <- justRMA(filenames = man$cel, celfile.path = celdir, verbose = FALSE)
cat("  elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 2), "min\n")

expr <- exprs(eset)                      # already log2
colnames(expr) <- man$sample_id[match(colnames(expr), man$cel)]
stopifnot(!any(is.na(colnames(expr))))

dir.create("data_processed/expression", recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(probe = rownames(expr), expr, check.names = FALSE),
            sprintf("data_processed/expression/%s_RMA_probe_matrix.tsv", gse),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  probe matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")

# ---- probe -> symbol ----
sym <- AnnotationDbi::mapIds(get(annpkg), keys = rownames(expr),
                             column = "SYMBOL", keytype = "PROBEID",
                             multiVals = "first")
keep <- !is.na(sym)
cat("  mapped probes:", sum(keep), "/", length(sym),
    " unique symbols:", length(unique(sym[keep])), "\n")
e <- expr[keep, , drop = FALSE]; g <- sym[keep]

collapse <- function(mat, genes, how) {
  sp <- split(seq_len(nrow(mat)), genes)
  out <- t(vapply(sp, function(ix) {
    if (length(ix) == 1) return(mat[ix, ])
    if (how == "best") mat[ix[which.max(rowMeans(mat[ix, , drop = FALSE]))], ]
    else colMeans(mat[ix, , drop = FALSE])
  }, numeric(ncol(mat))))
  colnames(out) <- colnames(mat); out
}

for (how in c("best", "mean")) {
  gm <- collapse(e, g, how)
  f <- sprintf("data_processed/expression/%s_gene_matrix_%sprobe.tsv", gse, how)
  write.table(data.frame(gene = rownames(gm), gm, check.names = FALSE),
              f, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  ", how, "-probe gene matrix:", nrow(gm), "genes ->", basename(f), "\n")
}

writeLines(capture.output(sessionInfo()),
           sprintf("logs/02_preprocess_%s_sessionInfo.txt", gse))
cat("[", gse, "] DONE\n")
