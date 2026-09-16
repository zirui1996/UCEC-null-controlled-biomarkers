#!/usr/bin/env Rscript
# 00_download_GEO.R -- Phase 1, Step 1
# download the raw CEL files (RAW.tar) of the three discovery cohorts
# protocol: UCEC_bioinformatics_replication_protocol_v2.1.md section 7.1
# note: GSE36389 was removed from this study per section 0.1 and is not downloaded

proj <- Sys.getenv("UCEC_PROJ", unset = getwd())
setwd(proj)

series <- list(
  GSE17025  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17025/suppl/GSE17025_RAW.tar",
  GSE63678  = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63678/suppl/GSE63678_RAW.tar",
  GSE115810 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115810/suppl/GSE115810_RAW.tar"
)

for (gse in names(series)) {
  dest <- file.path("data_raw/GEO", gse, paste0(gse, "_RAW.tar"))
  if (file.exists(dest) && file.info(dest)$size > 1e6) {
    cat("[skip]", gse, "already present\n"); next
  }
  cat("[download]", gse, "->", dest, "\n")
  # -C - allow resuming an interrupted download
  st <- system2("curl", c("-L", "-C", "-", "--retry", "5", "--retry-delay", "5",
                          "-o", shQuote(dest), shQuote(series[[gse]])))
  cat("  exit:", st, " size:", file.info(dest)$size, "bytes\n")
}

# unpack the CEL files
for (gse in names(series)) {
  d   <- file.path("data_raw/GEO", gse)
  tar <- file.path(d, paste0(gse, "_RAW.tar"))
  cel <- file.path(d, "CEL")
  dir.create(cel, showWarnings = FALSE)
  if (length(list.files(cel, pattern = "\\.CEL", ignore.case = TRUE)) > 0) {
    cat("[skip untar]", gse, "\n"); next
  }
  cat("[untar]", gse, "\n")
  untar(tar, exdir = cel)
  # keep only CEL.gz, discard vendor CHP files (protocol §7.1)
  chp <- list.files(cel, pattern = "\\.CHP", ignore.case = TRUE, full.names = TRUE)
  if (length(chp)) { file.remove(chp); cat("  removed", length(chp), "CHP files\n") }
  cat("  CEL files:", length(list.files(cel, pattern = "\\.CEL", ignore.case = TRUE)), "\n")
}

cat("\n===== DOWNLOAD DONE =====\n")
for (gse in names(series)) {
  n <- length(list.files(file.path("data_raw/GEO", gse, "CEL"),
                         pattern = "\\.CEL", ignore.case = TRUE))
  cat(sprintf("%-12s %3d CEL\n", gse, n))
}
