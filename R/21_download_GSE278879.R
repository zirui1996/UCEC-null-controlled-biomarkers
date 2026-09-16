#!/usr/bin/env Rscript
# 21_download_GSE278879.R -- per-sample download with exact byte-size verification (incomplete files are not unpacked)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
d <- "data_raw/GEO/GSE278879"
FILES <- c(
  GSM8556535_HEEC1=206781869, GSM8556536_HEEC2=145990315, GSM8556537_HEEC3=46793998,
  GSM8556538_MEEC1=16304393,  GSM8556539_MEEC2=38582548,
  GSM8556540_LEEC1=29186006,  GSM8556541_LEEC2=40539574,
  GSM8556542_NE1=66352823,    GSM8556543_NE2=152202608,
  GSM8556544_NE3=95469709,    GSM8556545_NE4=95420418)
for (nm in names(FILES)) {
  gsm <- sub("_.*","",nm); f <- sprintf("%s_matrix.tar.gz", nm); dest <- file.path(d, f)
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%snnn/%s/suppl/%s", substr(gsm,1,7), gsm, f)
  for (a in 1:8) {
    sz <- if (file.exists(dest)) file.info(dest)$size else 0
    if (sz == FILES[[nm]]) break
    system2("curl", c("-L","-C","-","--retry","5","--retry-delay","5",
                      "--speed-time","60","--speed-limit","10000","-s","-o", shQuote(dest), shQuote(url)))
  }
  sz <- file.info(dest)$size
  if (sz != FILES[[nm]]) { cat(sprintf("!! %s incomplete download %d / %d, not unpacking\n", f, sz, FILES[[nm]])); next }
  cat(sprintf("%-22s %6.1f MB  OK\n", nm, sz/1024/1024))
  od <- file.path(d, nm); dir.create(od, showWarnings=FALSE)
  untar(dest, exdir=od)
}
cat("\n=== unpacking results ===\n")
for (nm in names(FILES)) {
  fs <- list.files(file.path(d,nm), recursive=TRUE)
  cat(sprintf("%-22s %s\n", nm, if(length(fs)) paste(basename(fs), collapse=", ") else "EMPTY"))
}
cat("\n===== download complete =====\n")
