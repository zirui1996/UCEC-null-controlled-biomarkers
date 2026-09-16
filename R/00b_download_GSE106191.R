#!/usr/bin/env Rscript
# 00b_download_GSE106191.R -- replacement cohort for Discovery 3
# v2: add an integrity check; a truncated tar must never be unpacked
proj <- Sys.getenv("UCEC_PROJ", unset = getwd()); setwd(proj)
url      <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106191/suppl/GSE106191_RAW.tar"
dest     <- "data_raw/GEO/GSE106191/GSE106191_RAW.tar"
EXPECTED <- 471060480L          # from the GEO filelist.txt

for (attempt in 1:10) {
  sz <- if (file.exists(dest)) file.info(dest)$size else 0
  if (sz == EXPECTED) break
  cat(sprintf("[attempt %d] %.1f MB of %.1f MB so far, resuming...\n",
              attempt, sz/1e6, EXPECTED/1e6))
  system2("curl", c("-L","-C","-","--retry","5","--retry-delay","5",
                    "--speed-time","60","--speed-limit","10000","-s",
                    "-o", shQuote(dest), shQuote(url)))
}
sz <- file.info(dest)$size
if (sz != EXPECTED) stop(sprintf("incomplete download: %d / %d bytes, not unpacked", sz, EXPECTED))
cat(sprintf("[ok] tar is complete: %d bytes\n", sz))

cel <- "data_raw/GEO/GSE106191/CEL"
dir.create(cel, showWarnings = FALSE)
untar(dest, exdir = cel)
chp <- list.files(cel, pattern="\\.CHP", ignore.case=TRUE, full.names=TRUE)
if (length(chp)) { file.remove(chp); cat("  removed", length(chp), "CHP\n") }
n <- length(list.files(cel, pattern="\\.CEL", ignore.case=TRUE))
cat("  CEL files:", n, "\n")
if (n != 97) warning("the number of CEL files is not the expected 97")
cat("===== GSE106191 DOWNLOAD DONE =====\n")
