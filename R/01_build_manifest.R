#!/usr/bin/env Rscript
# 01_build_manifest.R -- Phase 1, Step 0
# build the approved analysis manifest from the raw manifest scraped from GEO
# protocol sections 6.1 / 6.2 / 6.3

proj <- Sys.getenv("UCEC_PROJ", unset = getwd()); setwd(proj)
raw <- read.delim("metadata/sample_manifest_raw.tsv", stringsAsFactors = FALSE,
                  check.names = FALSE, quote = "")

# ---- extract fields from characteristics ----
getf <- function(x, key) {
  out <- rep(NA_character_, length(x))
  for (i in seq_along(x)) {
    parts <- trimws(strsplit(x[i], "\\|")[[1]])
    hit <- grep(paste0("^", key, "\\s*:"), parts, value = TRUE, ignore.case = TRUE)
    if (length(hit)) out[i] <- trimws(sub("^[^:]+:\\s*", "", hit[1]))
  }
  out
}

m <- data.frame(
  sample_id = raw$sample_id, dataset = raw$dataset, platform = raw$platform,
  title = raw$title, stringsAsFactors = FALSE
)

m$histology <- NA_character_; m$grade <- NA_character_
m$stage <- NA_character_;     m$age   <- NA_character_

i <- raw$dataset == "GSE17025"
m$histology[i] <- getf(raw$characteristics[i], "histology")
m$grade[i]     <- getf(raw$characteristics[i], "grade")
m$stage[i]     <- getf(raw$characteristics[i], "stage")
m$age[i]       <- getf(raw$characteristics[i], "age \\(years\\)")

i <- raw$dataset == "GSE115810"
m$grade[i] <- getf(raw$characteristics[i], "grade")
m$histology[i] <- "endometrioid (not specified in GEO)"

i <- raw$dataset == "GSE21882"
m$grade[i] <- getf(raw$characteristics[i], "differentation grade")
m$stage[i] <- getf(raw$characteristics[i], "surgical stage")
m$age[i]   <- getf(raw$characteristics[i], "age at diagnosis")
m$histology[i] <- "type I endometrioid"

i <- raw$dataset == "GSE278879"
m$grade[i] <- getf(raw$characteristics[i], "grade")
m$age[i]   <- getf(raw$characteristics[i], "age")

# ---- grouping ----
m$group           <- raw$group
m$normal_or_tumor <- ifelse(raw$group %in% c("Tumor","HEEC","MEEC","LEEC",
                                             "Survivor","NonSurvivor"), "tumor",
                     ifelse(raw$group %in% c("Normal","NE"), "normal", NA))

# ---- survival fields: confirmed unavailable in GEO, see protocol section 17.1 ----
m$survival_time  <- NA
m$survival_event <- NA

m$source_file <- raw$raw_file

# ---- inclusion/exclusion (protocol section 6.2) ----
m$include_primary   <- "no"
m$exclusion_reason  <- NA_character_

disc <- raw$dataset %in% c("GSE17025","GSE63678","GSE115810") & raw$group %in% c("Tumor","Normal")
m$include_primary[disc] <- "yes"

m$exclusion_reason[raw$dataset == "GSE63678" & raw$group == "EXCLUDE"] <-
  "not endometrium (cervix/vulva) - protocol 4.2"
m$exclusion_reason[raw$dataset == "GSE36389"] <-
  "duplicate raw CEL of GSE115810 with conflicting phenotype labels - protocol 0.1"
m$include_primary[raw$dataset == "GSE278879"] <- "single_cell"
m$include_primary[raw$dataset == "GSE21882"]  <- "pending_mapping"
m$exclusion_reason[raw$dataset == "GSE21882"] <-
  "endpoint definition must be reconciled with PMID 20635577 before use - protocol 17.2"

write.table(m, "metadata/sample_manifest.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")

# ---- exclusion_log ----
ex <- m[!is.na(m$exclusion_reason), c("sample_id","dataset","exclusion_reason")]
names(ex)[3] <- "reason"
ex$decision_date  <- Sys.Date()
ex$decision_stage <- "Step 0 (pre-analysis)"
write.table(ex, "metadata/exclusion_log.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

# ---- dataset_manifest ----
ds <- do.call(rbind, lapply(split(m, m$dataset), function(d) {
  data.frame(dataset = d$dataset[1], platform = d$platform[1],
             n_total = nrow(d),
             n_tumor  = sum(d$normal_or_tumor == "tumor"  & d$include_primary != "no", na.rm = TRUE),
             n_normal = sum(d$normal_or_tumor == "normal" & d$include_primary != "no", na.rm = TRUE),
             n_excluded = sum(!is.na(d$exclusion_reason)),
             stringsAsFactors = FALSE)
}))
ds$role <- c(GSE115810 = "Discovery 3", GSE17025 = "Discovery 1", GSE21882 = "External prognosis (pending)",
             GSE278879 = "Single-cell localization", GSE36389 = "EXCLUDED (see protocol 0.1)",
             GSE63678 = "Discovery 2")[ds$dataset]
ds$download_date <- Sys.Date()
write.table(ds, "metadata/dataset_manifest.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("===== sample_manifest.tsv =====\n")
print(table(m$dataset, m$include_primary))
cat("\n===== dataset_manifest.tsv =====\n"); print(ds, row.names = FALSE)
cat("\n===== total included in the main analysis =====\n")
d <- m[m$include_primary == "yes", ]
print(table(d$dataset, d$normal_or_tumor))
