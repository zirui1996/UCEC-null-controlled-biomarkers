#!/usr/bin/env Rscript
# 40h_revision_sample_inclusion.R -- R1 minor 9: sample-level inclusion list for every screened GEO cohort
# (Online Resource 2). One row per GSM: cohort, role, group, whether it entered the analysis, and why not.
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
m  <- read.delim("metadata/sample_manifest.tsv", stringsAsFactors = FALSE, quote = "")
m6 <- read.delim("metadata/sample_manifest_GSE106191.tsv", stringsAsFactors = FALSE, quote = "")
m9 <- read.delim("metadata/sample_manifest_GSE146889.tsv", stringsAsFactors = FALSE)
role <- c(GSE17025 = "Discovery", GSE63678 = "Discovery", GSE106191 = "Discovery", GSE56087 = "Not admitted (supplementary direction check)",
          GSE115810 = "Excluded (failed admission check)", GSE36389 = "Excluded (duplicate raw data of GSE115810 with conflicting labels)",
          GSE146889 = "Held-out replication", GSE278879 = "Cellular localization (single-nucleus)")
lab <- function(x) ifelse(x == "tumor", "EC", ifelse(x == "normal", "Control", x))
a <- subset(m, dataset %in% c("GSE17025","GSE63678","GSE115810","GSE36389","GSE278879"))
a <- data.frame(cohort = a$dataset, sample_id = a$sample_id, platform = a$platform, title = a$title, group = lab(a$normal_or_tumor),
                analyzed = ifelse(a$include_primary %in% c("yes","single_cell"), "yes", "no"), reason_not_analyzed = a$exclusion_reason)
a$reason_not_analyzed[a$cohort == "GSE115810"] <- "cohort failed the label-permutation admission check"
b <- data.frame(cohort = "GSE106191", sample_id = m6$sample_id, platform = m6$platform, title = m6$title,
                group = ifelse(m6$normal_or_tumor == "tumor", "EC", "Control (hyperplasia)"), analyzed = "yes", reason_not_analyzed = "")
b$note <- paste("collection period", m6$era)
gsm56 <- c(GSM1354841="normal",GSM1354842="tumor",GSM1354843="normal",GSM1354844="tumor",GSM1354845="normal",GSM1354846="tumor",
           GSM1354847="normal",GSM1354848="tumor",GSM1354849="normal",GSM1354850="tumor",GSM1354851="tumor",GSM1354852="normal",
           GSM1354853="normal",GSM1354854="tumor",GSM1354855="normal",GSM1354856="tumor",GSM1354857="normal",GSM1354858="tumor")
c56 <- data.frame(cohort = "GSE56087", sample_id = names(gsm56), platform = "GPL11154", title = "", group = lab(unname(gsm56)),
                  analyzed = "yes (direction check only)", reason_not_analyzed = "")
d9 <- data.frame(cohort = "GSE146889", sample_id = m9$sample_id, platform = m9$platform, title = m9$title, group = lab(m9$group),
                 analyzed = ifelse(m9$tissue %in% "Endometrial", "yes", "no"),
                 reason_not_analyzed = ifelse(m9$tissue %in% "Endometrial", "", ifelse(is.na(m9$tissue), "tissue of origin not annotated as endometrium", paste("not endometrium:", m9$tissue))))
d9$note <- ifelse(m9$tissue %in% "Endometrial", paste("pair", m9$pair_id, "| MSI status", m9$msi), "")
a$note <- ""; c56$note <- ""
out <- rbind(a, b, c56, d9); out$role <- unname(role[out$cohort])
out <- out[order(match(out$cohort, names(role)), out$sample_id), c("cohort","role","sample_id","platform","title","group","analyzed","reason_not_analyzed","note")]
write.table(out, "results/revision/sample_inclusion_list.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(table(out$cohort, out$analyzed))
cat("===== 40h DONE =====\n")
