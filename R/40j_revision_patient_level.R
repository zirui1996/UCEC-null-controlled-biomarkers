#!/usr/bin/env Rscript
# 40j_revision_patient_level.R -- patient as the unit of analysis in TCGA-UCEC
# The 512 primary-tumour RNA-seq samples come from 504 patients: four patients contributed 2-3 replicate aliquots.
# The survival table built in R/16 joined clinical data to samples, so these four patients entered the Cox models three
# times each (508 rows, 500 patients; 80 events of which 2 duplicated). This script:
#   (1) reruns the survival analysis with one record per patient (replicate samples averaged), replacing Table 4 / Table S18
#   (2) repeats the four-gene discrimination analysis with one tumour per patient as a sensitivity analysis
#   (3) gives the events needed / detectable hazard ratio for 80% power (Hsieh & Lavori 2000)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(survival); library(pROC)})
set.seed(42)
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")

se  <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop = FALSE]); code <- substr(cdl$barcode, 14, 15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code == "01" & cdl$primary_diagnosis %in% c(EN, SR)) | code == "11"
cnt <- assay(se, "unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]; cnt <- cnt[i, ]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt) != "", ]
bc <- cdl$barcode[keep]; pid <- substr(bc, 1, 12)
grp <- factor(ifelse(code[keep] == "11", "N", "T")); isT <- grp == "T"
d <- DGEList(cnt); d <- d[filterByExpr(d, group = grp), , keep.lib.sizes = FALSE]; d <- calcNormFactors(d)
lc <- cpm(d, log = TRUE); z <- t(scale(t(lc)))
dup <- names(which(table(pid[isT]) > 1))
cat(sprintf("tumour samples %d from %d patients; patients with replicates: %s\n", sum(isT), length(unique(pid[isT])),
            paste(sprintf("%s (%d)", dup, table(pid[isT])[dup]), collapse = ", ")))

## ---------------------------------------------------------------- (1) survival, one record per patient
sv0 <- read.delim("results/supp/survival_core4_per_patient.tsv", check.names = FALSE)
tz  <- z[, isT]; tp <- pid[isT]
agg <- t(sapply(split(seq_along(tp), tp), function(ix) rowMeans(tz[CORE, ix, drop = FALSE])))   # patient x gene, mean of replicates
cl  <- sv0[!duplicated(sv0$pid), c("pid", "OS_time", "OS_event")]
stopifnot(all(tapply(sv0$OS_time, sv0$pid, function(v) length(unique(v))) == 1))      # replicate rows share follow-up
sv  <- merge(cl, data.frame(pid = rownames(agg), agg, check.names = FALSE), by = "pid")
sv$panel <- -rowMeans(sv[, CORE])
cat(sprintf("survival: %d patients, %d deaths (previously %d rows, %d deaths)\n", nrow(sv), sum(sv$OS_event), nrow(sv0), sum(sv0$OS_event)))
rk <- summary(survfit(Surv(OS_time, 1 - OS_event) ~ 1, data = sv))$table
fu <- quantile(sv$OS_time[sv$OS_event == 0], c(.25, .5, .75))
cens <- data.frame(n = nrow(sv), deaths = sum(sv$OS_event), censored = sum(!sv$OS_event),
                   reverse_KM_median_days = unname(rk["median"]),
                   censored_followup_q25 = fu[1], censored_followup_median = fu[2], censored_followup_q75 = fu[3],
                   censored_before_1y = sum(!sv$OS_event & sv$OS_time < 365), censored_before_3y = sum(!sv$OS_event & sv$OS_time < 3 * 365),
                   deaths_within_3y = sum(sv$OS_event & sv$OS_time <= 3 * 365), at_risk_5y = sum(sv$OS_time >= 5 * 365))
res <- do.call(rbind, lapply(c(CORE, "panel"), function(v) {
  sv$x <- as.numeric(scale(sv[[v]])); m <- coxph(Surv(OS_time, OS_event) ~ x, data = sv); s <- summary(m)
  zp <- cox.zph(m, transform = "km")
  sv$hi <- factor(ifelse(sv[[v]] > median(sv[[v]]), "High", "Low"), levels = c("Low", "High"))
  lr <- survdiff(Surv(OS_time, OS_event) ~ hi, data = sv)
  data.frame(variable = v, events = m$nevent, HR_per_SD = s$coef[, "exp(coef)"], CI_low = s$conf.int[, 3], CI_high = s$conf.int[, 4],
             P = s$coef[, "Pr(>|z|)"], C_index = unname(s$concordance[1]), PH_chisq = zp$table["x", "chisq"], PH_P = zp$table["x", "p"],
             logrank_P_median_split_visual_only = 1 - pchisq(lr$chisq, 1))
}))
res$BH_q <- p.adjust(res$P, "BH")
print(cens, row.names = FALSE); print(res, row.names = FALSE, digits = 3)
# detectable hazard ratio per SD with 80% power, two-sided alpha 0.05 (Hsieh & Lavori 2000): D = (z_a + z_b)^2 / log(HR)^2
lhr <- (qnorm(0.975) + qnorm(0.80)) / sqrt(sum(sv$OS_event))
cat(sprintf("with %d deaths, 80%% power at alpha 0.05 for HR per SD of %.2f (or %.2f)\n", sum(sv$OS_event), exp(lhr), exp(-lhr)))
write.table(cens, "results/revision/survival_censoring.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(res, "results/revision/survival_cox_PH.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sv[, c("pid", "OS_time", "OS_event", CORE, "panel")], "results/revision/survival_per_patient_deduplicated.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(sprintf("deaths %d; detectable HR per SD %.2f / %.2f", sum(sv$OS_event), exp(lhr), exp(-lhr)), "results/revision/survival_power.txt")

## ---------------------------------------------------------------- (2) discrimination with one tumour per patient
fast_auc <- function(x, isT) { r <- rank(x); nT <- sum(isT); nN <- sum(!isT); (sum(r[isT]) - nT * (nT + 1) / 2) / (nT * nN) }
keep_mean <- rep(TRUE, length(pid))
# mean of replicate tumours: replace each patient's tumour columns by their mean, drop the extras
zt <- z; for (p in dup) { ix <- which(isT & pid == p); zt[, ix[1]] <- rowMeans(z[, ix]); keep_mean[ix[-1]] <- FALSE }
sens <- data.frame(
  analysis = c("all 512 tumour samples (primary)", "one tumour per patient: mean of replicates"),
  n_tumour = c(sum(isT), sum(isT & keep_mean)),
  panel_AUC = c(fast_auc(-colMeans(z[CORE, ]), isT), fast_auc(-colMeans(zt[CORE, keep_mean]), isT[keep_mean])),
  NUMA1_AUC = c({ a <- fast_auc(lc["NUMA1", ], isT); max(a, 1 - a) }, { a <- fast_auc(zt["NUMA1", keep_mean], isT[keep_mean]); max(a, 1 - a) }))
print(sens, row.names = FALSE, digits = 4)
write.table(sens, "results/revision/TCGA_one_tumour_per_patient.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("===== 40j DONE =====\n")
