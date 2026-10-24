#!/usr/bin/env Rscript
# SUPERSEDED by R/40j: this version used the per-sample table of R/16 (508 rows, 500 patients; four patients counted three times).
# Kept for traceability; outputs carry the suffix _508rows_superseded.
# 40d_revision_survival.R -- R1 minor 8: event counts, censoring pattern and proportional-hazards checks for the
# univariable Cox models of the four core genes and the panel score (continuous, per SD). Median splits are kept
# for visualization only.
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages(library(survival))
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")
sv <- read.delim("results/supp/survival_core4_per_patient.tsv", check.names = FALSE)
cat(sprintf("n = %d | deaths = %d | censored = %d (%.1f%%)\n", nrow(sv), sum(sv$OS_event), sum(!sv$OS_event), 100 * mean(!sv$OS_event)))
rk <- summary(survfit(Surv(OS_time, 1 - OS_event) ~ 1, data = sv))$table
fu <- quantile(sv$OS_time[sv$OS_event == 0], c(.25, .5, .75))
cens <- data.frame(n = nrow(sv), deaths = sum(sv$OS_event), censored = sum(!sv$OS_event),
                   reverse_KM_median_days = unname(rk["median"]),
                   censored_followup_q25 = fu[1], censored_followup_median = fu[2], censored_followup_q75 = fu[3],
                   censored_before_1y = sum(!sv$OS_event & sv$OS_time < 365),
                   censored_before_3y = sum(!sv$OS_event & sv$OS_time < 3 * 365),
                   deaths_within_3y = sum(sv$OS_event & sv$OS_time <= 3 * 365),
                   at_risk_5y = sum(sv$OS_time >= 5 * 365))
print(cens, row.names = FALSE)

res <- do.call(rbind, lapply(c(CORE, "panel"), function(v) {
  sv$x <- as.numeric(scale(sv[[v]])); m <- coxph(Surv(OS_time, OS_event) ~ x, data = sv); s <- summary(m)
  zp <- cox.zph(m, transform = "km")
  data.frame(variable = v, events = m$nevent, HR_per_SD = s$coef[, "exp(coef)"], CI_low = s$conf.int[, 3], CI_high = s$conf.int[, 4],
             P = s$coef[, "Pr(>|z|)"], C_index = unname(s$concordance[1]),
             PH_chisq = zp$table["x", "chisq"], PH_P = zp$table["x", "p"])
}))
res$BH_q <- p.adjust(res$P, "BH")
print(res, row.names = FALSE, digits = 3)

pdf("results/revision/survival_schoenfeld_508rows_superseded.pdf", width = 9, height = 6); par(mfrow = c(2, 3))
for (v in c(CORE, "panel")) { sv$x <- as.numeric(scale(sv[[v]]))
  plot(cox.zph(coxph(Surv(OS_time, OS_event) ~ x, data = sv), transform = "km"), main = v, ylab = "Beta(t) for score per SD") ; abline(h = 0, lty = 3) }
dev.off()
write.table(cens, "results/revision/survival_censoring_508rows_superseded.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(res, "results/revision/survival_cox_PH_508rows_superseded.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("===== 40d DONE =====\n")
