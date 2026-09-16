#!/usr/bin/env Rscript
# 17_selfcheck_independent.R -- independent self-check of the section 16 results
# purpose: confirm that the negative NUMA1 overall-survival result is not an artefact of our processing or clinical annotation.
# recompute after changing two things:
#   expression  GDC STAR-counts (TMM-logCPM)  ->  UCSC Xena TCGA-UCEC HiSeqV2 (RNA-seqV2 RSEM, an independent quantification pipeline)
#   survival  GDC clinical fields        ->  TCGA-CDR (Liu et al. 2018 Cell), the curated survival standard
# this is not external validation (the patients are the same); it is a reproducibility self-check.
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(survival)})
GX <- "NUMA1"

x <- read.delim(gzfile("data_raw/Xena/TCGA_UCEC_HiSeqV2.gz"), row.names=1, check.names=FALSE)
cat("Xena matrix:", nrow(x), "genes x", ncol(x), "samples\n")
code <- substr(colnames(x), 14, 15)
cat("sample types: 01 (primary tumour)", sum(code=="01"), " 11 (adjacent normal)", sum(code=="11"), "\n")
xt <- x[, code=="01", drop=FALSE]
pid <- substr(colnames(xt), 1, 12)

cdr <- read.delim("data_raw/Xena/TCGA_CDR_Liu2018.tsv", check.names=FALSE, stringsAsFactors=FALSE)
cdr <- cdr[cdr[["cancer type abbreviation"]]=="UCEC", ]
cat("UCEC patients in the CDR:", nrow(cdr), "\n")
cdr$pid <- cdr[["_PATIENT"]]
m <- match(pid, cdr$pid)
cat("expression-matrix patients matched to the CDR:", sum(!is.na(m)), "/", length(pid), "\n\n")

d <- data.frame(pid=pid, expr=as.numeric(xt[GX, ]),
                OS=as.numeric(cdr$OS[m]), OS_time=as.numeric(cdr$OS.time[m]),
                age=as.numeric(cdr$age_at_initial_pathologic_diagnosis[m]),
                stage=cdr$clinical_stage[m], grade=cdr$histological_grade[m],
                histo=cdr$histological_type[m])
d <- d[!is.na(d$expr) & !is.na(d$OS) & !is.na(d$OS_time) & d$OS_time>0 & !duplicated(d$pid), ]
d$z <- as.numeric(scale(d$expr))
cat(sprintf("self-check analysis set: n=%d, OS events=%d, median follow-up %.0f days\n", nrow(d), sum(d$OS), median(d$OS_time)))

cat("\n=== univariable Cox (same definition as section 16.2) ===\n")
u <- coxph(Surv(OS_time, OS) ~ z, data=d); s <- summary(u)
cat(sprintf("  HR = %.3f (95%%CI %.3f-%.3f), P = %.3f | C-index = %.3f\n",
    s$coef[,"exp(coef)"], s$conf.int[,3], s$conf.int[,4], s$coef[,"Pr(>|z|)"], s$concordance[1]))

cat("\n=== Kaplan-Meier (median cutoff, same definition as section 16.5) ===\n")
d$g <- ifelse(d$z > median(d$z), "High", "Low")
lr <- survdiff(Surv(OS_time, OS) ~ g, data=d)
cat(sprintf("  log-rank P = %.3f\n", 1-pchisq(lr$chisq, length(lr$n)-1)))

cat("\n=== control: what an automatic optimal-cutpoint scan, as used by default in KM-Plotter, would give ===\n")
qs <- seq(0.2, 0.8, by=0.01)
ps <- sapply(qs, function(q){ g <- d$z > quantile(d$z, q)
  if (length(unique(g))<2) return(NA)
  lrx <- survdiff(Surv(OS_time,OS) ~ g, data=d); 1-pchisq(lrx$chisq, 1) })
best <- which.min(ps)
cat(sprintf("  scanned %d cutpoints (20th-80th percentile): smallest P = %.4f (at the %.0fth percentile)\n",
    length(qs), min(ps,na.rm=TRUE), 100*qs[best]))
cat(sprintf("  cutpoints reaching P<0.05: %d / %d (%.0f%%)\n", sum(ps<0.05,na.rm=TRUE), sum(!is.na(ps)),
    100*mean(ps<0.05,na.rm=TRUE)))
if (sum(ps<0.05,na.rm=TRUE)==0) {
  cat("  -> OBSERVED: not even an optimal-cutpoint scan yields any significance, so the NUMA1 OS result is negative\n")
  cat("     so robust that cutpoint hacking cannot rescue it. KM-Plotter automatic cutpoints will also be negative.\n")
} else {
  cat("  -> an optimal cutpoint does yield significance, but this is a multiple-testing artefact, so the protocol fixes the median cutpoint.\n")
}

pdf("results/survival/figures/selfcheck_cutpoint_scan.pdf", width=6, height=4.5)
plot(100*qs, ps, type="l", lwd=2, xlab="Cutoff percentile", ylab="log-rank P",
     main=sprintf("%s — cutoff scan (TCGA-UCEC, Xena+CDR)", GX), ylim=c(0,1))
abline(h=0.05, lty=2, col="red"); abline(v=50, lty=3, col="blue")
legend("topright", c("P = 0.05","median (protocol)"), lty=c(2,3), col=c("red","blue"), bty="n"); dev.off()

writeLines(capture.output(sessionInfo()), "logs/17_selfcheck_sessionInfo.txt")
cat("\n===== 17_selfcheck DONE =====\n")
