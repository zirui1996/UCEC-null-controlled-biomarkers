#!/usr/bin/env Rscript
# 16_TCGA_prognosis.R -- protocol section 16: prognostic analysis in TCGA-UCEC
# all modelling decisions were fixed in logs/16_prognosis_prespecification.md before any hazard ratio was seen
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(survival); library(survminer); library(timeROC)})
dir.create("results/survival/figures", showWarnings=FALSE, recursive=TRUE)
set.seed(42)
GX <- "NUMA1"

o <- readRDS("results/validation/TCGA_geneX_objects.rds")
lc <- o$lc; cd <- o$cd; grp <- o$grp
t <- cd[grp=="Tumor", ]; expr <- lc[, grp=="Tumor", drop=FALSE]

d <- data.frame(
  OS_time  = ifelse(t$vital_status=="Dead", t$days_to_death, t$days_to_last_follow_up),
  OS_event = as.integer(t$vital_status=="Dead"),
  age      = as.numeric(t$age_at_index),
  stage    = factor(sub("^(Stage I{1,3}V?|Stage IV).*$","\\1", t$figo_stage), ),
  grade_r  = t$tumor_grade,
  histo    = factor(ifelse(grepl("Endometrioid", t$primary_diagnosis), "endometrioid","serous")),
  subtype_r= as.character(t$paper_IntegrativeCluster),
  grade_chr= as.character(t$tumor_grade), stringsAsFactors=FALSE)
d[[paste0(GX,"_z")]] <- as.numeric(scale(expr[GX, ]))
# prespecified coding
d$stage <- factor(ifelse(grepl("Stage IV", t$figo_stage), "IV",
             ifelse(grepl("Stage III", t$figo_stage), "III",
             ifelse(grepl("Stage II", t$figo_stage), "II",
             ifelse(grepl("Stage I", t$figo_stage), "I", NA)))), levels=c("I","II","III","IV"))
d$grade <- factor(ifelse(d$grade_chr %in% c("G3","High Grade"), "G3", d$grade_chr), levels=c("G1","G2","G3"))
d$subtype <- factor(ifelse(d$subtype_r %in% c("CN high","CN low","MSI","POLE"), d$subtype_r, NA),
                    levels=c("CN low","CN high","MSI","POLE"))
d <- d[!is.na(d$OS_time) & d$OS_time>0 & !is.na(d$OS_event), ]
cat(sprintf("analysis set: n=%d, events=%d, median follow-up %.0f days\n\n", nrow(d), sum(d$OS_event), median(d$OS_time)))

f <- function(x) sprintf("%.3f (%.3f-%.3f), P=%.3g", exp(x[1]), exp(x[2]), exp(x[3]), x[4])
show <- function(fit, nm) {
  s <- summary(fit); cat(sprintf("\n--- %s | n=%d, events=%d, EPV=%.1f ---\n",
      nm, s$n, s$nevent, s$nevent/length(coef(fit))))
  cf <- cbind(coef(s)[,"coef"], confint(fit), coef(s)[,"Pr(>|z|)"])
  for (i in rownames(cf)) cat(sprintf("  %-24s HR = %s\n", i, f(cf[i,])))
  cat(sprintf("  C-index = %.3f\n", s$concordance[1]))
  invisible(fit)
}

## ---- section 16.2 univariable Cox (primary) ----
u <- coxph(as.formula(sprintf("Surv(OS_time,OS_event) ~ %s_z", GX)), data=d)
show(u, sprintf("section 16.2 univariable Cox (%s per 1 SD)", GX))

## ---- section 16.3 Model 1: primary clinical model ----
m1 <- coxph(as.formula(sprintf("Surv(OS_time,OS_event) ~ %s_z + age + stage + grade + histo", GX)), data=d)
show(m1, "section 16.3 Model 1, primary clinical model")

## ---- section 16.3 Model 2: sensitivity model adjusted for molecular subtype ----
d2 <- d[!is.na(d$subtype), ]
cat(sprintf("\n[handling of missing subtypes] samples with a subtype %d / %d (%.1f%% missing), events %d\n",
    nrow(d2), nrow(d), 100*(1-nrow(d2)/nrow(d)), sum(d2$OS_event)))
m2 <- coxph(as.formula(sprintf("Surv(OS_time,OS_event) ~ %s_z + age + stage + grade + histo + subtype", GX)), data=d2)
show(m2, "section 16.3 Model 2, sensitivity model adjusted for subtype")
if (summary(m2)$nevent/length(coef(m2)) < 10)
  cat("  WARNING: EPV < 10, the model is underpowered; the HR is indicative only and supports no independent conclusion (§16.3)\n")

## ---- section 16.4 proportional-hazards assumption ----
cat("\n--- section 16.4 proportional-hazards test (Model 1) ---\n"); print(cox.zph(m1))
pdf("results/survival/figures/PH_assumption_Model1.pdf", width=9, height=7)
par(mfrow=c(2,3)); plot(cox.zph(m1)); dev.off()

## ---- section 16.5 Kaplan-Meier (cutoff fixed at the median; optimal cutpoints are forbidden) ----
d$grp_med <- factor(ifelse(d[[paste0(GX,"_z")]] > median(d[[paste0(GX,"_z")]]), "High","Low"), levels=c("Low","High"))
sf <- survfit(Surv(OS_time,OS_event) ~ grp_med, data=d)
lr <- survdiff(Surv(OS_time,OS_event) ~ grp_med, data=d)
cat(sprintf("\n--- §16.5 KM（median cutoff）--- log-rank P = %.3g\n",
    1-pchisq(lr$chisq, length(lr$n)-1)))
p <- ggsurvplot(sf, data=d, pval=TRUE, risk.table=TRUE, conf.int=TRUE,
     legend.labs=c(paste0(GX," Low"), paste0(GX," High")), xlab="Days", palette=c("#4393C3","#D6604D"),
     title=sprintf("TCGA-UCEC OS by %s (median cutoff)", GX))
pdf("results/survival/figures/KM_NUMA1_median.pdf", width=7, height=7, onefile=FALSE); print(p); dev.off()

## ---- §16.6 time-dependent ROC + C-index ----
tr <- timeROC(T=d$OS_time, delta=d$OS_event, marker=d[[paste0(GX,"_z")]], cause=1,
              weighting="marginal", times=c(365,1095,1825), iid=TRUE)
cat("\n--- §16.6 time-dependent ROC ---\n")
for (i in seq_along(tr$times)) cat(sprintf("  %d-year AUC = %.3f (SE %.3f)\n",
    round(tr$times[i]/365), tr$AUC[i], tr$inference$vect_sd_1[i]))
cat(sprintf("  Harrell C-index (univariable) = %.3f\n", summary(u)$concordance[1]))
pdf("results/survival/figures/timeROC_NUMA1.pdf", width=5.5, height=5.5)
plot(tr, time=365, col="#2166AC", lwd=2, title=FALSE)
plot(tr, time=1095, col="#B2182B", lwd=2, add=TRUE); plot(tr, time=1825, col="#1B7837", lwd=2, add=TRUE)
legend("bottomright", sprintf("%d-year AUC=%.3f", c(1,3,5), tr$AUC), col=c("#2166AC","#B2182B","#1B7837"), lwd=2, bty="n")
title(sprintf("%s — time-dependent ROC", GX)); dev.off()

## ---- section 16.3 required supplementary figure: distribution of Gene X across the four molecular subtypes ----
cat("\n--- distribution of Gene X across molecular subtypes (Kruskal-Wallis) ---\n")
kw <- kruskal.test(d2[[paste0(GX,"_z")]] ~ d2$subtype)
cat(sprintf("  Kruskal-Wallis P = %.3g\n", kw$p.value))
pw <- pairwise.wilcox.test(d2[[paste0(GX,"_z")]], d2$subtype, p.adjust.method="BH")
print(signif(pw$p.value,3))
pdf("results/survival/figures/NUMA1_by_subtype.pdf", width=5.5, height=5)
boxplot(d2[[paste0(GX,"_z")]] ~ d2$subtype, col="grey85", outline=FALSE,
        ylab=sprintf("%s (z)", GX), xlab="TCGA molecular subtype",
        main=sprintf("Kruskal-Wallis P = %.3g", kw$p.value))
points(jitter(as.numeric(d2$subtype),.6), d2[[paste0(GX,"_z")]], pch=16, cex=.5, col=adjustcolor("black",.4)); dev.off()

saveRDS(list(d=d, d2=d2, uni=u, m1=m1, m2=m2, timeROC=tr), "results/survival/survival_objects.rds")
writeLines(capture.output(sessionInfo()), "logs/16_prognosis_sessionInfo.txt")
cat("\n===== 16_TCGA_prognosis DONE =====\n")
