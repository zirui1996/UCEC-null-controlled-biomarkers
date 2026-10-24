#!/usr/bin/env Rscript
# 40m_revision_mc_cis.R -- exact binomial Monte Carlo 95% CIs for the empirical P values of Tables S14-S16 (R1 minor 4).
# b is recovered from P = (b+1)/(n+1) with the number of draws n used by each script (40e: 500; 40c: 1,000; 40f: 1,000).
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
mc <- function(P, n) { b <- round(P * (n + 1) - 1); ci <- binom.test(b, n)$conf.int
  data.frame(n_draws = n, b = b, P = signif((b + 1) / (n + 1), 3), MC_CI_low = signif(ci[1], 3), MC_CI_high = signif(ci[2], 3)) }
mv <- read.delim("results/revision/multivariable_gene_sets.tsv", check.names = FALSE)
s14 <- do.call(rbind, lapply(seq_len(nrow(mv)), function(i) rbind(
  cbind(table = "S14", item = paste(mv$gene_set[i], mv$cohort[i], "recurrent-consensus null"), mc(mv$selmatched_P[i], 500)),
  cbind(table = "S14", item = paste(mv$gene_set[i], mv$cohort[i], "non-consensus null"), mc(mv$noncons_P[i], 500)))))
pm <- read.delim("results/revision/GSE106191_admission_permutation_schemes.tsv", check.names = FALSE)
s15 <- do.call(rbind, lapply(seq_len(nrow(pm)), function(i) cbind(table = "S15", item = pm$scheme[i], mc(pm$empirical_P[i], 1000))))
fc <- read.delim("results/revision/four_cohort_design.tsv", check.names = FALSE)
s16 <- cbind(table = "S16", item = fc$design[2], mc(fc$empirical_P[2], 1000))
out <- rbind(s14, s15, s16); print(out, row.names = FALSE)
inP <- c(as.vector(rbind(mv$selmatched_P, mv$noncons_P)), pm$empirical_P, fc$empirical_P[2])
stopifnot(all(abs(out$P - signif(inP, 3)) < 1e-9))   # b recovered exactly
write.table(out, "results/revision/monte_carlo_CIs_S14_S16.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
