#!/usr/bin/env Rscript
# 40c_revision_GSE106191_permutation.R -- R2.3: how the collection-period covariate was handled in the admission permutation
# Original procedure (R/04b, analysis S-b, which is the design used for consensus): limma ~ group + period on all 97 samples;
# group labels were permuted across all samples while the period covariate stayed fixed (unrestricted permutation).
# Under H0, labels are exchangeable only within period if period is associated with diagnosis or with expression,
# so the correct scheme permutes labels within period (exchangeability blocks; Winkler et al. 2014).
# This script reruns the admission check both ways, plus a Freedman-Lane permutation of residuals, for comparison.
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages(library(limma))
set.seed(42); NPERM <- 1000

man <- read.delim("metadata/sample_manifest_GSE106191.tsv", stringsAsFactors = FALSE, quote = "")
x   <- as.matrix(read.delim("data_processed/expression/GSE106191_all97_gene_matrix_bestprobe.tsv", check.names = FALSE, row.names = 1))
x   <- x[, man$sample_id]
grp <- man$normal_or_tumor; era <- man$era
cat("period x group:\n"); print(table(era, grp))
cat(sprintf("Fisher test of period vs diagnosis: P = %.3f\n", fisher.test(table(era, grp))$p.value))

count_p001 <- function(g) {
  f <- factor(g, levels = c("normal","tumor")); e <- factor(make.names(era))
  dsg <- model.matrix(~ 0 + f + e); colnames(dsg)[1:2] <- c("N","T"); colnames(dsg) <- make.names(colnames(dsg))
  tt <- topTable(eBayes(contrasts.fit(lmFit(x, dsg), makeContrasts(T - N, levels = dsg))), number = Inf)
  sum(tt$P.Value < 0.001)
}
obs <- count_p001(grp)

perm_free   <- replicate(NPERM, count_p001(sample(grp)))
within      <- function() { g <- grp; for (b in unique(era)) { ix <- which(era == b); g[ix] <- g[ix][sample(length(ix))] }; g }
perm_within <- replicate(NPERM, count_p001(within()))

# Freedman-Lane: permute residuals of the reduced model (~ period) and add back fitted values
e <- factor(make.names(era)); red <- model.matrix(~ e)
fitred <- lmFit(x, red); fitted_red <- fitted(fitred); res_red <- x - fitted_red
f_obs <- factor(grp, levels = c("normal","tumor"))
full_count <- function(y) {
  dsg <- model.matrix(~ 0 + f_obs + e); colnames(dsg)[1:2] <- c("N","T"); colnames(dsg) <- make.names(colnames(dsg))
  tt <- topTable(eBayes(contrasts.fit(lmFit(y, dsg), makeContrasts(T - N, levels = dsg))), number = Inf)
  sum(tt$P.Value < 0.001)
}
perm_fl <- replicate(NPERM, { pi <- sample(ncol(x)); full_count(fitted_red + res_red[, pi]) })

summ <- function(p, lab) data.frame(scheme = lab, observed = obs, null_median = median(p), null_q95 = quantile(p, .95),
                                    null_max = max(p), empirical_P = (sum(p >= obs) + 1) / (NPERM + 1))
out <- rbind(summ(perm_free, "unrestricted label permutation, period fixed (original)"),
             summ(perm_within, "label permutation within collection period (exchangeability blocks)"),
             summ(perm_fl, "Freedman-Lane residual permutation (reduced model ~ period)"))
print(out, row.names = FALSE)
write.table(out, "results/revision/GSE106191_admission_permutation_schemes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("===== 40c DONE =====\n")
