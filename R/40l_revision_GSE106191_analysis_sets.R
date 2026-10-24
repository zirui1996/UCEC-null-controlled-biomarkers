#!/usr/bin/env Rscript
# 40l_revision_GSE106191_analysis_sets.R -- GSE106191 under its three analysis sets (protocol section 0.6.3, R/04b):
#   prespecified main  2016 collection period only (54 samples), limma ~ scan day + group
#   S-a                the three 2016 scan days that contained both diagnoses (34 samples), ~ scan day + group
#   adopted (S-b)      all 97 samples, ~ collection period + group  (used for consensus after the first two failed admission)
# For each set: admission statistics (genes at raw P<0.001 vs 1,000 label permutations within the blocking factor, empirical
# P=(b+1)/(n+1) with exact binomial Monte Carlo CI), genes at FDR<0.05, the four core genes, and the consensus that would result.
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages(library(limma))
set.seed(4013); NPERM <- 1000; FDR <- 0.05; LFC <- 1
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")
man <- read.delim("metadata/sample_manifest_GSE106191.tsv", stringsAsFactors = FALSE, quote = "")
x16 <- as.matrix(read.delim("data_processed/expression/GSE106191_gene_matrix_bestprobe.tsv", check.names = FALSE, row.names = 1))
x97 <- as.matrix(read.delim("data_processed/expression/GSE106191_all97_gene_matrix_bestprobe.tsv", check.names = FALSE, row.names = 1))

fit_tt <- function(x, grp, blk) {
  f <- factor(grp, levels = c("normal","tumor")); b <- factor(make.names(blk))
  d <- model.matrix(~ 0 + f + b); colnames(d) <- make.names(colnames(d)); colnames(d)[1:2] <- c("N","T")
  topTable(eBayes(contrasts.fit(lmFit(x, d), makeContrasts(T - N, levels = d))), number = Inf)
}
perm_within <- function(grp, blk) { g <- grp; for (k in unique(blk)) { ix <- which(blk == k); g[ix] <- g[ix][sample.int(length(ix))] }; g }

sets <- list(
  `prespecified main: 2016 period, ~ scan day` = list(x = x16, s = subset(man, include_primary == "yes"), blk = "scan_day"),
  `S-a: 2016 mixed scan days, ~ scan day`      = list(x = x16, s = subset(man, mixed_day == "yes"), blk = "scan_day"),
  `adopted: all 97, ~ collection period`       = list(x = x97, s = man, blk = "era"))

# consensus that each analysis set would have produced (vote counting as R/05_consensus.R)
D <- lapply(c("GSE17025","GSE63678"), function(g) { t <- read.delim(sprintf("results/DEG/%s_all_genes.tsv", g)); t[!duplicated(t$gene), ] })
names(D) <- c("GSE17025","GSE63678")

res <- list(); core <- list()
for (nm in names(sets)) {
  S <- sets[[nm]]; x <- S$x[, S$s$sample_id]; grp <- S$s$normal_or_tumor; blk <- S$s[[S$blk]]
  tt <- fit_tt(x, grp, blk); obs <- sum(tt$P.Value < 0.001)
  pm <- replicate(NPERM, sum(fit_tt(x, perm_within(grp, blk), blk)$P.Value < 0.001))
  b <- sum(pm >= obs); ci <- binom.test(b, NPERM)$conf.int
  # informative samples: blocks containing both diagnoses
  inf <- tapply(grp, blk, function(g) length(unique(g)) == 2); inf_n <- sum(blk %in% names(inf)[inf])
  # consensus with this GSE106191 result
  uni <- Reduce(intersect, list(D$GSE17025$gene, D$GSE63678$gene, rownames(tt)))
  L <- list(D$GSE17025[match(uni, D$GSE17025$gene), ], D$GSE63678[match(uni, D$GSE63678$gene), ],
            data.frame(gene = uni, logFC = tt[uni, "logFC"], adj.P.Val = tt[uni, "adj.P.Val"]))
  sg <- sapply(L, function(t) t$adj.P.Val < FDR & abs(t$logFC) >= LFC); lf <- sapply(L, function(t) t$logFC); rownames(sg) <- uni
  ok <- apply(sg & lf > 0, 1, sum) == rowSums(sg) | apply(sg & lf < 0, 1, sum) == rowSums(sg)
  res[[nm]] <- data.frame(analysis_set = nm, n_tumour = sum(grp == "tumor"), n_hyperplasia = sum(grp == "normal"),
                          samples_in_blocks_with_both_groups = inf_n,
                          genes_FDR005 = sum(tt$adj.P.Val < FDR), min_adjP = signif(min(tt$adj.P.Val), 3),
                          observed_rawP001 = obs, perm_median = median(pm), n_perm_ge_obs = b,
                          empirical_P = signif((b + 1) / (NPERM + 1), 3), MC_CI_low = signif(ci[1], 3), MC_CI_high = signif(ci[2], 3),
                          strict_3of3_genes = sum(rowSums(sg) == 3 & ok), recurrent_2of3_genes = sum(rowSums(sg) >= 2 & ok),
                          core4_in_recurrent = sum(CORE %in% uni[rowSums(sg) >= 2 & ok]))
  core[[nm]] <- data.frame(analysis_set = nm, gene = CORE, log2FC = round(tt[CORE, "logFC"], 3), adjP = signif(tt[CORE, "adj.P.Val"], 3))
}
R <- do.call(rbind, res); C <- do.call(rbind, core)
print(R, row.names = FALSE); print(C, row.names = FALSE)
write.table(R, "results/revision/GSE106191_analysis_sets.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(C, "results/revision/GSE106191_analysis_sets_core_genes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("===== 40l DONE =====\n")
