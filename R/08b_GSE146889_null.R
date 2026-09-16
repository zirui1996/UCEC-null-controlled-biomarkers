#!/usr/bin/env Rscript
# 08b_GSE146889_null.R -- null distributions for the GSE146889 validation metrics (written to disk in v2.6; in v2.4 this computation was not saved as a script)
# random gene pool = genes shared by the three discovery cohorts that are outside the recurrent consensus and measurable in GSE146889
# direction taken from the sign of the mean discovery logFC; 1000 iterations
#   (1) direction concordance: random vs effect-size matched (binned into 20 quantiles of |mean discovery logFC|)
#   (2) Spearman rho between discovery and replication logFC: effect-size matched
#   (3) sample-level signature AUC: random signatures of the same size
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(edgeR); library(pROC)})
set.seed(42); NB <- 1000
V <- c("GSE17025","GSE63678","GSE106191")
L <- lapply(V, function(g){ t <- read.delim(sprintf("results/DEG/%s_all_genes.tsv", g), stringsAsFactors=FALSE)
  t <- t[!duplicated(t$gene),]; setNames(t$logFC, t$gene) })
uni <- Reduce(intersect, lapply(L, names))
mlfc <- rowMeans(sapply(L, function(x) x[uni])); names(mlfc) <- uni
rec <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors=FALSE)$gene
v <- read.delim("results/DEG/GSE146889_all_genes.tsv", stringsAsFactors=FALSE); v <- v[!duplicated(v$gene),]; rownames(v) <- v$gene
ok <- intersect(rec, rownames(v)); pool <- intersect(setdiff(uni, rec), rownames(v))
cat(sprintf("recurrent %d | measurable in GSE146889 %d | random pool %d\n", length(rec), length(ok), length(pool)))

conc <- function(g) mean(sign(mlfc[g]) == sign(v[g, "logFC"]))
rho  <- function(g) suppressWarnings(cor(mlfc[g], v[g, "logFC"], method = "spearman"))
obs_c <- conc(ok); obs_r <- rho(ok)

br <- quantile(abs(mlfc[c(ok, pool)]), probs = seq(0, 1, 0.05)); br[1] <- -Inf; br[length(br)] <- Inf
bin_ok <- cut(abs(mlfc[ok]), br); bin_pool <- cut(abs(mlfc[pool]), br)
split_pool <- split(pool, bin_pool)
matched <- function() unlist(lapply(names(table(bin_ok)), function(b) {
  k <- sum(bin_ok == b); if (!k) return(NULL); sample(split_pool[[b]], k, replace = length(split_pool[[b]]) < k) }))

n_rand_c <- replicate(NB, conc(sample(pool, length(ok))))
n_mat <- replicate(NB, { g <- matched(); c(conc(g), rho(g)) })
cat(sprintf("\n(1) direction concordance %d/%d = %.1f%% | random median %.1f%% | effect-size-matched median %.1f%% (95%% %.1f-%.1f%%) | empirical P = %s\n",
    sum(sign(mlfc[ok]) == sign(v[ok,"logFC"])), length(ok), 100*obs_c, 100*median(n_rand_c), 100*median(n_mat[1,]),
    100*quantile(n_mat[1,], .025), 100*quantile(n_mat[1,], .975),
    ifelse(mean(n_mat[1,] >= obs_c) == 0, sprintf("< %.3f", 1/NB), sprintf("%.3f", mean(n_mat[1,] >= obs_c)))))
cat(sprintf("(2) Spearman rho = %.3f | effect-size-matched median %.3f (95%% %.3f-%.3f) | empirical P = %s\n",
    obs_r, median(n_mat[2,]), quantile(n_mat[2,], .025), quantile(n_mat[2,], .975),
    ifelse(mean(n_mat[2,] >= obs_r) == 0, sprintf("< %.3f", 1/NB), sprintf("%.3f", mean(n_mat[2,] >= obs_r)))))
cat(sprintf("    replication at FDR<0.05 %.1f%% | genome-wide baseline %.1f%%\n", 100*mean(v[ok,"adj.P.Val"] < 0.05), 100*mean(v$adj.P.Val < 0.05)))

## (3) sample-level AUC
man <- subset(read.delim("metadata/sample_manifest_GSE146889.tsv", stringsAsFactors=FALSE), tissue=="Endometrial")
raw <- read.delim("data_raw/GEO/GSE146889/GSE146889_GeneCount.tsv.gz", check.names=FALSE)
cnt <- raw[, paste0(man$column,"_count")]; colnames(cnt) <- man$column
sym <- raw$GeneName; ord <- order(sym, -rowMeans(cnt)); idx <- ord[!duplicated(sym[ord])]
m <- as.matrix(cnt[idx,]); rownames(m) <- sym[idx]; m <- m[!is.na(rownames(m)) & rownames(m)!="",]
d <- DGEList(m); d <- d[filterByExpr(d, group=man$group),,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
z <- t(scale(t(cpm(d, log=TRUE)))); lab <- factor(man$group, levels=c("normal","tumor"))
sc <- function(g){ g <- intersect(g, rownames(z)); up <- g[mlfc[g] > 0]; dn <- g[mlfc[g] < 0]
  (if (length(up)) colMeans(z[up,,drop=FALSE]) else 0) - (if (length(dn)) colMeans(z[dn,,drop=FALSE]) else 0) }
okz <- intersect(ok, rownames(z)); poolz <- intersect(pool, rownames(z))
obs_a <- as.numeric(roc(lab, sc(okz), quiet=TRUE)$auc)
n_a <- replicate(NB, as.numeric(roc(lab, sc(sample(poolz, length(okz))), quiet=TRUE)$auc))
cat(sprintf("(3) signature AUC (%d genes) = %.3f | random median %.3f (95%% %.3f-%.3f) | empirical P = %.3f\n",
    length(okz), obs_a, median(n_a), quantile(n_a,.025), quantile(n_a,.975), mean(n_a >= obs_a)))

out <- data.frame(metric=c("direction_concordance","spearman_rho","signature_AUC"),
  observed=c(obs_c, obs_r, obs_a), null_median=c(median(n_mat[1,]), median(n_mat[2,]), median(n_a)),
  null_q025=c(quantile(n_mat[1,],.025), quantile(n_mat[2,],.025), quantile(n_a,.025)),
  null_q975=c(quantile(n_mat[1,],.975), quantile(n_mat[2,],.975), quantile(n_a,.975)),
  empirical_P=c(mean(n_mat[1,]>=obs_c), mean(n_mat[2,]>=obs_r), mean(n_a>=obs_a)),
  null_type=c("effect-size matched","effect-size matched","random same size"), n_genes=c(length(ok), length(ok), length(okz)))
write.table(out, "results/validation/GSE146889_null_controls.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("    (for reference, median direction concordance of purely random sets is %.1f%%)\n===== 08b DONE =====\n", 100*median(n_rand_c)))
