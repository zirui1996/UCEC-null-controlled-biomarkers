#!/usr/bin/env Rscript
# 40a_revision_discrimination.R -- reanalyses requested in peer review (FIG, major revision)
#   R1.1 / R2.2  four-gene panel AUC, CI and empirical P in the held-out GSE146889 cohort
#   R2.1 / R3.3  alternative null definitions: non-consensus (original), broad (all discovery-universe genes),
#                selection-matched (random four-gene subsets of the recurrent consensus) and effect-size-matched
#   R1.2         patient-clustered bootstrap CI (TCGA has 21 tumour-normal pairs), fully unpaired subset,
#                within-pair concordance, and restriction to tissue source sites that contributed normals
#   R1 minor 4   Monte Carlo 95% CIs for every empirical P
#   R1 minor 5   normalization robustness: TMM vs upper-quartile vs RLE (DESeq-type median-of-ratios) vs none
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(limma); library(pROC)})
set.seed(42)
dir.create("results/revision", showWarnings = FALSE, recursive = TRUE)
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1"); NP <- 10000; NB <- 2000

## ---------------------------------------------------------------- shared helpers
fast_auc <- function(x, isT) { r <- rank(x); nT <- sum(isT); nN <- sum(!isT)
  (sum(r[isT]) - nT * (nT + 1) / 2) / (nT * nN) }
mc_ci <- function(b, n) { ci <- binom.test(b, n)$conf.int; c(lo = ci[1], hi = ci[2]) }
emp <- function(null, obs) { b <- sum(null >= obs); n <- length(null)
  c(b = b, n = n, P = (b + 1) / (n + 1), mc_ci(b, n)) }

V <- c("GSE17025","GSE63678","GSE106191")
L <- lapply(V, function(x) { t <- read.delim(sprintf("results/DEG/%s_all_genes.tsv", x), stringsAsFactors = FALSE)
  t <- t[!duplicated(t$gene), ]; setNames(t$logFC, t$gene) })
uni  <- Reduce(intersect, lapply(L, names))
mlfc <- rowMeans(sapply(L, function(x) x[uni])); names(mlfc) <- uni
rec  <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors = FALSE)$gene

# signed score: mean z of up genes minus mean z of down genes, sign from mean discovery log2FC
signed_score <- function(z, gg) { gg <- intersect(gg, rownames(z)); up <- gg[mlfc[gg] > 0]; dn <- gg[mlfc[gg] < 0]
  (if (length(up)) colMeans(z[up, , drop = FALSE]) else 0) - (if (length(dn)) colMeans(z[dn, , drop = FALSE]) else 0) }

# the four null pools used in the revision
make_pools <- function(measurable) {
  m <- function(g) intersect(g, measurable)
  br <- quantile(abs(mlfc), probs = seq(0, 1, 0.05)); br[1] <- -Inf; br[length(br)] <- Inf
  bins <- cut(abs(mlfc), br); names(bins) <- names(mlfc)
  list(noncons  = m(setdiff(uni, rec)),                       # original analysis
       broad    = m(setdiff(uni, CORE)),                       # Venet-style: whole discovery universe, consensus included
       selmatch = m(setdiff(rec, CORE)),                       # panels built from genes passing the same consensus filter
       bins     = bins, pool_es = m(setdiff(uni, CORE)))
}
draw <- function(pools, type) {
  if (type == "esmatch") {                                     # effect-size matched: same |mean log2FC| bin as each core gene
    return(vapply(CORE, function(g) { cand <- pools$pool_es[pools$bins[pools$pool_es] == pools$bins[g]]
      sample(cand, 1) }, ""))
  }
  sample(pools[[type]], 4)
}
null_table <- function(z, isT, obs, pools, cohort) {
  types <- c(noncons = "non-consensus (original)", broad = "broad discovery universe",
             selmatch = "selection-matched (recurrent consensus)", esmatch = "effect-size-matched")
  do.call(rbind, lapply(names(types), function(ty) {
    nl <- replicate(NP, { a <- fast_auc(signed_score(z, draw(pools, ty)), isT); max(a, 1 - a) })
    e <- emp(nl, obs)
    data.frame(cohort = cohort, null = types[[ty]], pool_size = if (ty == "esmatch") NA else length(pools[[ty]]),
               draws = NP, observed_AUC = round(obs, 4), null_median = round(median(nl), 4),
               null_q025 = round(quantile(nl, .025), 4), null_q975 = round(quantile(nl, .975), 4),
               n_ge_obs = e[["b"]], empirical_P = signif(e[["P"]], 3),
               MC_CI_low = signif(e[["lo"]], 3), MC_CI_high = signif(e[["hi"]], 3))
  }))
}
cluster_boot <- function(score, isT, pid, B = NB) {               # resample patients, keep all of a patient's samples
  ids <- unique(pid); idx <- split(seq_along(pid), pid)
  out <- replicate(B, { s <- unlist(idx[sample(ids, length(ids), replace = TRUE)], use.names = FALSE)
    if (length(unique(isT[s])) < 2) NA else fast_auc(score[s], isT[s]) })
  quantile(out, c(.025, .975), na.rm = TRUE)
}

## ---------------------------------------------------------------- TCGA-UCEC
se  <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop = FALSE])
code <- substr(cdl$barcode, 14, 15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code == "01" & cdl$primary_diagnosis %in% c(EN, SR)) | code == "11"
cnt <- assay(se, "unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]
cnt <- cnt[i, ]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt) != "", ]
bc  <- cdl$barcode[keep]
grp <- factor(ifelse(code[keep] == "11", "Normal", "Tumor"), levels = c("Normal","Tumor")); isT <- grp == "Tumor"
pid <- substr(bc, 1, 12); tss <- substr(bc, 6, 7); plate <- substr(bc, 22, 25)
d0  <- DGEList(cnt); d0 <- d0[filterByExpr(d0, group = grp), , keep.lib.sizes = FALSE]
d   <- calcNormFactors(d0, method = "TMM")
lc  <- cpm(d, log = TRUE); z <- t(scale(t(lc)))
obsT <- fast_auc(-colMeans(z[CORE, ]), isT)
stopifnot(abs(obsT - 0.98839) < 1e-4)                         # reproduces the published 0.988
cat(sprintf("TCGA: %d tumours / %d normals | panel AUC %.4f\n", sum(isT), sum(!isT), obsT))

poolsT <- make_pools(rownames(z))
nullT  <- null_table(z, isT, obsT, poolsT, "TCGA-UCEC")
print(nullT[, c("null","pool_size","null_median","null_q975","n_ge_obs","empirical_P","MC_CI_low","MC_CI_high")], row.names = FALSE)

## ---- R1.2: dependence and sampling artefacts
sT <- -colMeans(z[CORE, ])
rocT <- roc(grp, sT, quiet = TRUE, direction = "<"); dl <- as.numeric(ci.auc(rocT))
cb  <- cluster_boot(sT, isT, pid)
paired_pids <- intersect(pid[isT], pid[!isT])
unp <- !(isT & pid %in% paired_pids)                           # drop tumours of patients who also gave a normal
a_unp <- fast_auc(sT[unp], isT[unp]); ci_unp <- as.numeric(ci.auc(roc(grp[unp], sT[unp], quiet = TRUE, direction = "<")))
wp <- sapply(paired_pids, function(p) mean(sT[pid == p & isT]) - mean(sT[pid == p & !isT]))
tss_n <- unique(tss[!isT]); sub_tss <- tss %in% tss_n
a_tss <- fast_auc(sT[sub_tss], isT[sub_tss])
ci_tss <- as.numeric(ci.auc(roc(grp[sub_tss], sT[sub_tss], quiet = TRUE, direction = "<")))
plate_n <- unique(plate[!isT]); sub_pl <- plate %in% plate_n
a_pl <- fast_auc(sT[sub_pl], isT[sub_pl])
dep <- data.frame(
  analysis = c("all samples, DeLong", "all samples, patient-clustered bootstrap",
               "fully unpaired (tumours of the 21 paired patients removed), DeLong",
               "21 matched pairs: tumour score above normal score",
               "tissue source sites that contributed normals only, DeLong",
               "sequencing plates that contained normals only"),
  n_tumour = c(sum(isT), sum(isT), sum(unp & isT), length(paired_pids), sum(sub_tss & isT), sum(sub_pl & isT)),
  n_normal = c(sum(!isT), sum(!isT), sum(unp & !isT), length(paired_pids), sum(sub_tss & !isT), sum(sub_pl & !isT)),
  estimate = c(obsT, obsT, a_unp, mean(wp > 0), a_tss, a_pl),
  CI_low = c(dl[1], cb[1], ci_unp[1], binom.test(sum(wp > 0), length(wp))$conf.int[1], ci_tss[1], NA),
  CI_high = c(dl[3], cb[2], ci_unp[3], binom.test(sum(wp > 0), length(wp))$conf.int[2], ci_tss[3], NA))
dep$note <- c("", sprintf("%d bootstrap resamples of %d patients", NB, length(unique(pid))), "",
              sprintf("%d of %d pairs; exact binomial CI; signed-rank P=%.2g", sum(wp > 0), length(wp),
                      wilcox.test(wp)$p.value),
              sprintf("%d sites contributed normals: %s", length(tss_n), paste(sort(tss_n), collapse = ",")),
              sprintf("%d plates contained normals", length(plate_n)))
print(dep, row.names = FALSE, digits = 3)

## ---- R1 minor 5: normalization robustness
norm_res <- do.call(rbind, lapply(c("TMM","upperquartile","RLE","none"), function(meth) {
  dn <- calcNormFactors(d0, method = meth); l <- cpm(dn, log = TRUE); zz <- t(scale(t(l)))
  v  <- voom(dn, model.matrix(~ grp)); fit <- eBayes(lmFit(v, model.matrix(~ grp)))
  tt <- topTable(fit, coef = 2, number = Inf)
  nf <- dn$samples$norm.factors
  data.frame(method = meth,
             norm_factor_median_tumour = round(median(nf[isT]), 3), norm_factor_median_normal = round(median(nf[!isT]), 3),
             pct_DE_FDR05 = round(100 * mean(tt$adj.P.Val < 0.05), 1),
             panel_AUC = round(fast_auc(-colMeans(zz[CORE, ]), isT), 4),
             NUMA1_AUC = round({ a <- fast_auc(l["NUMA1", ], isT); max(a, 1 - a) }, 4),
             median_logFC_nonDE = round(median(tt$logFC[tt$adj.P.Val >= 0.05]), 3))
}))
print(norm_res, row.names = FALSE)
# MA plot under TMM: per-gene mean log-CPM vs tumour-minus-normal difference, with loess trend
pdf("results/revision/TCGA_MA_TMM.pdf", width = 5.5, height = 4.5)
A <- rowMeans(lc); M <- rowMeans(lc[, isT]) - rowMeans(lc[, !isT])
smoothScatter(A, M, xlab = "Mean log2 CPM (TMM)", ylab = "Tumour minus normal (log2)", main = "TCGA-UCEC, TMM")
abline(h = 0, lty = 2); lines(lowess(A, M, f = 0.3), col = "red", lwd = 2)
dev.off()

## ---------------------------------------------------------------- GSE146889 (held out)
man <- subset(read.delim("metadata/sample_manifest_GSE146889.tsv", stringsAsFactors = FALSE), tissue == "Endometrial")
raw <- read.delim("data_raw/GEO/GSE146889/GSE146889_GeneCount.tsv.gz", check.names = FALSE)
gc  <- raw[, paste0(man$column, "_count")]; colnames(gc) <- man$column
sym <- raw$GeneName; ord <- order(sym, -rowMeans(gc)); idx <- ord[!duplicated(sym[ord])]
mm  <- as.matrix(gc[idx, ]); rownames(mm) <- sym[idx]; mm <- mm[!is.na(rownames(mm)) & rownames(mm) != "", ]
gv  <- factor(man$group, levels = c("normal","tumor")); isTv <- gv == "tumor"
dv  <- DGEList(mm); dv <- dv[filterByExpr(dv, group = gv), , keep.lib.sizes = FALSE]; dv <- calcNormFactors(dv)
zv  <- t(scale(t(cpm(dv, log = TRUE))))
stopifnot(all(CORE %in% rownames(zv)))
sV  <- -colMeans(zv[CORE, ]); obsV <- fast_auc(sV, isTv)
rocV <- roc(gv, sV, quiet = TRUE, direction = "<"); dlV <- as.numeric(ci.auc(rocV))
pidV <- if ("patient" %in% names(man)) man$patient else man$pair_id
cbV <- cluster_boot(sV, isTv, pidV)
pp  <- intersect(pidV[isTv], pidV[!isTv])
wpV <- sapply(pp, function(p) mean(sV[pidV == p & isTv]) - mean(sV[pidV == p & !isTv]))
cat(sprintf("\nGSE146889: %d tumours / %d normals, %d pairs | panel AUC %.4f (DeLong %.4f-%.4f; clustered %.4f-%.4f) | pairs with lower tumour score %d/%d\n",
            sum(isTv), sum(!isTv), length(pp), obsV, dlV[1], dlV[3], cbV[1], cbV[2], sum(wpV < 0), length(wpV)))
poolsV <- make_pools(rownames(zv))
nullV  <- null_table(zv, isTv, obsV, poolsV, "GSE146889")
# a second, paired statistic: mean within-pair score difference (tumour - normal), oriented so larger = better
pair_stat <- function(s) { w <- sapply(pp, function(p) mean(s[pidV == p & isTv]) - mean(s[pidV == p & !isTv])); abs(mean(w)) / sd(w) }
obs_ps <- pair_stat(sV)
nl_ps  <- replicate(NP, pair_stat(signed_score(zv, draw(poolsV, "noncons"))))
e_ps   <- emp(nl_ps, obs_ps)
print(nullV[, c("null","pool_size","null_median","null_q975","n_ge_obs","empirical_P","MC_CI_low","MC_CI_high")], row.names = FALSE)
cat(sprintf("paired effect size (mean/SD of within-pair difference) %.2f | non-consensus null median %.2f | empirical P %.3g (MC 95%% CI %.3g-%.3g)\n",
            obs_ps, median(nl_ps), e_ps[["P"]], e_ps[["lo"]], e_ps[["hi"]]))

## ---------------------------------------------------------------- Monte Carlo CIs for the published empirical P values
pub <- read.delim("results/validation/TCGA_null_controls.tsv"); sig <- read.delim("results/validation/TCGA_signature_null.tsv")
g146 <- read.delim("results/validation/GSE146889_null_controls.tsv")
mcp <- data.frame(
  quantity = c("NUMA1 single-gene AUC (TCGA)", "core genes single-gene AUC (TCGA)", "four-gene panel (TCGA)",
               "recurrent signature AUC (TCGA)", "direction concordance (GSE146889)", "Spearman rho (GSE146889)",
               "recurrent signature AUC (GSE146889)"),
  null = c("exhaustive over 11,800 genes", "exhaustive over 11,800 genes", "10,000 random panels",
           "1,000 random signatures", "1,000 effect-size-matched sets", "1,000 effect-size-matched sets", "1,000 random signatures"),
  empirical_P = c(pub$empirical_P[pub$feature == "NUMA1"], NA, nullT$empirical_P[1], sig$empirical_P,
                  g146$empirical_P[1], g146$empirical_P[2], g146$empirical_P[3]))
mcp$MC_CI <- c("none (exhaustive; no Monte Carlo error)", "none (exhaustive; no Monte Carlo error)",
               sprintf("%.4f-%.4f", nullT$MC_CI_low[1], nullT$MC_CI_high[1]),
               { b <- round(sig$empirical_P * 1000); ci <- mc_ci(b, 1000); sprintf("%.3f-%.3f", ci[1], ci[2]) },
               { b <- round(g146$empirical_P[1] * 1000); ci <- mc_ci(b, 1000); sprintf("%.4f-%.4f", ci[1], ci[2]) },
               { b <- round(g146$empirical_P[2] * 1000); ci <- mc_ci(b, 1000); sprintf("%.4f-%.4f", ci[1], ci[2]) },
               { b <- round(g146$empirical_P[3] * 1000); ci <- mc_ci(b, 1000); sprintf("%.3f-%.3f", ci[1], ci[2]) })
print(mcp, row.names = FALSE)

## ---------------------------------------------------------------- write
write.table(rbind(nullT, nullV), "results/revision/panel_null_definitions.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dep, "results/revision/TCGA_dependence_artefact_checks.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(norm_res, "results/revision/TCGA_normalization_sensitivity.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(mcp, "results/revision/monte_carlo_CIs.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(cohort = "GSE146889", n_tumour = sum(isTv), n_normal = sum(!isTv), n_pairs = length(pp),
                       AUC = obsV, DeLong_low = dlV[1], DeLong_high = dlV[3], cluster_low = cbV[1], cluster_high = cbV[2],
                       pairs_lower_tumour = sum(wpV < 0), paired_effect = obs_ps, paired_null_median = median(nl_ps),
                       paired_empirical_P = e_ps[["P"]], paired_MC_low = e_ps[["lo"]], paired_MC_high = e_ps[["hi"]]),
            "results/revision/GSE146889_panel.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(capture.output(sessionInfo()), "results/revision/40a_sessionInfo.txt")
cat("===== 40a DONE =====\n")
