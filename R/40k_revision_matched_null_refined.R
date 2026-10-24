#!/usr/bin/env Rscript
# 40k_revision_matched_null_refined.R -- sensitivity analysis for the effect-size-matched null of R/40a.
# In R/40a each core gene was replaced by a gene from the same 5% bin of |mean discovery log2FC|; all four core genes fall in the
# top bin (|log2FC| 0.87-3.39) and a panel could contain the same gene twice. Here each core gene is replaced by one of its k nearest
# neighbours on the signed mean discovery log2FC (same direction, closest magnitude), the four draws are distinct genes, and k is varied.
# Two candidate pools: the discovery universe (consensus genes included) and the recurrent consensus (genes significant in >=2 of
# the 3 discovery cohorts with concordant direction). No alternative four-gene set passes the strict 3-of-3 criterion, which
# selected exactly the four core genes, so these nulls approximate the selection process; they do not reproduce it.
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR)})
set.seed(4011)
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1"); NP <- 10000; KS <- c(10, 25, 50, 100)
fast_auc <- function(x, isT) { r <- rank(x); nT <- sum(isT); nN <- sum(!isT); (sum(r[isT]) - nT * (nT + 1) / 2) / (nT * nN) }
emp <- function(null, obs) { b <- sum(null >= obs); n <- length(null); ci <- binom.test(b, n)$conf.int
  c(b = b, P = (b + 1) / (n + 1), lo = ci[1], hi = ci[2]) }

V <- c("GSE17025","GSE63678","GSE106191")
L <- lapply(V, function(x) { t <- read.delim(sprintf("results/DEG/%s_all_genes.tsv", x), stringsAsFactors = FALSE)
  t <- t[!duplicated(t$gene), ]; setNames(t$logFC, t$gene) })
uni  <- Reduce(intersect, lapply(L, names))
mlfc <- rowMeans(sapply(L, function(x) x[uni])); names(mlfc) <- uni
rec  <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors = FALSE)$gene
cat("core mean discovery log2FC:", paste(sprintf("%s %.2f", CORE, mlfc[CORE]), collapse = "; "), "\n")
signed_score <- function(z, gg) { up <- gg[mlfc[gg] > 0]; dn <- gg[mlfc[gg] < 0]
  (if (length(up)) colMeans(z[up, , drop = FALSE]) else 0) - (if (length(dn)) colMeans(z[dn, , drop = FALSE]) else 0) }

# k nearest neighbours of each core gene on signed mean log2FC within a pool; one distinct draw per core gene
neighbours <- function(pool, k) setNames(lapply(CORE, function(g) { p <- setdiff(pool, CORE)
  p[order(abs(mlfc[p] - mlfc[g]))][seq_len(min(k, length(p)))] }), CORE)
draw_nn <- function(nb) { out <- character(0)
  for (g in sample(CORE)) out <- c(out, { c0 <- setdiff(nb[[g]], out); c0[sample.int(length(c0), 1)] }); out }

run <- function(z, isT, cohort) {
  obs <- fast_auc(-colMeans(z[CORE, ]), isT)
  pools <- list(`discovery universe` = intersect(uni, rownames(z)), `recurrent consensus` = intersect(rec, rownames(z)))
  do.call(rbind, lapply(names(pools), function(pn) do.call(rbind, lapply(KS, function(k) {
    if (pn == "recurrent consensus" && k > 50) return(NULL)
    nb <- neighbours(pools[[pn]], k)
    rng <- range(unlist(lapply(nb, function(x) mlfc[x])))
    nl <- replicate(NP, { a <- fast_auc(signed_score(z, draw_nn(nb)), isT); max(a, 1 - a) })
    e <- emp(nl, obs)
    data.frame(cohort = cohort, pool = pn, k_nearest = k, matched_log2FC_range = sprintf("%.2f to %.2f", rng[1], rng[2]),
               draws = NP, observed_AUC = round(obs, 4), null_median = round(median(nl), 4), null_q975 = round(quantile(nl, .975), 4),
               n_ge_obs = e[["b"]], empirical_P = signif(e[["P"]], 3), MC_CI_low = signif(e[["lo"]], 3), MC_CI_high = signif(e[["hi"]], 3))
  }))))
}

## TCGA-UCEC (as R/40a)
se  <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop = FALSE]); code <- substr(cdl$barcode, 14, 15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code == "01" & cdl$primary_diagnosis %in% c(EN, SR)) | code == "11"
cnt <- assay(se, "unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]; cnt <- cnt[i, ]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt) != "", ]
grp <- factor(ifelse(code[keep] == "11", "Normal", "Tumor"), levels = c("Normal","Tumor"))
d <- DGEList(cnt); d <- d[filterByExpr(d, group = grp), , keep.lib.sizes = FALSE]; d <- calcNormFactors(d)
zT <- t(scale(t(cpm(d, log = TRUE))))
rT <- run(zT, grp == "Tumor", "TCGA-UCEC")

## GSE146889 (as R/40a: all 72 annotated endometrial samples)
man <- subset(read.delim("metadata/sample_manifest_GSE146889.tsv", stringsAsFactors = FALSE), tissue == "Endometrial")
raw <- read.delim("data_raw/GEO/GSE146889/GSE146889_GeneCount.tsv.gz", check.names = FALSE)
gc  <- raw[, paste0(man$column, "_count")]; colnames(gc) <- man$column
sym <- raw$GeneName; ord <- order(sym, -rowMeans(gc)); idx <- ord[!duplicated(sym[ord])]
mm  <- as.matrix(gc[idx, ]); rownames(mm) <- sym[idx]; mm <- mm[!is.na(rownames(mm)) & rownames(mm) != "", ]
gv  <- factor(man$group, levels = c("normal","tumor"))
dv  <- DGEList(mm); dv <- dv[filterByExpr(dv, group = gv), , keep.lib.sizes = FALSE]; dv <- calcNormFactors(dv)
zV  <- t(scale(t(cpm(dv, log = TRUE))))
rV <- run(zV, gv == "tumor", "GSE146889")

R <- rbind(rT, rV); print(R, row.names = FALSE)
write.table(R, "results/revision/matched_null_nearest_neighbour.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("===== 40k DONE =====\n")
