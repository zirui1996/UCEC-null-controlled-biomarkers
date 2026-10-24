#!/usr/bin/env Rscript
# 40b_revision_correlation.R -- R1.3: the selection-frequency vs single-gene AUC correlation (rho = -0.34, P = 3.9e-6)
# treated genes as independent. This script quantifies two sources of uncertainty the original P ignored:
#   (1) Monte Carlo error of the selection frequencies (each from 1,000 bootstrap iterations): binomial propagation
#   (2) dependence among genes (co-expression, shared cell-cycle regulation, cell-type markers):
#       co-expression module block bootstrap, and the effective number of independent genes (Li & Ji 2005)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR)})
set.seed(42); B <- 5000; NIT <- 1000

st  <- read.delim("results/validation/TCGA_recurrent_single_gene_AUC.tsv", stringsAsFactors = FALSE)
f   <- setNames(st$freq_intersection, st$gene); a <- setNames(st$AUC, st$gene)
g   <- st$gene; n <- length(g)
obs <- cor(f, a, method = "spearman"); p_naive <- cor.test(f, a, method = "spearman", exact = FALSE)$p.value
cat(sprintf("genes %d | Spearman rho %.3f | naive P %.2g\n", n, obs, p_naive))

## TCGA expression of the 181 genes, to measure their co-expression
se  <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop = FALSE])
code <- substr(cdl$barcode, 14, 15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code == "01" & cdl$primary_diagnosis %in% c(EN, SR)) | code == "11"
cnt <- assay(se, "unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]
cnt <- cnt[i, ]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt) != "", ]
grp <- factor(ifelse(code[keep] == "11", "Normal", "Tumor"))
d <- DGEList(cnt); d <- d[filterByExpr(d, group = grp), , keep.lib.sizes = FALSE]; d <- calcNormFactors(d)
lc <- cpm(d, log = TRUE)[g, ]
R  <- cor(t(lc), method = "spearman")

## (1) Monte Carlo error of selection frequency
rho_mc <- replicate(B, cor(rbinom(n, NIT, f) / NIT, a, method = "spearman"))

## (2a) co-expression modules: average-linkage clustering on 1 - |rho|, cut where within-module |rho| >= 0.5
hc <- hclust(as.dist(1 - abs(R)), method = "average")
mod <- cutree(hc, h = 0.5); k <- length(unique(mod)); ms <- split(g, mod)
rho_blk <- replicate(B, { s <- unlist(ms[sample(names(ms), k, replace = TRUE)], use.names = FALSE)
  if (length(unique(s)) < 5) NA else cor(f[s], a[s], method = "spearman") })
# both sources together: module bootstrap with binomially perturbed frequencies
rho_both <- replicate(B, { s <- unlist(ms[sample(names(ms), k, replace = TRUE)], use.names = FALSE)
  fb <- rbinom(length(s), NIT, f[s]) / NIT; cor(fb, a[s], method = "spearman") })

## (2b) effective number of independent genes (Li & Ji 2005) and the corresponding t-test
ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values; ev <- pmax(ev, 0)
meff <- sum((ev >= 1) + (ev - floor(ev)))
tstat <- obs * sqrt((meff - 2) / (1 - obs^2)); p_meff <- 2 * pt(-abs(tstat), df = meff - 2)

## (2c) the same correlation within the 101 up- and 81 down-regulated genes, and without cell-cycle genes
cc <- setNames(st$cell_cycle, st$gene); dirn <- setNames(st$direction, st$gene)
sub_rho <- function(s) { r <- cor.test(f[s], a[s], method = "spearman", exact = FALSE); c(n = length(s), rho = unname(r$estimate), P = r$p.value) }
subs <- rbind(up = sub_rho(g[dirn == "up"]), down = sub_rho(g[dirn == "down"]), non_cell_cycle = sub_rho(g[!cc]))

q <- function(x) sprintf("%.3f to %.3f", quantile(x, .025, na.rm = TRUE), quantile(x, .975, na.rm = TRUE))
out <- data.frame(
  analysis = c("observed (genes treated as independent)", "Monte Carlo error of selection frequency (binomial, 1,000 iterations per gene)",
               sprintf("co-expression module bootstrap (%d modules at |rho|>=0.5)", k),
               "module bootstrap + Monte Carlo error", sprintf("effective number of independent genes, Meff = %.1f (Li & Ji)", meff)),
  rho = round(c(obs, median(rho_mc), median(rho_blk, na.rm = TRUE), median(rho_both), obs), 3),
  CI95 = c("", q(rho_mc), q(rho_blk), q(rho_both), ""),
  P = signif(c(p_naive, NA, NA, NA, p_meff), 3),
  share_below_zero = c(NA, mean(rho_mc < 0), mean(rho_blk < 0, na.rm = TRUE), mean(rho_both < 0), NA))
print(out, row.names = FALSE)
cat("\nsubsets:\n"); print(round(subs, 4))
cat(sprintf("\nmodule sizes: median %.0f, largest %d; median selection-frequency binomial SE %.4f\n",
            median(lengths(ms)), max(lengths(ms)), median(sqrt(f * (1 - f) / NIT))))
write.table(out, "results/revision/selection_vs_AUC_dependence.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(subset = rownames(subs), subs), "results/revision/selection_vs_AUC_subsets.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
cat("===== 40b DONE =====\n")
