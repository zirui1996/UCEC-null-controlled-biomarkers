#!/usr/bin/env Rscript
# 04_DEG.R -- Phase 2, Step 3
# independent limma differential expression per cohort. Protocol sections 9.2 / 9.5
# usage: Rscript R/04_DEG.R <GSE_ID> [scheme] [probe]
#   scheme: primary (default) | S1 | S2  -- S1/S2 apply to GSE115810 only (protocol section 9.4)
#   probe : best (default) | mean      -- sensitivity analysis of protocol section 9.1

args   <- commandArgs(trailingOnly = TRUE)
gse    <- args[1]
scheme <- ifelse(length(args) >= 2, args[2], "primary")
probe  <- ifelse(length(args) >= 3, args[3], "best")
proj   <- Sys.getenv("UCEC_PROJ", unset = getwd()); setwd(proj)
suppressPackageStartupMessages({ library(limma); library(ggplot2) })

cfg <- list(fdr = 0.05, lfc_primary = 1.0, lfc_high = 2.0)

man <- read.delim("metadata/sample_manifest.tsv", stringsAsFactors = FALSE, quote = "")
man <- man[man$dataset == gse & man$include_primary == "yes", ]
x   <- as.matrix(read.delim(sprintf("data_processed/expression/%s_gene_matrix_%sprobe.tsv",
                                    gse, probe), check.names = FALSE, row.names = 1))
man <- man[match(colnames(x), man$sample_id), ]
grp <- man$normal_or_tumor
lab <- sub("_.*$", "", man$title)

# ---- the two sensitivity schemes of protocol section 9.4 ----
disputed <- c("endometrium13", "endometrium16", "endometrium17", "endometrium24")
if (scheme == "S1") {                       # drop the 4 disputed samples
  stopifnot(gse == "GSE115810")
  keep <- !(lab %in% disputed)
  x <- x[, keep]; man <- man[keep, ]; grp <- grp[keep]; lab <- lab[keep]
} else if (scheme == "S2") {                # adopt the GSE36389 labels, using only the 20 shared samples
  stopifnot(gse == "GSE115810")
  g36_ctrl  <- c("endometrium4","endometrium5","endometrium6",
                 "endometrium13","endometrium16","endometrium17","endometrium24")
  g36_tumor <- c("endometrium2","endometrium3","endometrium7","endometrium8","endometrium9",
                 "endometrium12","endometrium14","endometrium15","endometrium18",
                 "endometrium19","endometrium20","endometrium21","endometrium22")
  keep <- lab %in% c(g36_ctrl, g36_tumor)
  x <- x[, keep]; man <- man[keep, ]; lab <- lab[keep]
  grp <- ifelse(lab %in% g36_ctrl, "normal", "tumor")
}

tag <- if (scheme == "primary" && probe == "best") gse else sprintf("%s_%s_%sprobe", gse, scheme, probe)
cat(sprintf("[%s] scheme=%s probe=%s | tumor=%d normal=%d genes=%d\n",
            gse, scheme, probe, sum(grp=="tumor"), sum(grp=="normal"), nrow(x)))
if (min(table(grp)) < 2) stop("a group has fewer than 2 samples; limma cannot be run")

# ---- limma (protocol section 9.2) ----
f <- factor(grp, levels = c("normal", "tumor"))
d <- model.matrix(~0 + f); colnames(d) <- c("Normal", "Tumor")
fit  <- lmFit(x, d)
fit2 <- eBayes(contrasts.fit(fit, makeContrasts(TumorVsNormal = Tumor - Normal, levels = d)))
tt   <- topTable(fit2, coef = "TumorVsNormal", number = Inf, adjust.method = "BH")
tt$gene <- rownames(tt)
tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

deg_p <- subset(tt, adj.P.Val < cfg$fdr & abs(logFC) >= cfg$lfc_primary)
deg_h <- subset(tt, adj.P.Val < cfg$fdr & abs(logFC) >= cfg$lfc_high)
cat(sprintf("  DEG primary (FDR<%.2f, |lfc|>=%.1f): %d  (up %d / down %d)\n",
            cfg$fdr, cfg$lfc_primary, nrow(deg_p),
            sum(deg_p$logFC > 0), sum(deg_p$logFC < 0)))
cat(sprintf("  DEG high-conf (|lfc|>=%.1f): %d\n", cfg$lfc_high, nrow(deg_h)))

od <- "results/DEG"; dir.create(od, recursive = TRUE, showWarnings = FALSE)
w <- function(df, suf) write.table(df, file.path(od, sprintf("%s_%s.tsv", tag, suf)),
                                   sep = "\t", row.names = FALSE, quote = FALSE)
w(tt, "all_genes"); w(deg_p, "DEG_primary"); w(deg_h, "DEG_high_conf")
w(subset(deg_p, logFC > 0), "DEG_up"); w(subset(deg_p, logFC < 0), "DEG_down")

# ---- volcano plot + top-50 heatmap (material for Figure 2) ----
tt$sig <- ifelse(tt$adj.P.Val < cfg$fdr & tt$logFC >=  cfg$lfc_primary, "Up",
          ifelse(tt$adj.P.Val < cfg$fdr & tt$logFC <= -cfg$lfc_primary, "Down", "NS"))
pdf(file.path("figures", sprintf("%s_volcano.pdf", tag)), width = 6, height = 5)
print(ggplot(tt, aes(logFC, -log10(adj.P.Val), colour = sig)) +
  geom_point(size = .6, alpha = .6) +
  scale_colour_manual(values = c(Up = "#C0392B", Down = "#2471A3", NS = "grey80")) +
  geom_vline(xintercept = c(-1, 1) * cfg$lfc_primary, linetype = 2, linewidth = .3) +
  geom_hline(yintercept = -log10(cfg$fdr), linetype = 2, linewidth = .3) +
  labs(title = tag, subtitle = sprintf("Up %d / Down %d (FDR<0.05, |log2FC|>=1)",
       sum(tt$sig=="Up"), sum(tt$sig=="Down")),
       x = "log2 fold change", y = "-log10 adjusted P") + theme_classic())
dev.off()

if (nrow(deg_p) >= 10) {
  top <- head(deg_p$gene[order(-abs(deg_p$logFC))], 50)
  hm  <- t(scale(t(x[top, , drop = FALSE])))
  hm[hm >  2] <-  2; hm[hm < -2] <- -2
  pdf(file.path("figures", sprintf("%s_heatmap_top50.pdf", tag)), width = 9, height = 8)
  heatmap(hm, scale = "none", margins = c(7, 6), cexRow = .5, cexCol = .45,
          col = colorRampPalette(c("#2471A3", "white", "#C0392B"))(64),
          ColSideColors = c(normal = "grey60", tumor = "#C0392B")[grp],
          main = sprintf("%s - top 50 DEG", tag))
  dev.off()
}
writeLines(capture.output(sessionInfo()), sprintf("logs/04_DEG_%s_sessionInfo.txt", tag))
cat("[", tag, "] DONE\n")
