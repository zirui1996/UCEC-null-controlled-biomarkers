#!/usr/bin/env Rscript
# 30_panel_immune.R -- rerun the section 18-19 immune analyses for the core 4 panel
# panel score = negative mean z-score of the four genes (all downregulated), the same definition used for the AUC in section 15
# the panel score is the primary variable; per-gene results go into a separate heatmap (avoiding the multiple-testing blow-up of 22x4 / 53x4 / 67x4 tests)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(ppcor); library(ggplot2); library(reshape2)})
set.seed(42)
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")

o  <- readRDS("results/immune/immune_objects.rds")      # lt = log2 TPM tumour matrix
pp <- readRDS("results/immune/immune_partial_objects.rds")
lt <- o$lt; ss <- pp$ss; prolif <- pp$prolif
g  <- intersect(CORE, rownames(lt))
cat("panel genes measurable in the TPM matrix:", paste(g, collapse=", "), "\n")
z  <- t(scale(t(lt[g,,drop=FALSE])))
panel <- -colMeans(z)                                   # all four genes are downregulated, so the sign is flipped: a high score is more tumour-like
cat(sprintf("panel score: n=%d, range %.2f to %.2f\n", length(panel), min(panel), max(panel)))
cat(sprintf("rho between the panel and the proliferation score = %.3f\n\n", cor(panel, prolif[names(panel)], method="spearman")))
per_gene <- lt[g,,drop=FALSE]

corblock <- function(mat, label) {
  r <- do.call(rbind, lapply(rownames(mat), function(f) {
    v <- as.numeric(mat[f,]); if (sd(v)==0) return(NULL)
    s  <- suppressWarnings(cor.test(panel, v, method="spearman"))
    pc <- pcor.test(rank(panel), rank(v), rank(prolif[names(panel)]))
    data.frame(feature=f, rho=unname(s$estimate), P=s$p.value,
               partial_rho=pc$estimate, partial_P=pc$p.value) }))
  r$BH_q <- p.adjust(r$P,"BH"); r$partial_BH_q <- p.adjust(r$partial_P,"BH")
  r$survives <- r$partial_BH_q<0.05 & sign(r$partial_rho)==sign(r$rho)
  r$block <- label; r[order(r$BH_q),] }

## ---- section 19.1 checkpoints ----
CP <- intersect(c("PDCD1","CD274","CTLA4","LAG3","TIGIT","HAVCR2"), rownames(lt))
rc <- corblock(lt[CP,,drop=FALSE], "checkpoint")
cat("=== section 19.1 immune checkpoints (panel score) ===\n")
print(data.frame(gene=rc$feature, rho=round(rc$rho,3), q=signif(rc$BH_q,3),
                 partial_rho=round(rc$partial_rho,3), partial_q=signif(rc$partial_BH_q,3),
                 survives=ifelse(rc$survives,"yes","no")), row.names=FALSE)

## ---- §19.2 ssGSEA ----
IMM <- grep("INTERFERON|INFLAMMATORY|ALLOGRAFT|IL6_JAK|COMPLEMENT|IL2_STAT5|TNFA_SIGNALING|antigen|cytotox|IFNg|T_cell_act",
            rownames(ss), value=TRUE)
rs <- corblock(ss[IMM,,drop=FALSE], "ssGSEA")
cat(sprintf("\n=== section 19.2 immune signatures (%d) ===\n", nrow(rs)))
print(data.frame(set=sub("^HALLMARK_","",rs$feature), rho=round(rs$rho,3), q=signif(rs$BH_q,3),
                 partial_rho=round(rs$partial_rho,3), partial_q=signif(rs$partial_BH_q,3),
                 survives=ifelse(rs$survives,"yes","no")), row.names=FALSE)
cat(sprintf("still significant with the same sign after controlling for proliferation: %d / %d\n", sum(rs$survives), nrow(rs)))

## ---- §18.5 xCell ----
xc <- pp$xcell[, names(panel), drop=FALSE]
rx <- corblock(xc, "xCell")
cat(sprintf("\n=== section 18.5 xCell: BH q<0.05 %d / %d | after controlling for proliferation %d ===\n",
    sum(rx$BH_q<0.05), nrow(rx), sum(rx$survives)))
print(head(data.frame(cell=rx$feature, rho=round(rx$rho,3), q=signif(rx$BH_q,3),
                      partial_rho=round(rx$partial_rho,3)), 15), row.names=FALSE)

## ---- §18.3 CIBERSORTx ----
cs <- read.delim("results/immune/cibersortx_out_Adjusted.txt", check.names=FALSE)
cs <- cs[cs$`P-value`<0.05, ]
CELLS <- setdiff(colnames(cs), c("Mixture","P-value","Correlation","RMSE","Absolute score (sig.score)"))
m <- match(cs$Mixture, names(panel)); ok <- !is.na(m)
csm <- t(as.matrix(cs[ok, CELLS])); colnames(csm) <- cs$Mixture[ok]
pan2 <- panel; prol2 <- prolif
panel <- panel[colnames(csm)]; prolif <- prolif[colnames(csm)]
rb <- corblock(csm[rowSums(csm>0)>=10,,drop=FALSE], "CIBERSORTx")
panel <- pan2; prolif <- prol2
cat(sprintf("\n=== section 18.3 CIBERSORTx (n=%d): BH q<0.05 %d / %d | after controlling for proliferation %d ===\n",
    sum(ok), sum(rb$BH_q<0.05), nrow(rb), sum(rb$survives)))
print(head(data.frame(cell=rb$feature, rho=round(rb$rho,3), q=signif(rb$BH_q,3),
                      partial_rho=round(rb$partial_rho,3)), 12), row.names=FALSE)

ALL <- rbind(rc, rs, rx, rb)
write.table(ALL, "results/immune/panel_immune_correlations.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## ---- per-gene heatmap (4 genes x immune features) ----
FE <- unique(c(CP, IMM, head(rx$feature[rx$BH_q<0.05],12)))
hm <- sapply(g, function(x) sapply(FE, function(f) {
  v <- if (f %in% rownames(lt)) as.numeric(lt[f,]) else
       if (f %in% rownames(ss)) as.numeric(ss[f,]) else as.numeric(xc[f,names(panel)])
  suppressWarnings(cor(as.numeric(lt[x,]), v, method="spearman")) }))
rownames(hm) <- sub("^HALLMARK_","",FE)
write.table(data.frame(feature=rownames(hm), hm), "results/immune/panel_per_gene_heatmap.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
d2 <- melt(hm); colnames(d2) <- c("feature","gene","rho")
p <- ggplot(d2, aes(gene, feature, fill=rho)) + geom_tile() +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0, limits=c(-0.45,0.45)) +
  labs(x=NULL, y=NULL, title="Per-gene immune correlations (core consensus panel)") +
  theme_minimal(base_size=9) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("results/immune/figures/panel_per_gene_immune_heatmap.pdf", p, width=6, height=max(5, nrow(hm)*0.22+1.5))
cat(sprintf("\nper-gene heatmap: %d immune features x %d genes\n", nrow(hm), ncol(hm)))
saveRDS(list(panel=panel, all=ALL, hm=hm), "results/immune/panel_immune_objects.rds")
cat("\n===== 30_panel_immune DONE =====\n")
