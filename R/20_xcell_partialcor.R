#!/usr/bin/env Rscript
# 20_xcell_partialcor.R
#  (a) section 18.5 xCell -- replication at the level of direction only; absolute values are not compared with CIBERSORTx
#  (b) proliferation confounding check -- partial correlation controlling for the proliferation score (the issue raised in section 9.3, carried into the immune module)
#  (c) restore the T_cell_activation signature that maxSize had filtered out
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(xCell); library(ppcor); library(GSVA); library(org.Hs.eg.db)})
set.seed(42); GX <- "NUMA1"
o  <- readRDS("results/immune/immune_objects.rds")
lt <- o$lt; gx <- o$gx; ss <- o$ss

## ---------- 1. proliferation score (defined before analysis: mean ssGSEA of the two Hallmark proliferation gene sets) ----------
PROLIF <- c("HALLMARK_MITOTIC_SPINDLE","HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS","HALLMARK_MYC_TARGETS_V1")
PROLIF <- intersect(PROLIF, rownames(ss))
prolif <- colMeans(ss[PROLIF, , drop=FALSE])
cat("proliferation score = mean ssGSEA of the following Hallmark sets:\n  ", paste(sub("^HALLMARK_","",PROLIF), collapse=", "), "\n")
cat(sprintf("Spearman rho between NUMA1 and the proliferation score = %.3f (P = %.2g)\n\n",
    cor(gx, prolif, method="spearman"), cor.test(gx, prolif, method="spearman")$p.value))

## ---------- 2. restore T_cell_activation (relaxing maxSize) ----------
tca <- unique(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys="GO:0042110",
        keytype="GOALL", columns="SYMBOL"))$SYMBOL)
tca <- intersect(tca, rownames(lt))
cat(sprintf("T_cell_activation (GO:0042110) has %d genes -- filtered out by maxSize=500 in the previous run, now relaxed to 2000\n", length(tca)))
p2 <- ssgseaParam(as.matrix(lt), list(T_cell_activation=tca), minSize=10, maxSize=2000, normalize=TRUE)
ss2 <- gsva(p2, verbose=FALSE)
ss <- rbind(ss, ss2)

## ---------- 3. partial correlation: do the immune associations survive after controlling for the proliferation score ----------
CP <- c("PDCD1","CD274","CTLA4","LAG3","TIGIT","HAVCR2")
IMM <- c("HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INFLAMMATORY_RESPONSE",
         "HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_COMPLEMENT","antigen_presentation","cytotoxicity",
         "IFNg_response","T_cell_activation")
IMM <- intersect(IMM, rownames(ss))
mk <- function(nm, v) {
  s <- suppressWarnings(cor.test(gx, v, method="spearman"))
  pc <- pcor.test(rank(gx), rank(v), rank(prolif), method="pearson")   # Spearman partial correlation = Pearson partial correlation on ranks
  data.frame(feature=nm, rho=unname(s$estimate), P=s$p.value,
             partial_rho=pc$estimate, partial_P=pc$p.value) }
res <- rbind(do.call(rbind, lapply(CP,  function(g) mk(g, lt[g,]))),
             do.call(rbind, lapply(IMM, function(g) mk(g, ss[g,]))))
res$BH_q <- p.adjust(res$P,"BH"); res$partial_BH_q <- p.adjust(res$partial_P,"BH")
res$survives <- res$partial_BH_q < 0.05 & sign(res$partial_rho)==sign(res$rho)
cat("\n=== immune associations before and after controlling for the proliferation score ===\n")
print(data.frame(feature=sub("^HALLMARK_","",res$feature), rho=round(res$rho,3), q=signif(res$BH_q,3),
                 partial_rho=round(res$partial_rho,3), partial_q=signif(res$partial_BH_q,3),
                 survives=ifelse(res$survives,"yes","no")), row.names=FALSE)
cat(sprintf("\nstill significant with the same sign after controlling for proliferation: %d / %d\n", sum(res$survives), nrow(res)))
write.table(res, "results/immune/partial_correlation_proliferation.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## ---------- 4. §18.5 xCell ----------
cat("\n=== section 18.5 xCell (replication at the level of direction only) ===\n")
xc <- xCellAnalysis(as.matrix(2^lt - 1), rnaseq=TRUE, parallel.sz=2)
cat("xCell output:", nrow(xc), "scores x", ncol(xc), "samples\n")
xr <- t(sapply(rownames(xc), function(r) {
  s <- suppressWarnings(cor.test(gx, xc[r,], method="spearman")); c(rho=unname(s$estimate), P=s$p.value)}))
xr <- as.data.frame(xr); xr$BH_q <- p.adjust(xr$P,"BH"); xr$cell <- rownames(xr)
xr$partial_rho <- sapply(rownames(xc), function(r) pcor.test(rank(gx), rank(xc[r,]), rank(prolif))$estimate)
xr$partial_P   <- sapply(rownames(xc), function(r) pcor.test(rank(gx), rank(xc[r,]), rank(prolif))$p.value)
xr$partial_BH_q <- p.adjust(xr$partial_P,"BH")
xr <- xr[order(xr$BH_q),]
write.table(xr[,c("cell","rho","P","BH_q","partial_rho","partial_P","partial_BH_q")],
            "results/immune/xCell_correlations.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("BH q<0.05: %d / %d | still significant after controlling for proliferation: %d\n", sum(xr$BH_q<0.05), nrow(xr),
    sum(xr$partial_BH_q<0.05 & sign(xr$partial_rho)==sign(xr$rho))))
cat("\n--- xCell Top 15 ---\n")
print(data.frame(cell=head(xr$cell,15), rho=round(head(xr$rho,15),3), q=signif(head(xr$BH_q,15),3),
                 partial_rho=round(head(xr$partial_rho,15),3), partial_q=signif(head(xr$partial_BH_q,15),3)),
      row.names=FALSE)

saveRDS(list(prolif=prolif, res=res, xcell=xc, xr=xr, ss=ss), "results/immune/immune_partial_objects.rds")
writeLines(capture.output(sessionInfo()), "logs/20_xcell_partial_sessionInfo.txt")
cat("\n===== 20_xcell_partialcor DONE =====\n")
