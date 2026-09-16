#!/usr/bin/env Rscript
# 31_stroma_adjusted.R -- stromal confounding check for the immune associations of the core 4 panel
# motivation: the strongest association of the panel is with the xCell StromaScore (-0.426), stronger than with any immune feature
#       -> immune associations may be only a downstream effect of "less stroma / more epithelium"
# approach: partial correlation controlling for StromaScore, and a second-order partial correlation controlling for StromaScore plus the proliferation score
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(ppcor); library(ggplot2)})
set.seed(42)
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")

o  <- readRDS("results/immune/immune_objects.rds")
pp <- readRDS("results/immune/immune_partial_objects.rds")
lt <- o$lt; ss <- pp$ss; xc <- pp$xcell; prolif <- pp$prolif
g <- intersect(CORE, rownames(lt))
panel <- -colMeans(t(scale(t(lt[g,,drop=FALSE]))))
S <- as.numeric(xc["StromaScore", names(panel)])
P <- as.numeric(prolif[names(panel)])
cat(sprintf("panel vs StromaScore rho = %.3f | panel vs proliferation rho = %.3f | stroma vs proliferation rho = %.3f\n\n",
    cor(panel,S,method="spearman"), cor(panel,P,method="spearman"), cor(S,P,method="spearman")))

CP  <- intersect(c("PDCD1","CD274","CTLA4","LAG3","TIGIT","HAVCR2"), rownames(lt))
IMM <- grep("INTERFERON|INFLAMMATORY|ALLOGRAFT|IL6_JAK|COMPLEMENT|IL2_STAT5|TNFA_SIGNALING|antigen|cytotox|IFNg|T_cell_act",
            rownames(ss), value=TRUE)
XC  <- c("ImmuneScore","CD8+ T-cells","CD8+ Tem","Th1 cells","NK cells","Macrophages M1",
         "Macrophages M2","DC","pDC","B-cells","Plasma cells","Tregs")
XC  <- intersect(XC, rownames(xc))

getv <- function(f) {
  if (f %in% rownames(lt)) as.numeric(lt[f,names(panel)])
  else if (f %in% rownames(ss)) as.numeric(ss[f,names(panel)])
  else as.numeric(xc[f,names(panel)]) }

res <- do.call(rbind, lapply(c(CP, IMM, XC), function(f) {
  v <- getv(f); if (all(is.na(v)) || sd(v,na.rm=TRUE)==0) return(NULL)
  r0 <- suppressWarnings(cor.test(panel, v, method="spearman"))
  p1 <- pcor.test(rank(panel), rank(v), rank(S))                       # controlling for stroma
  p2 <- pcor.test(rank(panel), rank(v), data.frame(rank(S), rank(P)))  # controlling for stroma + proliferation
  data.frame(feature=f,
             rho=unname(r0$estimate), P=r0$p.value,
             rho_S=p1$estimate, P_S=p1$p.value,
             rho_SP=p2$estimate, P_SP=p2$p.value) }))
res$q     <- p.adjust(res$P,"BH")
res$q_S   <- p.adjust(res$P_S,"BH")
res$q_SP  <- p.adjust(res$P_SP,"BH")
res$survives_S <- res$q_S<0.05 & sign(res$rho_S)==sign(res$rho)
res <- res[order(res$q),]

cat("=== before and after controlling for StromaScore (BH-adjusted) ===\n")
cat(sprintf("%-30s %7s %9s | %7s %9s | %7s %9s  %s\n",
    "feature","rho","q","rho|S","q|S","rho|S,P","q|S,P","survives"))
for (i in seq_len(nrow(res))) cat(sprintf("%-30s %+7.3f %9.2e | %+7.3f %9.2e | %+7.3f %9.2e  %s\n",
    substr(sub("^HALLMARK_","",res$feature[i]),1,30), res$rho[i], res$q[i],
    res$rho_S[i], res$q_S[i], res$rho_SP[i], res$q_SP[i],
    ifelse(res$survives_S[i],"yes","no")))

n0 <- sum(res$q<0.05); n1 <- sum(res$survives_S); n2 <- sum(res$q_SP<0.05 & sign(res$rho_SP)==sign(res$rho))
cat(sprintf("\nsignificant before adjustment: %d / %d\nstill significant with the same sign after controlling for stroma: %d\nafter controlling for stroma and proliferation: %d\n",
    n0, nrow(res), n1, n2))
cat(sprintf("\nverdict: %s\n", ifelse(n1 >= max(3, 0.4*n0),
  "immune associations are largely independent of stromal content -> the immune arm can be retained",
  "most immune associations vanish after controlling for stroma -> report as a difference in cell composition, not as an independent immune conclusion")))
write.table(res, "results/immune/panel_stroma_adjusted.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat("\n===== 31_stroma_adjusted DONE =====\n")
