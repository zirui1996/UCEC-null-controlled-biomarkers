#!/usr/bin/env Rscript
# 33_panel_singlecell.R -- single-cell localization of the core 4 panel (section 21, using the cache)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(Seurat); library(dplyr); library(ggplot2); library(patchwork)})
set.seed(42)
CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1")
obj <- readRDS("results/singlecell/_clustered.rds")

# reuse the lineage annotation logic of section 21
MK <- list(Epithelial=c("EPCAM","KRT8","KRT18","KRT19"), `T/NK`=c("CD3D","CD3E","NKG7","GNLY"),
           `B/plasma`=c("CD79A","MS4A1","MZB1","JCHAIN"), Myeloid=c("LST1","TYROBP","CD68"),
           Fibroblast=c("COL1A1","COL1A2","DCN"), Endothelial=c("PECAM1","VWF","EMCN"))
nc <- length(levels(obj))
sc <- sapply(MK, function(g){ g <- intersect(g, rownames(obj)); if(!length(g)) return(rep(NA_real_,nc))
  av <- log1p(as.matrix(AverageExpression(obj, features=g, group.by="seurat_clusters", layer="data")$RNA))
  m <- colMeans(t(scale(t(av))), na.rm=TRUE); if(length(m)!=nc) rep(NA_real_,nc) else m })
sc <- matrix(as.numeric(sc), nrow=nc, dimnames=list(levels(obj), names(MK)))
best <- apply(sc,1,function(r) colnames(sc)[which.max(r)]); bv <- apply(sc,1,max,na.rm=TRUE)
ann <- ifelse(bv>=0.25, best, "Unassigned"); names(ann) <- levels(obj)
obj$cell_type <- unname(ann[as.character(obj$seurat_clusters)])

g <- intersect(CORE, rownames(obj))
cat("core 4 genes measurable in the single-cell matrix:", paste(g, collapse=", "), sprintf("(%d/4)\n\n", length(g)))
ex <- FetchData(obj, vars=c(g,"cell_type","group","sample"))

## ---- per-gene distribution across lineages ----
cat("=== lineage distribution of the core 4 genes (percent positive / mean expression) ===\n")
tab <- ex %>% group_by(cell_type) %>%
  summarise(n=n(), across(all_of(g), list(pct=~round(100*mean(.x>0),1), mean=~round(mean(.x),3)))) %>%
  arrange(desc(n))
print(as.data.frame(tab), row.names=FALSE)
write.table(tab, "results/singlecell/core4_lineage.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## ---- key test: does the panel score still differ between normal and tumour within epithelial cells ----
cat("\n=== §21.5 patient-level pseudobulk（n=11）===\n")
for (x in g) {
  pb <- ex %>% group_by(sample, group) %>% summarise(m=mean(.data[[x]]), .groups="drop")
  w <- wilcox.test(m ~ group, data=pb)
  pbe <- ex %>% filter(cell_type=="Epithelial") %>% group_by(sample, group) %>%
         summarise(m=mean(.data[[x]]), n=n(), .groups="drop")
  we <- if (length(unique(pbe$group))==2) wilcox.test(m ~ group, data=pbe)$p.value else NA
  cat(sprintf("  %-11s all cells P=%.3f | epithelium only P=%.3f | normal median %.3f vs tumour median %.3f\n",
      x, w$p.value, we, median(pb$m[pb$group=="Normal"]), median(pb$m[pb$group=="Tumor"])))
}

## ---- cell composition statistics supporting section 3.8 ----
cat("\n=== cell composition (supporting the stromal interpretation in section 3.8) ===\n")
ct <- table(obj$cell_type, obj$group); print(ct)
cat("\nproportions within each group (%):\n"); print(round(100*prop.table(ct,2),1))
comp <- ex %>% group_by(sample, group, cell_type) %>% summarise(n=n(), .groups="drop") %>%
  group_by(sample) %>% mutate(pct=100*n/sum(n)) %>% ungroup()
cat("\nWilcoxon test on the patient-level fibroblast fraction: ")
fb <- comp %>% filter(cell_type=="Fibroblast")
cat(sprintf("P = %.3f (normal median %.1f%% vs tumour median %.1f%%, n=%d vs %d)\n",
    wilcox.test(pct~group, data=fb)$p.value,
    median(fb$pct[fb$group=="Normal"]), median(fb$pct[fb$group=="Tumor"]),
    sum(fb$group=="Normal"), sum(fb$group=="Tumor")))
write.table(comp, "results/singlecell/cell_composition.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## ---- plots ----
ggsave("results/singlecell/figures/core4_FeaturePlot.pdf",
       FeaturePlot(obj, features=g, ncol=2) & theme(plot.title=element_text(size=10)), width=9, height=8)
ggsave("results/singlecell/figures/core4_VlnPlot.pdf",
       VlnPlot(obj, features=g, group.by="cell_type", pt.size=0, ncol=2) & NoLegend(), width=11, height=7)
ggsave("results/singlecell/figures/core4_DotPlot.pdf",
       DotPlot(obj, features=g, group.by="cell_type") + coord_flip(), width=7, height=4)
# saveRDS(obj, "results/singlecell/_annotated.rds")  # retired: it differs from _clustered.rds only by the cell_type column; the annotation logic is above and rebuilds in seconds
cat("\n===== 33_panel_singlecell DONE =====\n")
