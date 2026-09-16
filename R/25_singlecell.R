#!/usr/bin/env Rscript
# 25_singlecell.R -- protocol section 21: single-cell localization in GSE278879
# section 21.1 per-sample QC (mitochondrial fraction is inherently low in snRNA-seq, so a 10% threshold is too permissive and must be set from the observed distribution)
# section 21.5: group comparisons must use patient-level pseudobulk (n=11); cells must not be treated as samples
# section 21.6: if NUMA1 is expressed mainly in fibroblasts/myeloid/endothelial cells, the conclusion becomes microenvironment-associated
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(Seurat); library(Matrix); library(dplyr); library(ggplot2)})
dir.create("results/singlecell/figures", showWarnings=FALSE, recursive=TRUE)
set.seed(42); options(future.globals.maxSize=3e9)
GX <- "NUMA1"

D <- "data_raw/GEO/GSE278879"
SAMP <- list.dirs(D, recursive=FALSE, full.names=FALSE)
SAMP <- SAMP[grepl("^GSM", SAMP)]
meta <- data.frame(dir=SAMP, sample=sub("^GSM\\d+_","",SAMP))
meta$group <- ifelse(grepl("^NE", meta$sample), "Normal", "Tumor")
meta$grade <- sub("\\d+$","", meta$sample)
cat("samples:\n"); print(meta)

## ---------- 1. read each sample and run QC ----------
qcstat <- data.frame()
objs <- list()
for (i in seq_len(nrow(meta))) {
  p <- file.path(D, meta$dir[i])
  f <- list.files(p, recursive=TRUE, full.names=TRUE)
  mtx <- grep("matrix.mtx", f, value=TRUE)[1]
  if (is.na(mtx)) { cat("!!", meta$sample[i], "no matrix.mtx, skipping\n"); next }
  dd <- dirname(mtx)
  m <- tryCatch(Read10X(dd), error=function(e) tryCatch(Read10X(dd, gene.column=1), error=function(e2) NULL))
  if (is.null(m)) { cat("!!", meta$sample[i], "failed to read\n"); next }
  if (is.list(m)) m <- m[[1]]
  o <- CreateSeuratObject(m, project=meta$sample[i], min.cells=3, min.features=200)
  o$sample <- meta$sample[i]; o$group <- meta$group[i]; o$grade <- meta$grade[i]
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern="^MT-")
  qcstat <- rbind(qcstat, data.frame(sample=meta$sample[i], cells_raw=ncol(o),
    med_nFeature=median(o$nFeature_RNA), med_nCount=median(o$nCount_RNA),
    med_mt=round(median(o$percent.mt),2), q95_mt=round(quantile(o$percent.mt,.95),2)))
  objs[[meta$sample[i]]] <- o
  cat(sprintf("  %-8s %6d cells | median nFeature %5.0f | median percent.mt %.2f%% (95%% %.2f%%)\n",
      meta$sample[i], ncol(o), median(o$nFeature_RNA), median(o$percent.mt), quantile(o$percent.mt,.95)))
}
write.table(qcstat, "results/singlecell/QC_per_sample.tsv", sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("\ntotal: %d samples, %d cells\n", length(objs), sum(sapply(objs, ncol))))

## section 21.1: set the mitochondrial threshold from the observed distribution (snRNA-seq is normally far below 10%)
allmt <- unlist(lapply(objs, function(o) o$percent.mt))
MT <- min(10, max(5, ceiling(quantile(allmt, .95))))
cat(sprintf("\npercent.mt across all cells: median %.2f%% | 95th percentile %.2f%% | 99th percentile %.2f%%\n",
    median(allmt), quantile(allmt,.95), quantile(allmt,.99)))
cat(sprintf("-> using the threshold percent.mt <= %g%% (protocol section 21.1 requires it to be set from the observed distribution and recorded)\n", MT))

objs <- lapply(objs, function(o) subset(o, nFeature_RNA>=200 & nFeature_RNA<=6000 & percent.mt<=MT))
cat(sprintf("%d cells after QC\n", sum(sapply(objs, ncol))))

## ---------- 2. merge and run the standard pipeline (with an intermediate cache so downstream changes need not rerun everything) ----------
CACHE <- "results/singlecell/_clustered.rds"
if (file.exists(CACHE)) {
  cat("\n[restoring clustering result from cache]\n"); obj <- readRDS(CACHE)
} else {
obj <- merge(objs[[1]], objs[-1], add.cell.ids=names(objs))
obj <- JoinLayers(obj)
obj <- NormalizeData(obj, verbose=FALSE) |> FindVariableFeatures(nfeatures=2000, verbose=FALSE)
obj <- ScaleData(obj, verbose=FALSE) |> RunPCA(npcs=30, verbose=FALSE)
## section 21.2: first check whether the unintegrated UMAP separates by sample
obj <- RunUMAP(obj, dims=1:20, verbose=FALSE, reduction.name="umap_raw")
p0 <- DimPlot(obj, reduction="umap_raw", group.by="sample") + ggtitle("Before integration")
## Harmony integration
suppressMessages(library(harmony))
obj <- RunHarmony(obj, group.by.vars="sample", verbose=FALSE)
obj <- RunUMAP(obj, dims=1:20, reduction="harmony", verbose=FALSE)
obj <- FindNeighbors(obj, dims=1:20, reduction="harmony", verbose=FALSE)
obj <- FindClusters(obj, resolution=0.5, verbose=FALSE)
cat("\nclusters:", length(levels(obj)), "\n")
p1 <- DimPlot(obj, group.by="sample") + ggtitle("After Harmony")
ggsave("results/singlecell/figures/UMAP_batch_check.pdf", p0|p1, width=13, height=5.5)
saveRDS(obj, CACHE); cat("[clustering result cached]\n")
}
cat("clusters:", length(levels(obj)), "\n")

## ---------- 3. section 21.3 primary lineage annotation ----------
MK <- list(Epithelial=c("EPCAM","KRT8","KRT18","KRT19"), `T/NK`=c("CD3D","CD3E","NKG7","GNLY"),
           `B/plasma`=c("CD79A","MS4A1","MZB1","JCHAIN"), Myeloid=c("LST1","TYROBP","CD68"),
           Fibroblast=c("COL1A1","COL1A2","DCN"), Endothelial=c("PECAM1","VWF","EMCN"))
nclust <- length(levels(obj))
sc <- sapply(MK, function(g) {
  g <- intersect(g, rownames(obj))
  if (!length(g)) return(rep(NA_real_, nclust))
  av <- AverageExpression(obj, features=g, group.by="seurat_clusters", layer="data")$RNA
  av <- log1p(as.matrix(av))                    # genes × clusters
  z  <- t(scale(t(av)))                          # standardise each gene across clusters
  m  <- colMeans(z, na.rm=TRUE)                  # one score per cluster
  if (length(m) != nclust) m <- rep(NA_real_, nclust)
  m })
sc <- matrix(as.numeric(sc), nrow=nclust, dimnames=list(levels(obj), names(MK)))
cat("\n=== lineage marker scores per cluster (standardised) ===\n"); print(round(sc,2))
# clusters whose top score falls below the threshold are called Unassigned rather than forced into a lineage (clusters 0, 21 and 23 score negatively on all six markers)
MINSCORE <- 0.25
best <- apply(sc, 1, function(r) if (all(is.na(r))) NA_character_ else colnames(sc)[which.max(r)])
bestv <- apply(sc, 1, function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm=TRUE))
ann <- ifelse(!is.na(bestv) & bestv >= MINSCORE, best, "Unassigned")
names(ann) <- levels(obj)
cat(sprintf("\nclusters called Unassigned for failing the score threshold %.2f: %s\n", MINSCORE,
    paste(names(ann)[ann=="Unassigned"], collapse=", ")))
obj$cell_type <- unname(ann[as.character(obj$seurat_clusters)])
cat("\n=== lineage annotation ===\n"); print(table(obj$cell_type))
cat("\nsample origin of each lineage:\n"); print(table(obj$cell_type, obj$group))
p2 <- DimPlot(obj, group.by="cell_type", label=TRUE) + ggtitle("Cell lineages")
ggsave("results/singlecell/figures/UMAP_celltype.pdf", p2, width=7.5, height=6)

## ---------- 4. section 21.4 localization of Gene X ----------
if (!GX %in% rownames(obj)) { cat("!!", GX, "is not in the matrix\n"); quit(save="no") }
ex <- FetchData(obj, vars=c(GX,"cell_type","group","sample"))
colnames(ex)[1] <- "expr"
loc <- ex |> group_by(cell_type) |>
  summarise(n_cells=n(), pct_expressing=round(100*mean(expr>0),1),
            mean_expr=round(mean(expr),3), mean_expr_in_pos=round(mean(expr[expr>0]),3)) |>
  arrange(desc(mean_expr))
cat("\n=== section 21.4 lineage distribution of NUMA1 (cell level, descriptive) ===\n"); print(as.data.frame(loc), row.names=FALSE)
write.table(loc, "results/singlecell/NUMA1_lineage_localization.tsv", sep="\t", row.names=FALSE, quote=FALSE)

ggsave("results/singlecell/figures/NUMA1_FeaturePlot.pdf", FeaturePlot(obj, features=GX), width=6.5, height=6)
ggsave("results/singlecell/figures/NUMA1_VlnPlot.pdf",
       VlnPlot(obj, features=GX, group.by="cell_type", pt.size=0) + NoLegend(), width=7, height=5)
ggsave("results/singlecell/figures/NUMA1_DotPlot.pdf",
       DotPlot(obj, features=GX, group.by="cell_type") + coord_flip(), width=6, height=3.5)

## ---------- 5. section 21.5 patient-level pseudobulk (the only valid basis for group comparison) ----------
cat("\n=== section 21.5 pseudobulk (patient level, n =", length(unique(ex$sample)), ") ===\n")
pb <- ex |> group_by(sample, group) |> summarise(mean_expr=mean(expr), .groups="drop")
print(as.data.frame(pb), row.names=FALSE)
w <- wilcox.test(mean_expr ~ group, data=pb)
cat(sprintf("\nall-cell pseudobulk, tumour vs normal: Wilcoxon P = %.3f (n=%d vs %d)\n",
    w$p.value, sum(pb$group=="Tumor"), sum(pb$group=="Normal")))
pbe <- ex |> filter(cell_type=="Epithelial") |> group_by(sample, group) |>
  summarise(mean_expr=mean(expr), n=n(), .groups="drop")
if (nrow(pbe) >= 6 && length(unique(pbe$group))==2) {
  we <- wilcox.test(mean_expr ~ group, data=pbe)
  cat(sprintf("epithelium-only pseudobulk, tumour vs normal: Wilcoxon P = %.3f (n=%d vs %d)\n",
      we$p.value, sum(pbe$group=="Tumor"), sum(pbe$group=="Normal")))
  print(as.data.frame(pbe), row.names=FALSE)
}
write.table(pb, "results/singlecell/NUMA1_pseudobulk.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## ---------- 6. section 21.6 decide how the conclusion must be phrased ----------
top <- loc$cell_type[1]
cat(sprintf("\n=== section 21.6 verdict ===\nlineage with the highest expression: %s\n", top))
if (top %in% c("Fibroblast","Myeloid","Endothelial")) {
  cat("-> the conclusion must be stated as a microenvironment-associated biomarker; tumor-intrinsic is NOT permitted\n")
} else if (top=="Epithelial") {
  cat("-> predominantly epithelial/tumour origin; tumor-intrinsic may be claimed, but the expression fraction in the other lineages must also be reported\n")
} else cat("-> predominantly immune/lymphoid origin; phrase the conclusion according to the observed lineage\n")

saveRDS(obj, "results/singlecell/GSE278879_seurat.rds")
writeLines(capture.output(sessionInfo()), "logs/25_singlecell_sessionInfo.txt")
cat("\n===== 25_singlecell DONE =====\n")
