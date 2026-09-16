#!/usr/bin/env Rscript
# 02d_GSE146889.R -- processing of the endometrial subset of GSE146889 and the section 8 admission test
# data: STAR + featureCounts raw counts (hg38), FFPE, Illumina HiSeq 2500
# keep only the 72 samples with source_name == "Endometrial" (excluding colorectal/ovarian/unlabelled tissue)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(limma); library(edgeR)})
set.seed(42)
dir.create("results/DEG", showWarnings=FALSE, recursive=TRUE)

## ---------- 1. sample table (parsed from SOFT, archived under metadata/) ----------
man <- read.delim("metadata/sample_manifest_GSE146889.tsv", stringsAsFactors=FALSE)
man <- subset(man, tissue=="Endometrial")
cat(sprintf("endometrial samples %d (tumour %d / normal %d); MSI %d / MSS %d\n", nrow(man),
    sum(man$group=="tumor"), sum(man$group=="normal"), sum(man$msi=="MSI"), sum(man$msi=="MSS")))

## ---------- 2. count matrix ----------
raw <- read.delim("data_raw/GEO/GSE146889/GSE146889_GeneCount.tsv.gz", check.names=FALSE)
cnt <- raw[, paste0(man$column, "_count"), drop=FALSE]
colnames(cnt) <- man$column
rownames(cnt) <- raw$GeneId
sym <- raw$GeneName
stopifnot(ncol(cnt)==nrow(man))

# for duplicate symbols keep the ENSG with the highest mean expression (consistent with the best-probe rule used in the other cohorts)
ord <- order(sym, -rowMeans(cnt))
keepg <- !duplicated(sym[ord]); idx <- ord[keepg]
m <- as.matrix(cnt[idx, ]); rownames(m) <- sym[idx]
m <- m[!is.na(rownames(m)) & rownames(m)!="", ]
cat("genes after collapsing duplicate symbols:", nrow(m), "\n")

## ---------- 3. low-expression filter + voom ----------
d <- DGEList(m, group=man$group)
keep <- filterByExpr(d, group=man$group)
d <- d[keep, , keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
cat("genes remaining after the low-expression filter:", nrow(d), "\n")

f <- factor(man$group, levels=c("normal","tumor"))
pat <- factor(man$pair_id)

run <- function(design, cm, label, permfun) {
  v  <- voom(d, design)
  fit<- eBayes(contrasts.fit(lmFit(v, design), cm))
  tt <- topTable(fit, number=Inf); tt$gene <- rownames(tt)
  n05 <- sum(tt$adj.P.Val<0.05)
  cat(sprintf("\n===== %s =====\n", label))
  cat(sprintf("criterion 1  FDR<0.05: %d | +|lfc|>=1: %d | +>=2: %d | smallest adj.P = %s\n",
      n05, sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>=1),
      sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>=2), signif(min(tt$adj.P.Val),3)))
  obs <- sum(tt$P.Value<0.001)
  perm <- replicate(1000, permfun())
  empP <- mean(perm >= obs)
  cat(sprintf("criterion 2  observed rawP<0.001 = %d | permuted median %.0f, 95%% percentile %.0f | empirical P = %.4f\n",
      obs, median(perm), quantile(perm,.95), empP))
  cat(sprintf("verdict   criterion 1 %s   criterion 2 %s\n", ifelse(n05>0,"pass","fail"), ifelse(empP<0.05,"pass","fail")))
  tt
}

## (a) paired analysis (the actual design of this dataset, 35 pairs)
paired <- man$pair_id %in% names(which(table(man$pair_id)==2))
cat("\ncomplete matched samples:", sum(paired), "(", sum(paired)/2, "pairs )\n")
dp <- d[, paired]; fp <- droplevels(f[paired]); pp <- droplevels(pat[paired])
Dp <- model.matrix(~0+fp+pp); colnames(Dp)[1:2] <- c("N","T")
cmP <- makeContrasts(T-N, levels=Dp)
permP <- function() {
  # paired permutation: swap tumor/normal labels at random within each pair
  g <- fp; for (p in levels(pp)) { i <- which(pp==p); if (runif(1)<0.5) g[i] <- rev(g[i]) }
  D <- model.matrix(~0+g+pp); colnames(D)[1:2] <- c("N","T")
  v <- voom(dp, D); ft <- eBayes(contrasts.fit(lmFit(v,D), makeContrasts(T-N, levels=D)))
  sum(topTable(ft, number=Inf)$P.Value < 0.001)
}
d_bak <- d; d <- dp
tt_paired <- run(Dp, cmP, "A paired analysis (35 pairs, the actual design of this dataset)", permP)
d <- d_bak

## (b) unpaired analysis (same definition as the other cohorts, for comparability)
Du <- model.matrix(~0+f); colnames(Du) <- c("N","T")
cmU <- makeContrasts(T-N, levels=Du)
permU <- function() {
  g <- sample(f); D <- model.matrix(~0+g); colnames(D) <- c("N","T")
  v <- voom(d, D); ft <- eBayes(contrasts.fit(lmFit(v,D), makeContrasts(T-N, levels=D)))
  sum(topTable(ft, number=Inf)$P.Value < 0.001)
}
tt_unpaired <- run(Du, cmU, "B unpaired analysis (37 vs 35, same definition as the other cohorts)", permU)

## (c) MSS-only sensitivity analysis (removing the MSI-H enrichment)
mss <- man$msi=="MSS"
cat(sprintf("\nMSS subset: tumour %d / normal %d\n", sum(mss & man$group=="tumor"), sum(mss & man$group=="normal")))
dm <- d[, mss]; fm <- droplevels(f[mss])
Dm <- model.matrix(~0+fm); colnames(Dm) <- c("N","T")
permM <- function() {
  g <- sample(fm); D <- model.matrix(~0+g); colnames(D) <- c("N","T")
  v <- voom(dm, D); ft <- eBayes(contrasts.fit(lmFit(v,D), makeContrasts(T-N, levels=D)))
  sum(topTable(ft, number=Inf)$P.Value < 0.001)
}
d_bak <- d; d <- dm
tt_mss <- run(Dm, makeContrasts(T-N, levels=Dm), "C MSS subset (9 vs 9, testing the effect of the MSI-H enrichment)", permM)
d <- d_bak

## ---------- 4. output ----------
sv <- function(tt, tag) {
  o <- tt[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val")]
  write.table(o, sprintf("results/DEG/GSE146889%s_all_genes.tsv", tag), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(subset(o, adj.P.Val<0.05 & abs(logFC)>=1),
              sprintf("results/DEG/GSE146889%s_DEG_primary.tsv", tag), sep="\t", row.names=FALSE, quote=FALSE)
}
sv(tt_paired, ""); sv(tt_unpaired, "_unpaired"); sv(tt_mss, "_MSSonly")
writeLines(capture.output(sessionInfo()), "logs/02d_GSE146889_sessionInfo.txt")
cat("\n===== GSE146889 DONE =====\n")
