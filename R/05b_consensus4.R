#!/usr/bin/env Rscript
# 05b_consensus4.R -- consensus over four voting cohorts (candidate rule of v2.4)
# voting: GSE17025 / GSE63678 / GSE106191 / GSE146889
# veto: GSE56087 (did not pass admission; contributes only the sign of its logFC)
# compute the >=4/4, >=3/4 and >=2/4 tiers together; the UpSet plot shows the full spectrum and the counts decide the final rule
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(UpSetR); library(org.Hs.eg.db)})
dir.create("results/consensus/figures", showWarnings=FALSE, recursive=TRUE)

VOTE <- c("GSE17025","GSE63678","GSE106191","GSE146889")
REF  <- "GSE56087"; FDR <- 0.05; LFC <- 1.0

load1 <- function(g){t <- read.delim(sprintf("results/DEG/%s_all_genes.tsv",g), stringsAsFactors=FALSE); t[!duplicated(t$gene),]}
D <- lapply(c(VOTE,REF), load1); names(D) <- c(VOTE,REF)
for (g in c(VOTE,REF)) cat(sprintf("  %-11s %6d genes | FDR<.05 %5d | +|lfc|>=1 %5d\n", g, nrow(D[[g]]),
    sum(D[[g]]$adj.P.Val<FDR), sum(D[[g]]$adj.P.Val<FDR & abs(D[[g]]$logFC)>=LFC)))

uni <- Reduce(intersect, lapply(D[VOTE], function(t) t$gene))
cat("\ngene space shared by the four cohorts:", length(uni), "(13,039 with three cohorts)\n")

sigmat <- sapply(VOTE, function(g){t<-D[[g]];rownames(t)<-t$gene;t<-t[uni,];(t$adj.P.Val<FDR)&(abs(t$logFC)>=LFC)})
lfcmat <- sapply(VOTE, function(g){t<-D[[g]];rownames(t)<-t$gene;t[uni,"logFC"]})
rownames(sigmat) <- rownames(lfcmat) <- uni
nvote <- rowSums(sigmat)
dir_ok <- sapply(seq_len(nrow(sigmat)), function(i){
  s <- which(sigmat[i,]); if(!length(s)) return(NA)
  length(unique(sign(lfcmat[i,s])))==1 })

cat("\n=== UpSet counts (distribution of the number of significant cohorts) ===\n"); print(table(n_significant_cohorts=nvote))
cat("\n=== direction concordance filter ===\n")
for (k in 4:1) cat(sprintf("  %d/4: same direction %5d | removed for conflicting direction %3d\n", k,
    sum(nvote==k & dir_ok %in% TRUE), sum(nvote==k & dir_ok %in% FALSE)))

# GSE56087 direction veto
r <- D[[REF]]; rownames(r) <- r$gene
veto <- function(gs){ ok <- intersect(gs, rownames(r))
  ok[sign(r[ok,"logFC"]) != sign(rowMeans(lfcmat[ok,,drop=FALSE]))] }

tiers <- list()
for (k in c(4,3,2)) {
  g0 <- uni[nvote>=k & dir_ok %in% TRUE]
  v  <- veto(g0); g1 <- setdiff(g0, v)
  nm <- intersect(g0, rownames(r))
  cat(sprintf("\n>=%d/4 : same direction %d -> measurable in GSE56087 %d, vetoed for opposite sign %d (concordance %.1f%%) -> final %d\n",
      k, length(g0), length(nm), length(v), 100*(length(nm)-length(v))/max(length(nm),1), length(g1)))
  tiers[[as.character(k)]] <- g1
}

core <- tiers[["4"]]; recur <- tiers[["3"]]
cat(sprintf("\n>>> core consensus (4/4)      : %d\n", length(core)))
cat(sprintf(">>> recurrent consensus (>=3/4): %d (up %d / down %d)\n", length(recur),
    sum(rowMeans(lfcmat[recur,,drop=FALSE])>0), sum(rowMeans(lfcmat[recur,,drop=FALSE])<0)))
cat(sprintf(">>> relaxing to >=2/4 would give: %d (for reference only, not adopted)\n", length(tiers[["2"]])))

cat("\n=== composition of the significant combinations in the >=3/4 set ===\n")
combo <- apply(sigmat[recur,,drop=FALSE], 1, function(x) paste(VOTE[x], collapse="+"))
print(sort(table(combo), decreasing=TRUE))

# cell-cycle fraction (same GO definition as sections 9.3 and 35.7)
cc <- unique(suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
        keys=c("GO:0000278","GO:0007049","GO:0051301"), keytype="GOALL", columns="SYMBOL"))$SYMBOL)
cat("\n=== cell-cycle fraction ===\n")
for (nm in c("4","3","2")) {
  g <- tiers[[nm]]; up <- g[rowMeans(lfcmat[g,,drop=FALSE])>0]; dn <- g[rowMeans(lfcmat[g,,drop=FALSE])<0]
  cat(sprintf("  >=%s/4  cell-cycle genes: %d of %d upregulated (%.1f%%) | %d of %d downregulated (%.1f%%)\n", nm,
      length(up), sum(up %in% cc), 100*mean(up %in% cc), length(dn), sum(dn %in% cc), 100*mean(dn %in% cc)))
}

# output
out <- data.frame(gene=recur, n_significant=nvote[recur], combo=combo[recur],
                  lfcmat[recur,,drop=FALSE],
                  GSE56087_logFC=r[match(recur, rownames(r)),"logFC"],
                  direction=ifelse(rowMeans(lfcmat[recur,,drop=FALSE])>0,"up","down"),
                  cell_cycle=recur %in% cc,
                  tier=ifelse(recur %in% core, "core_4of4","recurrent_3of4"))
out <- out[order(-out$n_significant, out$gene), ]
write.table(out, "results/consensus/v4_recurrent_consensus_3of4.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(subset(out, tier=="core_4of4"), "results/consensus/v4_core_consensus_4of4.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# UpSet
sets <- lapply(VOTE, function(g) uni[sigmat[,g]]); names(sets) <- VOTE
pdf("results/consensus/figures/UpSet_consensus_4cohort.pdf", width=10, height=6.5, onefile=FALSE)
print(upset(fromList(sets), nsets=4, order.by="freq", nintersects=15,
      mainbar.y.label="Genes in intersection", sets.x.label="Significant genes per cohort",
      text.scale=1.2))
dev.off()
saveRDS(list(sigmat=sigmat,lfcmat=lfcmat,nvote=nvote,dir_ok=dir_ok,tiers=tiers), "results/consensus/v4_objects.rds")
writeLines(capture.output(sessionInfo()), "logs/05b_consensus4_sessionInfo.txt")
cat("\n===== 05b_consensus4 DONE =====\n")
