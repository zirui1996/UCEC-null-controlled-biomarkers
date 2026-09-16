#!/usr/bin/env Rscript
# 05_consensus.R -- reproducibility consensus (>=2/3 rule) + UpSet
# v2.6: the GSE56087 direction veto was withdrawn (it rested on noise-level signs with |log2FC|<=0.2 and adj.P>=0.88);
#       GSE56087 is used only for a supplementary direction-concordance check and affects no gene set
# voting cohorts: GSE17025 / GSE63678 / GSE106191 (all passed admission criteria 1 and 2)
# GSE56087 did not pass admission: it does not vote and does not veto; only its direction concordance is reported (Supplementary)
proj <- Sys.getenv("UCEC_PROJ", unset = getwd()); setwd(proj)

VOTE <- c("GSE17025","GSE63678","GSE106191")
REF  <- "GSE56087"
FDR  <- 0.05; LFC <- 1.0

load1 <- function(g) {
  t <- read.delim(sprintf("results/DEG/%s_all_genes.tsv", g), stringsAsFactors=FALSE)
  t[!duplicated(t$gene), ]
}
D <- lapply(c(VOTE, REF), load1); names(D) <- c(VOTE, REF)

# shared gene space (voting cohorts only)
uni <- Reduce(intersect, lapply(D[VOTE], function(t) t$gene))
cat("gene space shared by the voting cohorts:", length(uni), "\n")
for (g in VOTE) cat(sprintf("  %-11s %d genes\n", g, nrow(D[[g]])))

sigmat <- sapply(VOTE, function(g) {
  t <- D[[g]]; rownames(t) <- t$gene; t <- t[uni,]
  (t$adj.P.Val < FDR) & (abs(t$logFC) >= LFC)
})
lfcmat <- sapply(VOTE, function(g) { t <- D[[g]]; rownames(t) <- t$gene; t[uni,"logFC"] })
rownames(sigmat) <- rownames(lfcmat) <- uni

nvote <- rowSums(sigmat)
# same direction: every significant effect has the same sign
dir_ok <- sapply(seq_len(nrow(sigmat)), function(i) {
  s <- which(sigmat[i,]); if (!length(s)) return(NA)
  length(unique(sign(lfcmat[i, s]))) == 1
})

cat("\n=== UpSet counts (distribution of the number of significant cohorts) ===\n")
print(table(n_significant_cohorts=nvote))
cat("\n=== after the direction concordance filter ===\n")
for (k in 3:1) {
  ok <- which(nvote==k & dir_ok %in% TRUE)
  bad <- sum(nvote==k & dir_ok %in% FALSE)
  cat(sprintf("  significant in %d/3: %4d with the same direction | %d removed for conflicting direction\n", k, length(ok), bad))
}

core <- uni[nvote==3 & dir_ok %in% TRUE]                 # 3/3
recur<- uni[nvote>=2 & dir_ok %in% TRUE]                 # >=2/3
cat(sprintf("\ncore consensus  (3/3)  : %d genes\n", length(core)))
cat(sprintf("recurrent consensus (>=2/3): %d genes (up %d / down %d)\n", length(recur),
            sum(rowMeans(lfcmat[recur,,drop=FALSE])>0), sum(rowMeans(lfcmat[recur,,drop=FALSE])<0)))

# contribution of each pairwise combination of cohorts
cat("\n=== composition of the significant combinations in the >=2/3 set ===\n")
combo <- apply(sigmat[recur,,drop=FALSE], 1, function(r) paste(VOTE[r], collapse="+"))
print(sort(table(combo), decreasing=TRUE))

# supplementary direction-concordance check with GSE56087 (v2.6: no gene is removed)
r <- D[[REF]]; rownames(r) <- r$gene
ok <- intersect(recur, rownames(r))
conc <- sign(r[ok,"logFC"]) == sign(rowMeans(lfcmat[ok,,drop=FALSE]))
cat("\n=== supplementary direction concordance with GSE56087 (no vote, no veto) ===\n")
cat(sprintf("  of the %d genes significant in >=2/3 cohorts, %d are measurable in GSE56087; %d share the direction (%.1f%%); binomial P = %.2g\n",
            length(recur), length(ok), sum(conc), 100*mean(conc), binom.test(sum(conc), length(ok), 0.5)$p.value))
cat(sprintf("  discordant direction: %s\n", paste(sprintf("%s(%.2f)", ok[!conc], r[ok[!conc],"logFC"]), collapse=", ")))
cat(sprintf("  final: core (3/3) = %d genes, recurrent (>=2/3) = %d genes\n", length(core), length(recur)))
write.table(data.frame(gene=ok, discovery_mean_logFC=rowMeans(lfcmat[ok,,drop=FALSE]), GSE56087_logFC=r[ok,"logFC"],
                       GSE56087_adjP=r[ok,"adj.P.Val"], concordant=conc),
            "results/consensus/GSE56087_direction_check.tsv", sep="\t", row.names=FALSE, quote=FALSE)

combo <- apply(sigmat[recur,,drop=FALSE], 1, function(rr) paste(VOTE[rr], collapse="+"))
out <- data.frame(gene=recur, n_significant=nvote[recur],
                  combo=combo[recur],
                  GSE17025_logFC=lfcmat[recur,"GSE17025"],
                  GSE63678_logFC=lfcmat[recur,"GSE63678"],
                  GSE106191_logFC=lfcmat[recur,"GSE106191"],
                  GSE56087_logFC=r[match(recur, rownames(r)),"logFC"],
                  direction=ifelse(rowMeans(lfcmat[recur,,drop=FALSE])>0,"up","down"),
                  tier=ifelse(recur %in% core, "core_3of3","recurrent_2of3"),
                  stringsAsFactors=FALSE)
out <- out[order(-out$n_significant, out$direction, out$gene), ]
dir.create("results/consensus", showWarnings=FALSE, recursive=TRUE)
write.table(out, "results/consensus/recurrent_consensus.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(subset(out, tier=="core_3of3"), "results/consensus/core_consensus_3of3.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

# UpSet plot
if (requireNamespace("UpSetR", quietly=TRUE)) {
  lst <- lapply(VOTE, function(g) uni[sigmat[,g]]); names(lst) <- VOTE
  pdf("figures/UpSet_consensus.pdf", width=8, height=5.5)
  print(UpSetR::upset(UpSetR::fromList(lst), order.by="freq", nsets=3,
                      mainbar.y.label="Genes", sets.x.label="DEGs per cohort"))
  dev.off(); cat("\nUpSet plot -> figures/UpSet_consensus.pdf\n")
} else cat("\n[note] UpSetR is not installed, skipping the UpSet plot\n")
cat("\n===== CONSENSUS DONE =====\n")
