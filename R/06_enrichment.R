#!/usr/bin/env Rscript
# 06_enrichment.R -- Phase 3 / protocol sections 11.1-11.2: ORA (GO / KEGG / Reactome)
# background = the gene space shared by the three voting cohorts (13,039 genes), not the whole genome (section 11.1)
# foreground = recurrent consensus (v2.6: 182) / up / down / cell-cycle removed
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)})
dir.create("results/enrichment", showWarnings=FALSE, recursive=TRUE)
dir.create("results/enrichment/figures", showWarnings=FALSE, recursive=TRUE)
set.seed(42)

VOTE <- c("GSE17025","GSE63678","GSE106191")

## ---------- 1. background gene space (section 11.1) ----------
load1 <- function(g){t <- read.delim(sprintf("results/DEG/%s_all_genes.tsv",g), stringsAsFactors=FALSE); t[!duplicated(t$gene),]}
bg_sym <- Reduce(intersect, lapply(lapply(VOTE, load1), function(t) t$gene))
cat("background gene space (symbols):", length(bg_sym), "\n")

## ---------- 2. foreground sets ----------
rec  <- read.delim("results/consensus/recurrent_consensus.tsv", stringsAsFactors=FALSE)
core <- read.delim("results/consensus/core_consensus_3of3.tsv", stringsAsFactors=FALSE)

# cell-cycle genes (GO definition, identical to sections 9.3 and 35.7)
gos <- c("GO:0000278","GO:0007049","GO:0051301")
cc  <- unique(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=gos,
        keytype="GOALL", columns="SYMBOL"))$SYMBOL)

rec$cell_cycle <- rec$gene %in% cc
noCC <- rec[!rec$cell_cycle, ]
cat(sprintf("recurrent %d (up %d / down %d); cell-cycle %d; after removing cell-cycle genes %d (up %d / down %d)\n",
    nrow(rec), sum(rec$direction=="up"), sum(rec$direction=="down"), sum(rec$cell_cycle),
    nrow(noCC), sum(noCC$direction=="up"), sum(noCC$direction=="down")))

# output file required by section 10.1.4 (not previously generated)
write.table(rec,  "results/consensus/recurrent_consensus_annotated.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(noCC, "results/consensus/recurrent_noCellCycle.tsv",         sep="\t", row.names=FALSE, quote=FALSE)

SETS <- list(
  recurrent_all  = rec$gene,
  recurrent_up   = rec$gene[rec$direction=="up"],
  recurrent_down = rec$gene[rec$direction=="down"],
  noCellCycle    = noCC$gene,
  noCellCycle_up = noCC$gene[noCC$direction=="up"]
)

## ---------- 3. symbol -> ENTREZID (the mapping rate must be recorded faithfully) ----------
map <- function(sym) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=unique(sym),
        keytype="SYMBOL", columns="ENTREZID"))
  m <- m[!is.na(m$ENTREZID) & !duplicated(m$SYMBOL), ]
  m$ENTREZID
}
bg_ez <- map(bg_sym)
cat(sprintf("background mapping: %d / %d = %.1f%%\n", length(bg_ez), length(bg_sym), 100*length(bg_ez)/length(bg_sym)))
maprate <- data.frame()
SETS_EZ <- lapply(names(SETS), function(n) {
  e <- map(SETS[[n]])
  maprate <<- rbind(maprate, data.frame(set=n, n_symbol=length(SETS[[n]]), n_entrez=length(e),
                                        rate=round(100*length(e)/length(SETS[[n]]),1)))
  e
})
names(SETS_EZ) <- names(SETS)
print(maprate)
write.table(maprate, "results/enrichment/id_mapping_rates.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## ---------- 4. ORA ----------
run_ora <- function(ez, tag) {
  out <- list()
  for (ont in c("BP","CC","MF")) {
    r <- try(enrichGO(ez, universe=bg_ez, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont=ont,
                      pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE), silent=TRUE)
    out[[paste0("GO_",ont)]] <- if (inherits(r,"try-error")) NULL else r
  }
  r <- try(enrichKEGG(ez, universe=bg_ez, organism="hsa", pAdjustMethod="BH", pvalueCutoff=0.05), silent=TRUE)
  if (!inherits(r,"try-error") && !is.null(r)) r <- setReadable(r, org.Hs.eg.db, "ENTREZID")
  out[["KEGG"]] <- if (inherits(r,"try-error")) NULL else r
  if (requireNamespace("ReactomePA", quietly=TRUE) && requireNamespace("reactome.db", quietly=TRUE)) {
    r <- try(ReactomePA::enrichPathway(ez, universe=bg_ez, organism="human",
             pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE), silent=TRUE)
    out[["Reactome"]] <- if (inherits(r,"try-error")) NULL else r
  } else cat("  [Reactome skipped: reactome.db is not installed]\n")

  for (db in names(out)) {
    r <- out[[db]]
    n <- if (is.null(r)) 0 else sum(r@result$p.adjust < 0.05)
    cat(sprintf("  %-10s %-9s pathways at FDR<0.05: %d\n", tag, db, n))
    if (!is.null(r) && n > 0) {
      write.table(as.data.frame(r), sprintf("results/enrichment/%s_%s.tsv", tag, db),
                  sep="\t", row.names=FALSE, quote=FALSE)
      p <- dotplot(r, showCategory=min(20,n)) + ggtitle(sprintf("%s — %s", tag, db)) +
           theme(axis.text.y=element_text(size=8))
      ggsave(sprintf("results/enrichment/figures/%s_%s_dotplot.pdf", tag, db), p,
             width=9, height=max(4, min(20,n)*0.32+2), limitsize=FALSE)
    }
  }
  out
}

ALL <- list()
for (n in names(SETS_EZ)) { cat("\n=== ", n, " (n=", length(SETS_EZ[[n]]), ") ===\n", sep="")
  ALL[[n]] <- run_ora(SETS_EZ[[n]], n) }

saveRDS(list(ora=ALL, bg_ez=bg_ez, sets=SETS, maprate=maprate), "results/enrichment/ORA_objects.rds")
writeLines(capture.output(sessionInfo()), "logs/06_enrichment_sessionInfo.txt")
cat("\nKEGG access date:", format(Sys.Date()), "(KEGG is updated over time, so the Methods must record this)\n")
cat("===== 06_enrichment DONE =====\n")
