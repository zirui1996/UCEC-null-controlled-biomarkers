#!/usr/bin/env Rscript
# 26_mutation.R -- protocol section 22.1 mutation analysis (Supplementary)
# mandatory precondition: the POLE ultramutated subtype has extremely high TMB and would dominate any mutation-frequency comparison, so stratification or adjustment is required
# forbidden: picking genes by eye from a waterfall plot (every gene gets a Fisher test with BH-FDR)
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(maftools); library(dplyr)})
dir.create("results/mutation/figures", showWarnings=FALSE, recursive=TRUE)
set.seed(42); GX <- "NUMA1"

maf <- readRDS("data_raw/TCGA/mut_cnv/TCGA_UCEC_maf.rds")
maf <- as.data.frame(maf)
maf$pid <- substr(maf$Tumor_Sample_Barcode, 1, 12)
cat(sprintf("MAF: %d rows, %d patients\n", nrow(maf), length(unique(maf$pid))))

## ---- NUMA1 expression groups (same analysis set as sections 16 and 20) ----
o  <- readRDS("results/immune/immune_objects.rds"); gx <- o$gx
sv <- readRDS("results/survival/survival_objects.rds")
names(gx) <- substr(names(gx), 1, 12)
gx <- gx[!duplicated(names(gx))]
cat("patients with a NUMA1 expression value:", length(gx), "\n")

## ---- molecular subtypes (used for the mandatory stratification of section 22.1) ----
cd <- read.delim("metadata/TCGA_UCEC_colData.tsv", stringsAsFactors=FALSE)
cd$pid <- substr(cd$barcode,1,12)
sub <- cd$paper_IntegrativeCluster[match(names(gx), cd$pid)]
sub <- ifelse(sub %in% c("CN high","CN low","MSI","POLE"), sub, NA)
cat("samples with a molecular subtype:", sum(!is.na(sub)), "/", length(gx), "\n")
cat("  "); print(table(sub, useNA="no"))

## ---- relationship between TMB and NUMA1 (first assess how severe the POLE confounding is) ----
tmb <- table(maf$pid)
common <- intersect(names(gx), names(tmb))
d <- data.frame(pid=common, gx=gx[common], tmb=as.numeric(tmb[common]),
                subtype=sub[match(common, names(gx))])
cat(sprintf("\n=== TMB distribution (mutations per patient) ===\nmedian %.0f | range %d-%d\n",
    median(d$tmb), min(d$tmb), max(d$tmb)))
print(d %>% filter(!is.na(subtype)) %>% group_by(subtype) %>%
      summarise(n=n(), median_TMB=median(tmb), max_TMB=max(tmb)) %>% as.data.frame(), row.names=FALSE)
ct <- cor.test(d$gx, log10(d$tmb+1), method="spearman")
cat(sprintf("\nNUMA1 vs log10(TMB): rho = %.3f, P = %.3g\n", ct$estimate, ct$p.value))

## ---- section 22.1 Fisher test with BH adjustment for every gene (stratified by molecular subtype) ----
d$grp <- ifelse(d$gx > median(d$gx), "High", "Low")
GENES <- names(which(table(unique(maf[,c("pid","Hugo_Symbol")])$Hugo_Symbol) >= 0.05*nrow(d)))
cat(sprintf("\ngenes tested (mutated in at least 5%% of patients): %d\n", length(GENES)))
mut <- unique(maf[,c("pid","Hugo_Symbol")])

fisher_by <- function(sub_d, tag) {
  if (nrow(sub_d) < 20 || length(unique(sub_d$grp)) < 2) return(NULL)
  r <- do.call(rbind, lapply(GENES, function(g) {
    has <- sub_d$pid %in% mut$pid[mut$Hugo_Symbol==g]
    tb <- table(factor(sub_d$grp, c("Low","High")), factor(has, c(FALSE,TRUE)))
    if (any(dim(tb)!=2) || sum(tb[,2])<3) return(NULL)
    ft <- fisher.test(tb)
    data.frame(stratum=tag, gene=g, n=nrow(sub_d),
               pct_low=round(100*tb[1,2]/sum(tb[1,]),1), pct_high=round(100*tb[2,2]/sum(tb[2,]),1),
               OR=unname(ft$estimate), P=ft$p.value) }))
  if (is.null(r)) return(NULL)
  r$BH_q <- p.adjust(r$P,"BH"); r[order(r$BH_q),] }

cat("\n=== A. whole cohort (unstratified, for reference only; POLE confounding is not controlled) ===\n")
ra <- fisher_by(d, "all")
print(head(data.frame(gene=ra$gene, pct_low=ra$pct_low, pct_high=ra$pct_high,
                      OR=round(ra$OR,2), q=signif(ra$BH_q,3)), 12), row.names=FALSE)
cat(sprintf("BH q<0.05: %d / %d\n", sum(ra$BH_q<0.05), nrow(ra)))

cat("\n=== B. stratified by molecular subtype (required by section 22.1) ===\n")
rs <- do.call(rbind, lapply(c("CN low","CN high","MSI","POLE"),
       function(st) fisher_by(d[which(d$subtype==st),], st)))
if (!is.null(rs)) {
  rs$BH_q <- ave(rs$P, rs$stratum, FUN=function(p) p.adjust(p,"BH"))
  for (st in unique(rs$stratum)) {
    x <- rs[rs$stratum==st,]; x <- x[order(x$BH_q),]
    cat(sprintf("\n-- %s (n=%d): genes at BH q<0.05: %d / %d\n", st, x$n[1], sum(x$BH_q<0.05), nrow(x)))
    print(head(data.frame(gene=x$gene, pct_low=x$pct_low, pct_high=x$pct_high,
                          OR=round(x$OR,2), q=signif(x$BH_q,3)), 5), row.names=FALSE) }
  write.table(rs, "results/mutation/fisher_by_subtype.tsv", sep="\t", row.names=FALSE, quote=FALSE)
}
write.table(ra, "results/mutation/fisher_all.tsv", sep="\t", row.names=FALSE, quote=FALSE)

## ---- the directed hypothesis raised in section 20: ARID1A / SWI-SNF ----
cat("\n=== C. directed test: SWI/SNF members (motivated by the SMARCC2 co-expression of 0.610 in section 20) ===\n")
SWI <- c("ARID1A","ARID1B","ARID2","SMARCA4","SMARCB1","SMARCC1","SMARCC2","PBRM1")
for (g in SWI) {
  has <- d$pid %in% mut$pid[mut$Hugo_Symbol==g]
  if (sum(has) < 5) { cat(sprintf("  %-9s only %d mutated patients, too few, skipping\n", g, sum(has))); next }
  w <- wilcox.test(d$gx[has], d$gx[!has])
  cat(sprintf("  %-9s mutated %3d / %3d (%.1f%%) | median NUMA1, mutated %.3f vs wild-type %.3f | Wilcoxon P = %.3g\n",
      g, sum(has), nrow(d), 100*mean(has), median(d$gx[has]), median(d$gx[!has]), w$p.value))
}

saveRDS(list(d=d, ra=ra, rs=if(exists("rs")) rs else NULL), "results/mutation/mutation_objects.rds")
writeLines(capture.output(sessionInfo()), "logs/26_mutation_sessionInfo.txt")
cat("\n===== 26_mutation DONE =====\n")
