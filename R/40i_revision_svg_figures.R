#!/usr/bin/env Rscript
# 40i_revision_svg_figures.R -- editable SVG versions of Figs. S7 and S8 (revision) for redrawing in PowerPoint
#   Fig. S7  mean-difference (MA) plot of TCGA-UCEC under TMM: point density drawn as vector tiles, lowess trend, zero line
#   Fig. S8  scaled Schoenfeld residuals for the five univariable Cox models
# Text is written as SVG text (svglite, fix_text_size = FALSE) so that labels stay editable after import.
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(SummarizedExperiment); library(edgeR); library(survival); library(svglite)})
OUT <- "results/revision"; MM <- 1 / 25.4
FONT <- "Arial"

## ---------------------------------------------------------------- Fig. S7
se  <- readRDS("data_raw/TCGA/TCGA_UCEC_STARcounts_SE.rds")
cdl <- as.data.frame(colData(se)[, !sapply(as.data.frame(colData(se)), is.list), drop = FALSE]); code <- substr(cdl$barcode, 14, 15)
EN <- c("Endometrioid adenocarcinoma, NOS","Endometrioid adenocarcinoma, secretory variant")
SR <- c("Serous cystadenocarcinoma, NOS","Serous surface papillary carcinoma","Papillary serous cystadenocarcinoma")
keep <- (code == "01" & cdl$primary_diagnosis %in% c(EN, SR)) | code == "11"
cnt <- assay(se, "unstranded")[, keep]; rownames(cnt) <- rowData(se)$gene_name
o <- order(rownames(cnt), -rowMeans(cnt)); i <- o[!duplicated(rownames(cnt)[o])]; cnt <- cnt[i, ]; cnt <- cnt[!is.na(rownames(cnt)) & rownames(cnt) != "", ]
grp <- factor(ifelse(code[keep] == "11", "N", "T")); isT <- grp == "T"
d <- DGEList(cnt); d <- d[filterByExpr(d, group = grp), , keep.lib.sizes = FALSE]; d <- calcNormFactors(d)
lc <- cpm(d, log = TRUE)
A <- rowMeans(lc); M <- rowMeans(lc[, isT]) - rowMeans(lc[, !isT])
cat(sprintf("genes %d | median M %.3f | lowess at A quartiles: %s\n", length(A), median(M),
            paste(round(approx(lowess(A, M, f = 0.3), xout = quantile(A, c(.25, .5, .75)))$y, 3), collapse = ", ")))

# Panel a: distribution of tumor-minus-normal differences across expression levels, summarized per expression bin
# (median, IQR and 5-95% range, loess-smoothed) -- conveys centering without plotting 18,000 points.
# Panel b: normalization factors of tumors and normals under TMM, upper quartile and RLE.
o <- order(A); bins <- split(o, cut(seq_along(o), 60, labels = FALSE))
st <- t(sapply(bins, function(ix) c(A = median(A[ix]), quantile(M[ix], c(.05, .25, .5, .75, .95)))))
colnames(st) <- c("A", "q05", "q25", "q50", "q75", "q95")
xs <- seq(min(st[, "A"]), max(st[, "A"]), length.out = 200)
sm <- sapply(c("q05", "q25", "q50", "q75", "q95"), function(k) predict(loess(st[, k] ~ st[, "A"], span = 0.3), xs))
nf <- sapply(c(TMM = "TMM", `Upper quartile` = "upperquartile", RLE = "RLE"), function(m) calcNormFactors(d, method = m)$samples$norm.factors)
dA <- density(A, from = min(xs), to = max(xs), n = 200)

INK <- "#333333"; GRID <- "#E6E6E6"; BLUE <- "#1F4E8C"
band_out <- adjustcolor(BLUE, alpha.f = 0.13); band_in <- adjustcolor(BLUE, alpha.f = 0.32)
C_N <- "#C9C9C9"; C_T <- "#5B86BE"
xl <- range(xs); yl <- c(-3, 2.5)

svglite(file.path(OUT, "FigS7_MA_TMM.svg"), width = 174 * MM, height = 78 * MM, fix_text_size = FALSE)
layout(matrix(c(1, 2, 3, 3), 2, 2), heights = c(0.2, 1), widths = c(1.7, 1))
par(family = FONT, ps = 10, col.axis = INK, col.lab = INK, fg = INK, las = 1, lwd = 0.6, cex.axis = 0.9, tcl = -0.22, mgp = c(2.0, 0.45, 0))

## a, top strip: where the genes are along the expression axis
par(mar = c(0.2, 3.6, 1.5, 4.2))
plot(NA, xlim = xl, ylim = c(0, max(dA$y) * 1.08), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
polygon(c(dA$x, rev(dA$x)), c(dA$y, rep(0, length(dA$y))), col = "#D9D9D9", border = NA)
segments(xl[1], 0, xl[2], 0, col = "#BDBDBD", lwd = 0.5)
text(xl[1], max(dA$y) * 0.95, "Genes", adj = c(0, 1), cex = 0.8, col = "#7A7A7A")
mtext("a", side = 3, line = 0.15, adj = -0.14, font = 2, cex = 1.1, col = "black")

## a, main panel
par(mar = c(3.3, 3.6, 0.4, 4.2))
plot(NA, xlim = xl, ylim = yl, axes = FALSE, xaxs = "i", yaxs = "i",
     xlab = expression("Mean " * log[2] * " CPM (TMM)"), ylab = expression("Tumor " * "\u2212" * " normal (" * log[2] * ")"))
abline(h = setdiff(-3:2, 0), col = GRID, lwd = 0.5)
polygon(c(xs, rev(xs)), c(sm[, "q05"], rev(sm[, "q95"])), col = band_out, border = NA)
polygon(c(xs, rev(xs)), c(sm[, "q25"], rev(sm[, "q75"])), col = band_in, border = NA)
abline(h = 0, col = "#4D4D4D", lty = "22", lwd = 0.6)
lines(xs, sm[, "q50"], col = "white", lwd = 3.6)
lines(xs, sm[, "q50"], col = BLUE, lwd = 1.8)
mlab <- function(v) sub("^-", "\u2212", format(v, trim = TRUE))
xt <- axTicks(1); axis(1, at = xt, labels = mlab(xt), lwd = 0, lwd.ticks = 0.6)
axis(2, at = -3:2, labels = mlab(-3:2), lwd = 0, lwd.ticks = 0.6)
segments(xl[1], yl[1], xl[2], yl[1], lwd = 0.6, xpd = NA); segments(xl[1], yl[1], xl[1], yl[2], lwd = 0.6, xpd = NA)
# direct labels at the right end of each band boundary
ye <- sm[nrow(sm), ]; lab_y <- c(ye["q05"], ye["q25"], ye["q50"], ye["q75"], ye["q95"])
lab_y[2] <- min(lab_y[2], lab_y[3] - 0.32); lab_y[4] <- max(lab_y[4], lab_y[3] + 0.32)
text(xl[2] + 0.18, lab_y, c("5th", "25th", "Median", "75th", "95th"), adj = c(0, 0.5), xpd = NA, cex = 0.8,
     col = c("#7A7A7A", "#7A7A7A", BLUE, "#7A7A7A", "#7A7A7A"), font = c(1, 1, 2, 1, 1))
text(xl[2] + 0.18, yl[2] - 0.05, "Percentile", adj = c(0, 1), xpd = NA, cex = 0.8, col = "#7A7A7A", font = 3)

## b
par(mar = c(3.3, 3.4, 1.5, 0.6))
pos <- c(1, 2, 4, 5, 7, 8)
vals <- list(nf[!isT, 1], nf[isT, 1], nf[!isT, 2], nf[isT, 2], nf[!isT, 3], nf[isT, 3])
plot(NA, xlim = c(0.3, 8.7), ylim = c(0.6, 1.6), axes = FALSE, xaxs = "i", yaxs = "i", xlab = "", ylab = "Normalization factor")
abline(h = c(0.8, 1.2, 1.4), col = GRID, lwd = 0.5); abline(h = 1, col = "#4D4D4D", lty = "22", lwd = 0.6)
boxplot(vals, at = pos, col = rep(c(C_N, C_T), 3), border = "#2B2B2B", outline = FALSE, axes = FALSE, add = TRUE,
        boxwex = 0.72, staplewex = 0, lty = 1, lwd = 0.6, whisklwd = 0.6, medlwd = 1.4, medcol = "#1A1A1A")
axis(2, at = seq(0.6, 1.6, 0.2), lwd = 0, lwd.ticks = 0.6)
segments(0.3, 0.6, 0.3, 1.6, lwd = 0.6, xpd = NA); segments(0.3, 0.6, 8.7, 0.6, lwd = 0.6, xpd = NA)
axis(1, at = c(1.5, 4.5, 7.5), labels = c("TMM", "UQ", "RLE"), lwd = 0, line = -0.5)
rect(5.9, c(1.52, 1.43) - 0.03, 6.25, c(1.52, 1.43) + 0.03, col = c(C_N, C_T), border = "#2B2B2B", lwd = 0.6)
text(6.45, c(1.52, 1.43), c("Normal", "Tumor"), adj = c(0, 0.5), cex = 0.9)
mtext("b", side = 3, line = 0.15, adj = -0.32, font = 2, cex = 1.1, col = "black")
dev.off()

## ---------------------------------------------------------------- Fig. S8
# one record per patient (500 patients, 78 deaths; R/40j)
sv <- read.delim(file.path(OUT, "survival_per_patient_deduplicated.tsv"), check.names = FALSE)
vars <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1","panel")
ph <- read.delim(file.path(OUT, "survival_cox_PH.tsv"), check.names = FALSE)
stopifnot(nrow(sv) == 500, all(ph$events == sum(sv$OS_event)))
draw_s8 <- function() {
  par(family = FONT, mfrow = c(2, 3), mar = c(3.0, 3.4, 1.8, 0.6), mgp = c(1.9, 0.5, 0), tcl = -0.3, cex = 0.7, las = 1)
  for (v in vars) {
    sv$x <- as.numeric(scale(sv[[v]]))
    z <- cox.zph(coxph(Surv(OS_time, OS_event) ~ x, data = sv), transform = "km")
    plot(z, resid = TRUE, se = TRUE, col = c("#1F5A99", "#1F5A99"), lwd = 1.4, pch = 16, cex = 0.45,
         xlab = "Time (Kaplan–Meier scale)", ylab = "Beta(t), per SD", main = "")
    abline(h = 0, lty = 1, col = "grey40")
    title(main = if (v == "panel") "Panel score" else v, font.main = if (v == "panel") 2 else 4, cex.main = 1.05, line = 0.5, adj = 0)
    mtext(sprintf("PH test P = %.3f", ph$PH_P[ph$variable == v]), side = 3, line = 0.45, adj = 1, cex = 0.66)
  }
  plot.new()
}
svglite(file.path(OUT, "FigS8_Schoenfeld.svg"), width = 174 * MM, height = 110 * MM, fix_text_size = FALSE); draw_s8(); dev.off()
ragg::agg_png(file.path(OUT, "FigS8_Schoenfeld.png"), width = 174 * MM, height = 110 * MM, units = "in", res = 600, background = "white"); draw_s8(); dev.off()
cat("wrote FigS7_MA_TMM.svg and FigS8_Schoenfeld.svg\n")
