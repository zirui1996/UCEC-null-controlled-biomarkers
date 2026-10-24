#!/usr/bin/env Rscript
# 40n_revision_FigS3_patient_level.R -- Fig. S3 redrawn with one record per patient (500 patients, 78 deaths; R/40j).
# Same layout as R/supp/FS3_survival.R. Median-split curves are for visualization only, so the log-rank P values of the
# original figure are no longer printed; panel labels give the Cox estimates that are the basis of inference.
setwd(Sys.getenv("UCEC_PROJ", unset = getwd()))
suppressMessages({library(survival); library(ggplot2); library(patchwork); library(dplyr)})
# figure theme (same settings as the other figures of the paper)
COL <- c(core = "#E28E2C", mid = "#767676"); CORE <- c("HAND2-AS1","OSR2","PEG3","RUNX1T1"); BASE <- 6.3
theme_set(theme_classic(base_size = BASE, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.3, colour = "black"), axis.ticks = element_line(linewidth = 0.3, colour = "black"),
        axis.ticks.length = unit(1.2, "pt"), axis.text = element_text(size = BASE - 0.5, colour = "black"),
        axis.title = element_text(size = BASE), legend.title = element_text(size = BASE - 0.3),
        legend.text = element_text(size = BASE - 0.7), legend.key.size = unit(6, "pt"), legend.background = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = BASE, face = "bold"),
        plot.title = element_text(size = BASE + 0.3, face = "bold", hjust = 0), plot.title.position = "plot",
        plot.subtitle = element_text(size = BASE - 0.5, colour = COL[["mid"]]), plot.margin = margin(3, 3, 3, 3),
        panel.grid = element_blank()))
update_geom_defaults("text", list(family = "Arial", size = (BASE - 0.7) / .pt))
tag_theme <- theme(plot.tag = element_text(size = 8, face = "bold", family = "Arial"))
sv <- read.delim("results/revision/survival_per_patient_deduplicated.tsv", check.names = FALSE, stringsAsFactors = FALSE)
cx <- read.delim("results/revision/survival_cox_PH.tsv", stringsAsFactors = FALSE)
stopifnot(nrow(sv) == 500, sum(sv$OS_event) == 78, all(cx$events == 78))
cx <- transmute(cx, gene = variable, HR = HR_per_SD, lo = CI_low, hi = CI_high, P, BH_q)
FEAT <- c(CORE, "panel"); LABF <- c(setNames(CORE, CORE), panel = "Panel score")
km <- bind_rows(lapply(FEAT, function(f) {
  sv$hg <- factor(ifelse(sv[[f]] > median(sv[[f]]), "High", "Low"), levels = c("Low", "High"))
  fit <- survfit(Surv(OS_time / 365.25, OS_event) ~ hg, data = sv)
  s <- summary(fit, censored = TRUE)
  data.frame(feature = f, time = c(0, 0, s$time), surv = c(1, 1, s$surv),
             strata = c("hg=Low", "hg=High", as.character(s$strata)), n.censor = c(0, 0, s$n.censor)) %>%
    mutate(group = sub("hg=", "", strata))
}))
km <- km %>% arrange(feature, group, time)
lab <- cx %>% transmute(feature = gene, txt = sprintf("HR per SD %.2f (%.2f–%.2f)\nCox P = %.3f", HR, lo, hi, P))
km$feature <- factor(km$feature, levels = FEAT, labels = LABF[FEAT]); lab$feature <- factor(lab$feature, levels = FEAT, labels = LABF[FEAT])
a <- ggplot(km, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.45) +
  geom_point(data = filter(km, n.censor > 0), shape = 3, size = 0.6, stroke = 0.3) +
  geom_text(data = lab, aes(x = 0.2, y = 0.08, label = txt), inherit.aes = FALSE, hjust = 0, vjust = 0, size = 4.8/.pt, lineheight = 0.9) +
  facet_wrap(~feature, nrow = 2) +
  scale_colour_manual(values = c(Low = "#5B6C8F", High = COL[["core"]]), name = "Split at median\n(visualization only)") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) +
  labs(x = "Years", y = "Overall survival", title = "Kaplan–Meier curves (median split)",
       subtitle = sprintf("TCGA-UCEC, %d patients, %d deaths; + censored", nrow(sv), sum(sv$OS_event))) +
  theme(strip.text = element_text(face = "bold.italic"), legend.position = c(0.84, 0.22), legend.justification = c(0.5, 0.5))
fp <- cx %>% mutate(label = factor(LABF[gene], levels = rev(LABF[FEAT])),
                    txt = sprintf("%.2f (%.2f–%.2f)   q = %.2f", HR, lo, hi, BH_q))
b <- ggplot(fp, aes(HR, label)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, linewidth = 0.35, orientation = "y") +
  geom_point(shape = 22, fill = COL[["core"]], size = 2, stroke = 0.3) +
  geom_text(aes(x = 1.75, label = txt), hjust = 0, size = 5/.pt) +
  scale_x_log10(limits = c(0.6, 3.3), breaks = c(0.6, 0.8, 1, 1.25, 1.5)) +
  labs(x = "Hazard ratio per SD (log scale, 95% CI)", y = NULL, title = "Univariable Cox regression",
       subtitle = "BH q across five models") +
  theme(axis.text.y = element_text(face = "italic"))
fig <- a / b + plot_layout(heights = c(1, 0.42)) + plot_annotation(tag_levels = "a") & tag_theme
f <- "results/revision/FigS3_patient_level"; w <- 183 / 25.4; h <- 150 / 25.4
svglite::svglite(paste0(f, ".svg"), width = w, height = h, system_fonts = list(sans = "Arial"), fix_text_size = FALSE); print(fig); dev.off()
grDevices::quartz(type = "pdf", file = paste0(f, ".pdf"), width = w, height = h, family = "Arial"); print(fig); dev.off()
ragg::agg_png(paste0(f, ".png"), width = w, height = h, units = "in", res = 600, background = "white"); print(fig); dev.off()
cat("Fig. S3 (patient level) done\n")
