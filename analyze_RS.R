### =========================================================
### analyze_RS.R — Figures and statistics from 100-rep RS run
### =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

data.table::setDTthreads(1)  

t0 <- Sys.time()

rds_path <- "results/stats_RS_v2_full100.rds"
if (!file.exists(rds_path)) {
  stop("Results file not found: ", rds_path)
}

stats_all <- readRDS(rds_path)
setDT(stats_all)

nreps <- length(unique(stats_all$rep))
ncyc  <- max(stats_all$cycle)
cat(sprintf("Loaded %d reps x %d cycles x 2 encodings (%d rows)\n",
            nreps, ncyc, nrow(stats_all)))

theme_pub <- function() {
  theme_classic(base_size = 13) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      plot.margin = margin(8, 8, 8, 8, "pt"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.4),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 14)
    )
}

cols <- c("C" = "#D73027", "Bmid" = "#1A9850")
labs_calling <- c("C" = "Triploid encoding", "Bmid" = "Diploidised Encoding")

stats_all[, cumH := mean_hyb_TBV - first(mean_hyb_TBV), by = .(rep, calling)]

summ_cycle <- stats_all[, .(
  meanH  = mean(mean_hyb_TBV),
  seH    = sd(mean_hyb_TBV) / sqrt(.N),
  meanCumH = mean(cumH),
  seCumH   = sd(cumH) / sqrt(.N),
  N = .N
), by = .(cycle, calling)]

cat("\n=== Per-cycle means across reps ===\n")
print(summ_cycle)

## Figure 1: Hybrid genetic value trajectory
p1 <- ggplot(summ_cycle, aes(x = cycle, y = meanH, color = calling, fill = calling)) +
  geom_line(linewidth = 1.3) +
  geom_ribbon(aes(ymin = meanH - seH, ymax = meanH + seH),
              alpha = 0.18, color = NA) +
  geom_point(size = 2.4) +
  scale_color_manual(values = cols, labels = labs_calling) +
  scale_fill_manual(values = cols, labels = labs_calling) +
  scale_x_continuous(breaks = 1:ncyc) +
  labs(x = "Selection cycle",
       y = "Mean hybrid genetic value",
       color = "Genotype encoding",
       fill  = "Genotype encoding",
       title = sprintf("Hybrid genetic gain across cycles", nreps)) +
  theme_pub()
ggsave("results/fig_RS_gain_trajectory.png", p1, width = 8, height = 5.5, dpi = 300)
cat("  saved fig_RS_gain_trajectory.png\n")

## Figure 2: Cumulative gain
p2 <- ggplot(summ_cycle, aes(x = cycle, y = meanCumH, color = calling, fill = calling)) +
  geom_line(linewidth = 1.3) +
  geom_ribbon(aes(ymin = meanCumH - seCumH, ymax = meanCumH + seCumH),
              alpha = 0.18, color = NA) +
  geom_point(size = 2.4) +
  scale_color_manual(values = cols, labels = labs_calling) +
  scale_fill_manual(values = cols, labels = labs_calling) +
  scale_x_continuous(breaks = 1:ncyc) +
  labs(x = "Selection cycle",
       y = "Cumulative gain relative to cycle 1",
       color = "Genotype encoding",
       fill  = "Genotype encoding",
       title = sprintf("Cumulative hybrid gain", nreps)) +
  theme_pub()
ggsave("results/fig_RS_cumulative_gain.png", p2, width = 8, height = 5.5, dpi = 300)
cat("  saved fig_RS_cumulative_gain.png\n")

## Figure 3: Final-cycle boxplot
final_stats <- stats_all[cycle == ncyc]

p3 <- ggplot(final_stats, aes(x = calling, y = mean_hyb_TBV, fill = calling)) +
  geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = 21) +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = cols, labels = labs_calling) +
  scale_x_discrete(labels = labs_calling) +
  labs(x = NULL,
       y = sprintf("Hybrid genetic value at cycle %d", ncyc),
       fill = NULL,
       title = sprintf("Final-cycle hybrid genetic value distribution (%d reps)", nreps)) +
  theme_pub() +
  theme(legend.position = "none")
ggsave("results/fig_RS_final_boxplot.png", p3, width = 6, height = 5.5, dpi = 300)
cat("  saved fig_RS_final_boxplot.png\n")

## Figure 4: Paired scatter
end_wide <- dcast(final_stats, rep ~ calling, value.var = "mean_hyb_TBV")

p4 <- ggplot(end_wide, aes(x = Bmid, y = C)) +
  geom_point(size = 2.5, alpha = 0.7, color = "#333333") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Diploidised kernel hybrid genetic value",
       y = "Triploid Kernel hybrid genetic value",
       title = sprintf("Paired final genetic value triploid vs diploid kernel", nreps)) +
  theme_pub()
ggsave("results/fig_RS_paired_scatter.png", p4, width = 6.5, height = 6, dpi = 300)
cat("  saved fig_RS_paired_scatter.png\n")

## Headline stats (unchanged)
end_wide[, diff := C - Bmid]
cat("\n=== HEADLINE STATISTICS ===\n")
cat(sprintf("N replicates: %d\n", nrow(end_wide)))
cat(sprintf("Mean C final: %.4f (SD %.4f)\n", mean(end_wide$C), sd(end_wide$C)))
cat(sprintf("Mean B final: %.4f (SD %.4f)\n", mean(end_wide$Bmid), sd(end_wide$Bmid)))
cat(sprintf("Mean diff (C - B): %.4f\n", mean(end_wide$diff)))
cat(sprintf("SE of diff: %.4f\n", sd(end_wide$diff) / sqrt(nrow(end_wide))))
cat(sprintf("95%% CI: [%.4f, %.4f]\n",
            mean(end_wide$diff) - 1.96 * sd(end_wide$diff) / sqrt(nrow(end_wide)),
            mean(end_wide$diff) + 1.96 * sd(end_wide$diff) / sqrt(nrow(end_wide))))
cat(sprintf("Relative gain: %.2f%%\n", 100 * mean(end_wide$diff) / mean(end_wide$Bmid)))

tt <- t.test(end_wide$C, end_wide$Bmid, paired = TRUE)
cat(sprintf("Paired t-test: t = %.3f, df = %d, p = %.4g\n",
            tt$statistic, tt$parameter, tt$p.value))
wt <- wilcox.test(end_wide$C, end_wide$Bmid, paired = TRUE)
cat(sprintf("Wilcoxon: V = %.0f, p = %.4g\n", wt$statistic, wt$p.value))
cat(sprintf("Reps where C > B: %d / %d (%.1f%%)\n",
            sum(end_wide$diff > 0), nrow(end_wide), 100 * mean(end_wide$diff > 0)))

fwrite(summ_cycle,  "results/RS_summary_by_cycle.csv")
fwrite(end_wide,    "results/RS_final_paired.csv")

cat(sprintf("\nTotal time: %.1f sec\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat("=== ANALYSIS DONE ===\n")
EOF

getwd()

