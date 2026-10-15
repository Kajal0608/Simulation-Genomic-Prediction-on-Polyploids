### =========================================================
### 07_figures_50reps.R — Regenerate all sensitivity figures
### from the 50-rep sensitivity analysis
### =========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

data.table::setDTthreads(1)

## Load 50-rep results (adjust filename)
res <- fread("results_dosage_resolution_sensitivity_50reps.csv")

## ---------------------------------------------------------
## Aggregate to scenario level (mean and SE across 50 reps)
## ---------------------------------------------------------

## Add scenario id
res[, scenario := paste(meanDD, varDD, h2_base, n_cross, n_prog, sep="_")]

## Compute realized dominance proportion per rep
res[, prop_D := Vd_trip / (Va_trip + Vd_trip)]

## Aggregate to scenario means with SE
agg <- res[, .(
  # accuracies
  mean_C_AD = mean(C_AD, na.rm = TRUE),
  se_C_AD   = sd(C_AD, na.rm = TRUE) / sqrt(.N),
  mean_B_AD = mean(Bmid_AD, na.rm = TRUE),
  se_B_AD   = sd(Bmid_AD, na.rm = TRUE) / sqrt(.N),
  
  # gaps
  mean_gap_AD = mean(gap_AD, na.rm = TRUE),
  se_gap_AD   = sd(gap_AD, na.rm = TRUE) / sqrt(.N),
  
  # dominance proportion
  mean_prop_D = mean(prop_D, na.rm = TRUE),
  
  # factors
  meanDD, varDD, h2_base, n_cross, n_prog,
  n_train = n_cross * n_prog,
  n_reps  = .N
), by = scenario]
agg <- unique(agg, by = "scenario")


theme_pub <- function() {
  theme_classic(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 13),
      plot.margin = margin(6, 6, 6, 6, "pt")
    )
}

cols <- c("C" = "#D73027", "B" = "#1A9850")
lab_C <- "C (proper triploid dosage)"
lab_B <- "B (diploidized encoding)"


fig2 <- ggplot(agg, aes(x = mean_B_AD, y = mean_C_AD)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = mean_C_AD - se_C_AD,
                    ymax = mean_C_AD + se_C_AD),
                color = "gray70", width = 0, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = mean_B_AD - se_B_AD,
                     xmax = mean_B_AD + se_B_AD),
                 color = "gray70", height = 0, linewidth = 0.3) +
  geom_point(aes(color = factor(meanDD)), size = 3, alpha = 0.85) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(
    values = c("0.2" = "#D73027", "0.6" = "#4575B4", "1" = "#1A9850"),
    name = "Mean dominance degree"
  ) +
  labs(
    x = "Diplodised encoding A+D accuracy",
    y = "Triploid encoding A+D accuracy",
    title = "Triploid vs diploidised encoding across different scenarios"
  ) +
  theme_pub() +
  coord_fixed()

ggsave("results/fig2_paired_scatter_50reps.png", fig2,
       width = 8, height = 8, dpi = 300, bg = "white")


marg_by <- function(factor_col) {
  x <- agg[, .(
    C = mean(mean_C_AD, na.rm = TRUE),
    C_se = sd(mean_C_AD, na.rm = TRUE) / sqrt(.N),
    B = mean(mean_B_AD, na.rm = TRUE),
    B_se = sd(mean_B_AD, na.rm = TRUE) / sqrt(.N)
  ), by = factor_col]
  
  x_long_C <- x[, .(x = get(factor_col), encoding = "C",
                    accuracy = C, se = C_se)]
  x_long_B <- x[, .(x = get(factor_col), encoding = "B",
                    accuracy = B, se = B_se)]
  rbind(x_long_C, x_long_B)
}

## Panel A: training size
P_A <- marg_by("n_train")
pA <- ggplot(P_A, aes(x = factor(x), y = accuracy,
                      color = encoding, group = encoding)) +
  geom_ribbon(aes(ymin = accuracy - se, ymax = accuracy + se,
                  fill = encoding, group = encoding),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = cols, labels = c("C" = lab_C, "B" = lab_B)) +
  scale_fill_manual(values = cols, labels = c("C" = lab_C, "B" = lab_B)) +
  labs(x = "Training population size", y = "Prediction accuracy",
       color = NULL, fill = NULL, title = "A. Training population size") +
  theme_pub()

## Panel B: meanDD
P_B <- marg_by("meanDD")
pB <- ggplot(P_B, aes(x = x, y = accuracy,
                      color = encoding, group = encoding)) +
  geom_ribbon(aes(ymin = accuracy - se, ymax = accuracy + se,
                  fill = encoding, group = encoding),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = sort(unique(agg$meanDD))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = cols, labels = c("C" = lab_C, "B" = lab_B)) +
  scale_fill_manual(values = cols, labels = c("C" = lab_C, "B" = lab_B)) +
  labs(x = "Mean dominance degree", y = "Prediction accuracy",
       color = NULL, fill = NULL, title = "B. Mean dominance degree") +
  theme_pub()

## Panel C: h2_base
P_C <- marg_by("h2_base")
pC <- ggplot(P_C, aes(x = x, y = accuracy,
                      color = encoding, group = encoding)) +
  geom_ribbon(aes(ymin = accuracy - se, ymax = accuracy + se,
                  fill = encoding, group = encoding),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = sort(unique(agg$h2_base))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = cols, labels = c("C" = lab_C, "B" = lab_B)) +
  scale_fill_manual(values = cols, labels = c("C" = lab_C, "B" = lab_B)) +
  labs(x = "Baseline heritability", y = "Prediction accuracy",
       color = NULL, fill = NULL, title = "C. Baseline heritability") +
  theme_pub()


fig3 <- (pA | pB | pC) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

ggsave("results/fig3_composite_50reps.png", fig3,
       width = 14, height = 5, dpi = 300, bg = "white")


fig4 <- ggplot(agg, aes(x = mean_prop_D, y = mean_gap_AD)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = mean_gap_AD - se_gap_AD,
                    ymax = mean_gap_AD + se_gap_AD),
                color = "gray70", width = 0, linewidth = 0.3) +
  geom_point(aes(color = factor(h2_base)), size = 3, alpha = 0.85) +
  geom_smooth(method = "loess", color = "black", se = FALSE, linewidth = 0.9) +
  scale_color_manual(
    values = c("0.05" = "#D73027", "0.12" = "#4575B4", "0.3" = "#1A9850"),
    name = "Baseline heritability"
  ) +
  labs(
    x = "Realized dominance proportion V_D / (V_A + V_D)",
    y = "Encoding advantage (C accuracy - B accuracy)",
    title = "Encoding advantage decreases with dominance proportion"
  ) +
  theme_pub()

ggsave("results/fig4_gap_vs_dominance_50reps.png", fig4,
       width = 8, height = 6, dpi = 300, bg = "white")


cat("\n===== 50-REP SENSITIVITY SUMMARY =====\n")
cat("Total scenarios:", nrow(agg), "\n")
cat("Reps per scenario:", agg$n_reps[1], "\n\n")

cat("Under encoding C (AD model):\n")
cat("  Range:", round(min(agg$mean_C_AD), 3), "to", round(max(agg$mean_C_AD), 3), "\n")
cat("  Mean:", round(mean(agg$mean_C_AD), 3), "\n\n")

cat("Under encoding B (AD model):\n")
cat("  Range:", round(min(agg$mean_B_AD), 3), "to", round(max(agg$mean_B_AD), 3), "\n")
cat("  Mean:", round(mean(agg$mean_B_AD), 3), "\n\n")

cat("Encoding advantage (gap_AD):\n")
cat("  Mean:", round(mean(agg$mean_gap_AD), 3), "\n")
cat("  SD:", round(sd(agg$mean_gap_AD), 3), "\n")
cat("  Range:", round(min(agg$mean_gap_AD), 3), "to", round(max(agg$mean_gap_AD), 3), "\n")
cat("  Scenarios favoring C:", sum(agg$mean_gap_AD > 0), "of", nrow(agg), "\n")

cat("\nSaved figures to results/\n")

cat("\n===== A-MODEL SUMMARY =====\n")
cat("Encoding advantage (gap_A):\n")
cat("  Mean:", round(mean(agg$mean_gap_A, na.rm = TRUE), 3), "\n")
cat("  SD:", round(sd(agg$mean_gap_A, na.rm = TRUE), 3), "\n")
cat("  Range:", round(min(agg$mean_gap_A, na.rm = TRUE), 3),
    "to", round(max(agg$mean_gap_A, na.rm = TRUE), 3), "\n")
cat("  Scenarios favoring C (A model):", 
    sum(agg$mean_gap_A > 0), "of", nrow(agg), "\n")

