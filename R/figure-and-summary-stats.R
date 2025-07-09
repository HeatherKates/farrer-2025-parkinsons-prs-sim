library(readxl)
library(ggplot2)

# --- 1. Load excel file ---
# Adjust 'mysheet' as needed; sheet = 1 for first sheet, or character string for sheet name.
file <- "../data/Leonard_2025_file04_S3.xlsx"
mysheet <- 1
d0 <- as.data.frame(read_excel(file, sheet = mysheet))

# --- 2. Remove rs34637584 ---
d0 <- d0[as.character(d0$rsID) != "rs34637584", ]

# --- 3. Use default values for Beta and AF ---
beta_col <- grep("Beta, all studies", names(d0), ignore.case = TRUE, value = TRUE)[1]
af_col   <- grep("Effect allele frequency", names(d0), ignore.case = TRUE, value = TRUE)[1]
rsid_col <- grep("rsid", names(d0), ignore.case = TRUE, value = TRUE)[1]
pop_size <- 100000

# --- Simulation function (returns list with PRS, OR, stats) ---
simulate_prs <- function(table, rsid_col, beta_col, af_col, pop_size) {
  table <- table[complete.cases(table[, c(rsid_col, beta_col, af_col)]), ]
  AF <- as.numeric(table[[af_col]])
  Beta <- as.numeric(table[[beta_col]])
  keep <- !is.na(AF) & !is.na(Beta) & AF >= 0 & AF <= 1
  AF <- AF[keep]; Beta <- Beta[keep]; rsids <- table[[rsid_col]][keep]
  genos <- sapply(AF, function(freq) rbinom(pop_size, 2, freq))
  colnames(genos) <- rsids
  PRS <- as.vector(genos %*% Beta)
  OR <- exp(PRS)
  list(
    PRS = PRS,
    OR = OR,
    mean_PRS = mean(PRS),
    sd_PRS = sd(PRS),
    min_PRS = min(PRS),
    max_PRS = max(PRS),
    mean_OR = mean(OR),
    sd_OR = sd(OR),
    min_OR = min(OR),
    max_OR = max(OR)
  )
}

# --- 4. Simulation 1: "Known" only ---
main_cols <- c(rsid_col, beta_col, af_col)
d0 <- d0[complete.cases(d0[, main_cols]), , drop = FALSE]
if ("Known or Novel" %in% names(d0)) {
  d_known <- d0[d0[["Known or Novel"]] == "Known", ]
  d_known <- d_known[complete.cases(d_known[, main_cols]), , drop = FALSE]
} else {
  stop("Column 'Known or Novel' not found in file!")
}
sim1 <- simulate_prs(d_known, rsid_col, beta_col, af_col, pop_size)

# --- 5. Simulation 2: All data after filtering ---
sim2 <- simulate_prs(d0, rsid_col, beta_col, af_col, pop_size)

# --- 6. Plot both curves with marked means, SDs, threshold ---
df1 <- data.frame(PRS = sim1$PRS, Set = "Known")
df2 <- data.frame(PRS = sim2$PRS, Set = "All")
df <- rbind(df1, df2)

mean1 <- sim1$mean_PRS
mean2 <- sim2$mean_PRS
sd1 <- sim1$sd_PRS
sd2 <- sim2$sd_PRS
xcut <- 2.2

cols <- c(Known = "#3182bd", All = "#de2d26")
plt <- ggplot(df, aes(x = PRS, color = Set, fill = Set)) +
  geom_density(alpha = 0.25, adjust = 2) +
  geom_vline(aes(xintercept = mean1), color = cols["Known"], linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = mean1 + sd1), color = cols["Known"], linetype = 3, linewidth = 0.9) +
  geom_vline(aes(xintercept = mean1 + 2*sd1), color = cols["Known"], linetype = 3, linewidth = 0.7) +
  geom_vline(aes(xintercept = mean1 + 3*sd1), color = cols["Known"], linetype = 3, linewidth = 0.5) +
  geom_vline(aes(xintercept = mean1 - sd1), color = cols["Known"], linetype = 3, linewidth = 0.9) +
  geom_vline(aes(xintercept = mean1 - 2*sd1), color = cols["Known"], linetype = 3, linewidth = 0.7) +
  geom_vline(aes(xintercept = mean1 - 3*sd1), color = cols["Known"], linetype = 0.5) +
  geom_vline(aes(xintercept = mean2), color = cols["All"], linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = mean2 + sd2), color = cols["All"], linetype = 3, linewidth = 0.9) +
  geom_vline(aes(xintercept = mean2 + 2*sd2), color = cols["All"], linetype = 3, linewidth = 0.7) +
  geom_vline(aes(xintercept = mean2 + 3*sd2), color = cols["All"], linetype = 3, linewidth = 0.5) +
  geom_vline(aes(xintercept = mean2 - sd2), color = cols["All"], linetype = 3, linewidth = 0.9) +
  geom_vline(aes(xintercept = mean2 - 2*sd2), color = cols["All"], linetype = 3, linewidth = 0.7) +
  geom_vline(aes(xintercept = mean2 - 3*sd2), color = cols["All"], linetype = 0.5) +
  geom_vline(xintercept = xcut, color = "purple", linetype = 1, linewidth = 1.2) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(title = "Simulated PRS Distributions",
       subtitle = "Dashed = mean; dotted = +/-1,2,3 SD; purple = 2.2",
       x = "PRS", y = "Density") +
  theme_minimal(base_size = 16)

# --- 7. Summary stats for both ---
cstat <- function(x, mu, sd) {
  list(
    n1 = sum(x >= mu - sd & x <= mu + sd),
    n2 = sum(x >= mu - 2*sd & x <= mu + 2*sd),
    n3 = sum(x >= mu - 3*sd & x <= mu + 3*sd),
    n_cut = sum(x >= xcut),
    pct1 = mean(x >= mu - sd & x <= mu + sd),
    pct2 = mean(x >= mu - 2*sd & x <= mu + 2*sd),
    pct3 = mean(x >= mu - 3*sd & x <= mu + 3*sd),
    pct_cut = mean(x >= xcut)
  )
}
sum1 <- cstat(sim1$PRS, mean1, sd1)
sum2 <- cstat(sim2$PRS, mean2, sd2)

summary_text <- paste(
  "=== KNOWN ===\n",
  sprintf("Mean PRS: %.4f | SD: %.4f | Min: %.4f | Max: %.4f", mean1, sd1, sim1$min_PRS, sim1$max_PRS), "\n",
  sprintf("Mean OR: %.4f | SD: %.4f | Min: %.4f | Max: %.4f", sim1$mean_OR, sim1$sd_OR, sim1$min_OR, sim1$max_OR), "\n",
  sprintf("N within 1SD: %d (%.2f%%)", sum1$n1, 100*sum1$pct1), "\n",
  sprintf("N within 2SD: %d (%.2f%%)", sum1$n2, 100*sum1$pct2), "\n",
  sprintf("N within 3SD: %d (%.2f%%)", sum1$n3, 100*sum1$pct3), "\n",
  sprintf("N with PRS >= %.2f: %d (%.2f%%)", xcut, sum1$n_cut, 100*sum1$pct_cut), "\n\n",
  "=== ALL ===\n",
  sprintf("Mean PRS: %.4f | SD: %.4f | Min: %.4f | Max: %.4f", mean2, sd2, sim2$min_PRS, sim2$max_PRS), "\n",
  sprintf("Mean OR: %.4f | SD: %.4f | Min: %.4f | Max: %.4f", sim2$mean_OR, sim2$sd_OR, sim2$min_OR, sim2$max_OR), "\n",
  sprintf("N within 1SD: %d (%.2f%%)", sum2$n1, 100*sum2$pct1), "\n",
  sprintf("N within 2SD: %d (%.2f%%)", sum2$n2, 100*sum2$pct2), "\n",
  sprintf("N within 3SD: %d (%.2f%%)", sum2$n3, 100*sum2$pct3), "\n",
  sprintf("N with PRS >= %.2f: %d (%.2f%%)", xcut, sum2$n_cut, 100*sum2$pct_cut)
)

# --- 8. Export image (high quality) ---
ggsave("PRS_distributions.png", plt, width = 11, height = 7, dpi = 300)

# --- 9. Save stats as text ---
writeLines(summary_text, "PRS_simulation_summary.txt")

cat(summary_text)
cat("\nFigure saved as PRS_distributions.png\n")
cat("Summary saved as PRS_simulation_summary.txt\n")

## Standardized PRS

## Standardized PRS analysis and plot

# 1. Standardize PRS vectors to 'All' mean/SD
std_mean <- sim2$mean_PRS
std_sd   <- sim2$sd_PRS
xcut <- 2.2
xcut_std <- (xcut - std_mean) / std_sd

PRS1_std <- (sim1$PRS - std_mean) / std_sd
PRS2_std <- (sim2$PRS - std_mean) / std_sd

df_std <- rbind(
  data.frame(PRS = PRS1_std, Set = "Known (Standardized to All)"),
  data.frame(PRS = PRS2_std, Set = "All (Standardized to All)")
)

cols_std <- c("Known (Standardized to All)" = "#3182bd", "All (Standardized to All)" = "#de2d26")

# 2. Plot, with one set of reference mean/SD lines
plt_std <- ggplot(df_std, aes(x = PRS, color = Set, fill = Set)) +
  geom_density(alpha = 0.25, adjust = 2) +
  geom_vline(xintercept = 0, color = "gray30", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = c(-1, 1), color = "gray50", linetype = 3, linewidth = 0.9) +
  geom_vline(xintercept = c(-2, 2), color = "gray50", linetype = 3, linewidth = 0.7) +
  geom_vline(xintercept = c(-3, 3), color = "gray50", linetype = 3, linewidth = 0.5) +
  geom_vline(xintercept = xcut_std, color = "purple", linetype = 1, linewidth = 1.2) +
  scale_color_manual(values = cols_std) +
  scale_fill_manual(values = cols_std) +
  labs(
    title = "Standardized PRS Distributions",
    subtitle = "Standardized to All SNPs: Dashed = Mean (All), Gray Dotted = ±1,2,3 SD (All), Purple = 2.2 threshold",
    x = "Standardized PRS (All SNPs)", y = "Density"
  ) +
  theme_bw(base_size = 16)
ggsave("Standardized_PRS_distributions.png", plt_std, width = 11, height = 7, dpi = 300)
cat("Standardized PRS plot saved as Standardized_PRS_distributions.png\n")

# 3. Summary stats: use reference mean=0, SD=1, and threshold from All for both groups
cstat_std <- function(x, cut) {
  n <- length(x)
  list(
    mean = mean(x),
    sd = sd(x),
    min = min(x),
    max = max(x),
    n1 = sum(abs(x) < 1),    pct1 = mean(abs(x) < 1),
    n2 = sum(abs(x) < 2),    pct2 = mean(abs(x) < 2),
    n3 = sum(abs(x) < 3),    pct3 = mean(abs(x) < 3),
    n_cut = sum(x >= cut), pct_cut = mean(x >= cut)
  )
}
sum1s <- cstat_std(PRS1_std, xcut_std)
sum2s <- cstat_std(PRS2_std, xcut_std)
summary_text_std <- paste(
  "=== KNOWN (Standardized to All SNPs) ===\n",
  sprintf("Mean (standardized): %.4f | SD (standardized): %.4f | Min: %.4f | Max: %.4f",
          sum1s$mean, sum1s$sd, sum1s$min, sum1s$max), "\n",
  sprintf("Within 1 SD: %d (%.2f%%)", sum1s$n1, 100*sum1s$pct1), "\n",
  sprintf("Within 2 SD: %d (%.2f%%)", sum1s$n2, 100*sum1s$pct2), "\n",
  sprintf("Within 3 SD: %d (%.2f%%)", sum1s$n3, 100*sum1s$pct3), "\n",
  sprintf("PRS >= %.2f (std): %d (%.2f%%)", xcut_std, sum1s$n_cut, 100*sum1s$pct_cut), "\n\n",
  "=== ALL (Standardized to All SNPs) ===\n",
  sprintf("Mean (standardized): %.4f | SD (standardized): %.4f | Min: %.4f | Max: %.4f",
          sum2s$mean, sum2s$sd, sum2s$min, sum2s$max), "\n",
  sprintf("Within 1 SD: %d (%.2f%%)", sum2s$n1, 100*sum2s$pct1), "\n",
  sprintf("Within 2 SD: %d (%.2f%%)", sum2s$n2, 100*sum2s$pct2), "\n",
  sprintf("Within 3 SD: %d (%.2f%%)", sum2s$n3, 100*sum2s$pct3), "\n",
  sprintf("PRS >= %.2f (std): %d (%.2f%%)", xcut_std, sum2s$n_cut, 100*sum2s$pct_cut)
)
writeLines(summary_text_std, "Standardized_PRS_stats.txt")
cat(summary_text_std)
cat("\nStandardized PRS stats written to: Standardized_PRS_stats.txt\n")
