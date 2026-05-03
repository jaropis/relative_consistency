# Master analysis script for SampEn curve intersection study
#
# This script orchestrates the full analysis pipeline:
#   1. Load SampEn data for all RR-interval recordings.
#   2. Detect pairwise intersections of the SampEn curves.
#   3. Export summary figures and tabular results.

source("samp_ens.R")

# --- Data loading and intersection detection ----------------------------------

master_analysis <- get_data_paths("data/RR-SampEn-data") %>%
  get_figures_and_results(., sep = "\t")

# --- Figure: all SampEn curves ------------------------------------------------

png("all_RR_samp.png", width = 2200, height = 1200, res = 300)
print(master_analysis$plot)
dev.off()

# --- Crossing-point summary table ---------------------------------------------

tabelka_crosses <- summarize_crosses(master_analysis$num_results)

# --- Figure: histogram of intersection counts per recording pair --------------

png("hist_RR_crosses_no.png", width = 1400, height = 1200, res = 300)
cross_hist <- sapply(master_analysis$num_results, function(elem) {
  length(elem[elem != "pusta"])
}) %>%
  hist(
    xlab = "number of intersections between two curves",
    ylab = "number of events",
    main = "",
    xlim = c(1, 85),
    breaks = seq(0, 85)
  )
dev.off()
print(cross_hist$counts)

# --- Figure: histogram of intersection tolerance values -----------------------

png("hist_RR_crosses_r.png", width = 1400, height = 1200, res = 300)
r_values <- hist(
  c(tabelka_crosses[["r1"]], tabelka_crosses[["r2"]]),
  xlab = expression(r),
  ylab = "number of intersections",
  main = ""
)
dev.off()

# --- Figure: combined panel (tolerance histogram + intersection-count histogram)

png("hist_rr_crosses_all.png", width = 2800, height = 1200, res = 300)
par(mfrow = c(1, 2))

r_breaks <- seq(0, ceiling(max(c(tabelka_crosses[["r1"]], tabelka_crosses[["r2"]]), na.rm = TRUE) / 0.25) * 0.25, by = 0.25)
r_values <- hist(
  c(tabelka_crosses[["r1"]], tabelka_crosses[["r2"]]),
  breaks = r_breaks,
  xlab = expression(rho),
  ylab = "number of intersections",
  main = "",
  col = c("grey80", rep("white", length(r_breaks)))
)

cross_hist <- sapply(master_analysis$num_results, function(elem) {
  length(elem[elem != "pusta"])
}) %>%
  hist(
    xlab = "number of intersections between two curves",
    ylab = "number of events",
    main = "",
    xlim = c(1, 85),
    breaks = seq(0, 85)
  )

dev.off()
print(cross_hist$counts)
