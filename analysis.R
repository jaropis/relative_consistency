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

# --- Figure: enlarged view of the recommended tolerance range ----------------

zoom_limits <- c(0.1, 0.25)
zoom_data <- master_analysis$plot$data
zoom_y_range <- range(
  zoom_data$SampEn[zoom_data$r >= zoom_limits[1] &
                     zoom_data$r <= zoom_limits[2]],
  finite = TRUE
)

# Keep the complete profiles so lines are clipped, rather than truncated,
# at the boundaries of the enlarged view.
rr_zoom <- ggplot(zoom_data, aes(r, SampEn, group = group)) +
  geom_line(linewidth = 0.25) +
  scale_x_continuous(
    name = expression(rho),
    breaks = seq(0.1, 0.25, by = 0.025),
    labels = scales::label_number(accuracy = 0.001),
    expand = expansion(mult = 0)
  ) +
  coord_cartesian(xlim = zoom_limits, ylim = zoom_y_range) +
  master_analysis$plot$theme +
  theme(plot.margin = margin(5.5, 16, 5.5, 5.5))

png("../submission2/all_RR_samp_zoom.png", width = 2200, height = 1200, res = 300)
print(rr_zoom)
dev.off()

# --- Crossing-point summary table ---------------------------------------------

tabelka_crosses <- summarize_crosses(master_analysis$num_results)

# --- Histogram input: all retained crossing events ----------------------------

# The summary table above keeps at most two crossings per pair. Histograms of
# crossing locations must instead use every event returned by the detector.
crossing_rho <- unlist(lapply(master_analysis$num_results, function(pair) {
  vapply(Filter(is.list, pair), function(crossing) crossing$r, numeric(1))
}), use.names = FALSE)

plot_crossing_locations <- function(crossing_rho) {
  r_breaks <- sort(unique(c(
    0.1,
    seq(0, ceiling(max(crossing_rho) / 0.25) * 0.25, by = 0.25)
  )))
  r_values <- hist(crossing_rho, breaks = r_breaks, plot = FALSE)
  bin_left <- head(r_values$breaks, -1)
  bin_right <- tail(r_values$breaks, -1)
  recommended_bin <- bin_left == 0.1 & bin_right == 0.25
  band_count <- sum(crossing_rho >= 0.1 & crossing_rho <= 0.25)

  # Draw counts as bar heights explicitly: unequal bin widths isolate the
  # recommended interval, so the vertical axis is not a probability density.
  plot(
    NA_real_, NA_real_, type = "n",
    xlim = range(r_breaks), ylim = c(0, max(r_values$counts) * 1.12),
    xlab = expression(rho), ylab = "number of crossing events per bin",
    xaxt = "n", yaxs = "i"
  )
  rect(
    bin_left, 0, bin_right, r_values$counts,
    col = ifelse(recommended_bin, "grey80", "white")
  )
  axis(1, at = seq(0, max(r_breaks), by = 0.5), cex.axis = 0.9)
  axis(1, at = c(0.1, 0.25), labels = FALSE)
  # Put the narrow-bin labels on a separate line, aligned away from each other.
  mtext("0.10", side = 1, at = 0.1, line = 1.7, adj = 1, cex = 0.8)
  mtext("0.25", side = 1, at = 0.25, line = 1.7, adj = 0, cex = 0.8)
  legend(
    "topright", fill = "grey80", bty = "n", cex = 0.9,
    legend = paste0("0.1-0.25: ", band_count, " crossing events")
  )
  invisible(r_values)
}

# --- Frequencies of retained crossings per recording pair ---------------------

crossings_per_pair <- vapply(master_analysis$num_results, function(pair) {
  length(Filter(is.list, pair))
}, integer(1))

plot_crossing_counts <- function(crossings_per_pair, max_crossings = 40L) {
  frequency <- tabulate(crossings_per_pair,
                        nbins = max(max_crossings, crossings_per_pair))
  shown <- seq_len(max_crossings)
  nonzero <- frequency[shown] > 0
  y_ticks <- 10^seq.int(0, ceiling(log10(max(frequency[shown]))))
  plot(
    shown[nonzero], frequency[shown][nonzero],
    log = "y", pch = 16, cex = 0.8,
    xlim = c(1, max_crossings), ylim = range(y_ticks),
    xlab = "number of crossings per pair",
    ylab = "number of recording pairs (log scale)",
    xaxt = "n", yaxt = "n"
  )
  axis(1, at = c(1, seq(10, max_crossings, by = 10)))
  axis(2, at = y_ticks, labels = y_ticks)
  invisible(frequency)
}

# --- Figure: crossing-count frequencies on a logarithmic vertical axis --------

png("hist_RR_crosses_no.png", width = 1400, height = 1200, res = 300)
cross_frequency <- plot_crossing_counts(crossings_per_pair)
dev.off()
print(cross_frequency)

# --- Figure: histogram of intersection tolerance values -----------------------

png("hist_RR_crosses_r.png", width = 1400, height = 1200, res = 300)
r_values <- plot_crossing_locations(crossing_rho)
dev.off()

# --- Figure: combined panel (tolerance histogram + crossing-count frequencies)

png("../submission2/hist_RR_crosses_all.png", width = 2800, height = 1200, res = 300)
par(mfrow = c(1, 2))

r_values <- plot_crossing_locations(crossing_rho)

cross_frequency <- plot_crossing_counts(crossings_per_pair)

dev.off()
print(cross_frequency)
