#!/usr/bin/env Rscript


# ---- Load libraries ----
library(gtrellis)
library(readr)
library(grid)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)


# ---- Parse command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_bed_tracks.R <normal_bed_file> <tumor_bed_file>")
}
normal_file <- args[1]
tumor_file  <- args[2]

# ---- Read BED files ----
normal_bed <- read_tsv(normal_file, col_names = FALSE)
tumor_bed  <- read_tsv(tumor_file, col_names = FALSE)

# ---- Assign column names ----
colnames(normal_bed) <- c("chr", "start", "end", "score")
colnames(tumor_bed)  <- c("chr", "start", "end", "score")

# ---- Ensure proper types ----
normal_bed <- as.data.frame(normal_bed)
tumor_bed  <- as.data.frame(tumor_bed)

normal_bed$chr   <- as.character(normal_bed$chr)
tumor_bed$chr    <- as.character(tumor_bed$chr)
normal_bed$start <- as.numeric(normal_bed$start)
tumor_bed$start  <- as.numeric(tumor_bed$start)
normal_bed$end   <- as.numeric(normal_bed$end)
tumor_bed$end    <- as.numeric(tumor_bed$end)
normal_bed$score <- as.numeric(normal_bed$score)
tumor_bed$score  <- as.numeric(tumor_bed$score)

# ---- Compute log10-transformed scores ----
normal_bed_log_score <- log10(normal_bed$score + 1)
normal_log_score_finite <- normal_bed_log_score[is.finite(normal_bed_log_score)]
tumor_bed_log_score  <- log10(tumor_bed$score + 1)
tumor_log_score_finite <- tumor_bed_log_score[is.finite(tumor_bed_log_score)]

# ---- Get y-limits from both datasets ----
ylim_log <- range(tumor_bed_log_score, normal_bed_log_score, na.rm = TRUE)

# Prepare track_ylim with one row per track
track_ylim_mat <- matrix(rep(ylim_log, 2), nrow = 2, byrow = TRUE)

# ---- Open PDF for the result ----
output_pdf <- "bed_tracks_plot.pdf"
pdf(output_pdf, width = 10, height = 6)

# ---- Create gtrellis layout ----
gtrellis_layout(
  n_track = 2,
  nrow = 2,
  compact = TRUE,
  track_ylim = track_ylim_mat,
  track_ylab = c("log10(score) - Tumor", "log10(score) - Normal"),
  add_name_track = TRUE,
  add_ideogram_track = TRUE
)

# ---- Add tumor data (first track) ----
tumor_finite <- is.finite(tumor_bed_log_score)
add_points_track(
  tumor_bed[tumor_finite, ], tumor_bed_log_score[tumor_finite],
  pch = 16,
  size = grid::unit(2, "bigpts"),
  gp = grid::gpar(col = "#FF000050")
)

# ---- Add normal data (second track) ----
normal_finite <- is.finite(normal_bed_log_score)
add_points_track(
  normal_bed[normal_finite, ], normal_bed_log_score[normal_finite],
  pch = 16,
  size = grid::unit(2, "bigpts"),
  gp = grid::gpar(col = "#0000FF50")
)


# ---- Close the PDF ----
dev.off()
cat("Plot saved to", output_pdf, "\n")
