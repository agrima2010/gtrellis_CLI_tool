#For automation - all BED files in one folder

library(gtrellis)
library(readr)
library(grid)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)

# Path to your BED files
bed_dir <- "path/to/file"

# Find all normal BED files
normal_files <- list.files(bed_dir, pattern = "_normal\\.regions\\.bed$", full.names = TRUE)

# Extract sample IDs
sample_ids <- gsub("^.*abh_subject_|_normal\\.regions\\.bed$", "", normal_files)

# Function to process one sample
process_sample <- function(sample_id) {
  normal_path <- file.path(bed_dir, paste0("abh_subject_", sample_id, "_normal.regions.bed"))
  tumor_path  <- file.path(bed_dir, paste0("abh_subject_", sample_id, "_tumor.regions.bed"))
  
  # Read files
  normal_bed <- read_tsv(normal_path, col_names = FALSE)
  tumor_bed  <- read_tsv(tumor_path,  col_names = FALSE)
  
  colnames(normal_bed) <- c("chr", "start", "end", "score")
  colnames(tumor_bed)  <- c("chr", "start", "end", "score")
  
  # Ensure proper types
  normal_bed$chr   <- as.character(normal_bed$chr)
  tumor_bed$chr    <- as.character(tumor_bed$chr)
  normal_bed$start <- as.numeric(normal_bed$start)
  tumor_bed$start  <- as.numeric(tumor_bed$start)
  normal_bed$end   <- as.numeric(normal_bed$end)
  tumor_bed$end    <- as.numeric(tumor_bed$end)
  normal_bed$score <- as.numeric(normal_bed$score)
  tumor_bed$score  <- as.numeric(tumor_bed$score)
  
  # Compute log10 scores
  normal_log <- log10(normal_bed$score + 1)
  tumor_log  <- log10(tumor_bed$score + 1)
  
  ylim_log <- range(tumor_log, normal_log, na.rm = TRUE)
  track_ylim_mat <- matrix(rep(ylim_log, 2), nrow = 2, byrow = TRUE)
  
  # Layout and plotting
  gtrellis_layout(
    n_track = 2,
    nrow = 2,
    compact = TRUE,
    track_ylim = track_ylim_mat,
    track_ylab = c(paste0("log10(score) - Tumor (", sample_id, ")"),
                   paste0("log10(score) - Normal (", sample_id, ")")),
    add_name_track = TRUE,
    add_ideogram_track = TRUE
  )
  
  # Tumor
  tumor_finite <- is.finite(tumor_log)
  add_points_track(
    tumor_bed[tumor_finite, ], tumor_log[tumor_finite],
    pch = 16,
    size = grid::unit(2, "bigpts"),
    gp = grid::gpar(col = "#FF000050")
  )
  
  # Normal
  normal_finite <- is.finite(normal_log)
  add_points_track(
    normal_bed[normal_finite, ], normal_log[normal_finite],
    pch = 16,
    size = grid::unit(2, "bigpts"),
    gp = grid::gpar(col = "#0000FF50")
  )
  
  message("Added plot for sample: ", sample_id)
}

# One PDF for all samples
pdf(file.path(bed_dir, "all_samples_plots.pdf"), width = 8, height = 6)

for (id in sample_ids) {
  process_sample(id)
}

dev.off()
message("✅ All plots saved in all_samples_plots.pdf")
