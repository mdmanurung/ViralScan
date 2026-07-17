#!/usr/bin/env Rscript
# emptydrops.R — call cells from a raw kb-count matrix using DropletUtils::emptyDrops.
#
# ViralScan (Python) shells out to this script so viral rates can be reported over
# *called cells* (real, non-empty droplets), not over all barcodes — the latter
# produces meaningless "% infected" denominators (see .living finding F-005; the
# HSV-1 P22.5 denominator artifact is the same trap).
#
# Input is the kb-count `counts_unfiltered` directory, whose MatrixMarket matrix is
# barcodes (rows) x genes (cols). emptyDrops expects genes x cells, so we transpose.
#
# Usage:
#   Rscript emptydrops.R <counts_unfiltered_dir> <out_tsv> [fdr] [lower] [niters] [seed]
#
# Output TSV columns: barcode  total  FDR  is_cell  knee  inflection
#   is_cell = (FDR <= fdr) & !is.na(FDR), OR total >= knee (always-cell above the knee).

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
}))

if (length(args) < 2) {
  stop("Usage: emptydrops.R <counts_unfiltered_dir> <out_tsv> [fdr] [lower] [niters] [seed]")
}

mtx_dir   <- args[[1]]
out_tsv   <- args[[2]]
fdr_thr   <- if (length(args) >= 3) as.numeric(args[[3]]) else 0.01
lower     <- if (length(args) >= 4) as.numeric(args[[4]]) else 100
niters    <- if (length(args) >= 5) as.integer(args[[5]]) else 10000L
seed      <- if (length(args) >= 6) as.integer(args[[6]]) else 100L

# Optional custom R library path (e.g. where DropletUtils is installed). Set
# VIRALSCAN_R_LIBS to prepend a location; R's standard R_LIBS_USER is honoured otherwise.
user_lib <- Sys.getenv("VIRALSCAN_R_LIBS", unset = "")
if (nzchar(user_lib) && dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressWarnings(suppressMessages({
  library(Matrix)
  library(DropletUtils)
}))

# --- locate the matrix triplet (kb-python naming) --------------------------------
find_one <- function(dir, patterns) {
  for (p in patterns) {
    hit <- list.files(dir, pattern = p, full.names = TRUE)
    if (length(hit) >= 1) return(hit[[1]])
  }
  stop(sprintf("No file matching %s in %s", paste(patterns, collapse="/"), dir))
}
mtx_f  <- find_one(mtx_dir, c("\\.mtx$", "\\.mtx\\.gz$"))
bc_f   <- find_one(mtx_dir, c("\\.barcodes\\.txt$", "barcodes\\.txt$", "barcodes\\.tsv"))

message(sprintf("[emptydrops] matrix=%s", mtx_f))
m <- readMM(mtx_f)                    # barcodes x genes (kb cells_x_genes)
barcodes <- readLines(bc_f)
if (nrow(m) != length(barcodes)) {
  stop(sprintf("barcode count (%d) != matrix rows (%d)", length(barcodes), nrow(m)))
}
gxc <- t(m)                          # genes x cells for emptyDrops
gxc <- as(gxc, "CsparseMatrix")
message(sprintf("[emptydrops] %d barcodes x %d genes; lower=%g fdr=%g niters=%d",
                ncol(gxc), nrow(gxc), lower, fdr_thr, niters))

set.seed(seed)
ed <- emptyDrops(gxc, lower = lower, niters = niters)
br <- barcodeRanks(gxc, lower = lower)
knee <- metadata(br)$knee
infl <- metadata(br)$inflection

total <- Matrix::colSums(gxc)
FDR   <- ed$FDR
is_cell <- (!is.na(FDR) & FDR <= fdr_thr) | (total >= knee)

out <- data.frame(
  barcode = barcodes,
  total = as.numeric(total),
  FDR = as.numeric(FDR),
  is_cell = is_cell,
  knee = knee,
  inflection = infl,
  stringsAsFactors = FALSE
)
write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[emptydrops] called cells: %d / %d  (knee=%.0f, inflection=%.0f)",
                sum(is_cell, na.rm = TRUE), nrow(out), knee, infl))
cat("EMPTYDROPS_DONE:", sum(is_cell, na.rm = TRUE), "cells\n")
