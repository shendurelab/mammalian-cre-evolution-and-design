#!/usr/bin/env Rscript
# Preprocess large external chromBPNet HDF5 files into lightweight TSVs.
#
# The prediction H5 files are ~1.5 GB each (4 cells x 2 peak sets), far too
# large to include in the repo. This script reads only the predictions/logcounts
# dataset (small) and merges it with the *_logcts_observed.bed columns,
# producing per-cell TSVs that the figure script can load directly.
#
# Run once from the repo root:
#   Rscript scripts/utils/preprocess_ext_data_fig1_h5.R
#
# Outputs (gzipped; ~a few MB each):
#   data/ext_data_fig1/consensus_peak_preds/{cell}_cons_peaks_obs_pred.txt.gz
#   data/ext_data_fig1/da_peak_preds/{cell}_da_peaks_obs_pred.txt.gz
#   data/ext_data_fig1/gkm_svm/{endo_pos,endo_neg}_chrombpnet_logcounts.txt.gz

suppressPackageStartupMessages({
  library(rhdf5)
})

# ---- Source paths (external; required only during preprocessing) ----
CONS_DIR <- "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet/chrombpnet_overlap_peaks/peaks_prediction"
DA_DIR   <- "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/EB_chrombpnet/DA_peaks_prediction"
GKM_DIR  <- "/net/shendure/vol10/projects/tli/nobackup/EB_dev_project/run_gkmSVM/chrombpnet_matched_gkmSVM_predictions"

OUT_BASE <- "data/ext_data_fig1"
dir.create(file.path(OUT_BASE, "consensus_peak_preds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_BASE, "da_peak_preds"),        recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_BASE, "gkm_svm"),              recursive = TRUE, showWarnings = FALSE)

CELLS   <- c("endo", "pluri", "meso", "ecto")
BED_COLS <- c("chr","start","end","logcts","pred_logcts","V6","V7","V8","qval","summit")

extract_peak_preds <- function(h5_file, bed_file, out_file) {
  cat("  ", basename(h5_file), "\n")
  logcounts <- h5read(h5_file, "predictions/logcounts")
  obs <- read.delim(bed_file, sep = "\t", header = FALSE)
  obs$V5 <- logcounts
  colnames(obs) <- BED_COLS
  write.table(obs, gzfile(out_file), sep = "\t", row.names = FALSE, quote = FALSE)
}

# --- Consensus peaks ---
cat("Extracting consensus peak predictions...\n")
for (cell in CELLS) {
  extract_peak_preds(
    h5_file  = file.path(CONS_DIR, paste0(cell, "_chrombpnet_nobias_predictions.h5")),
    bed_file = file.path(CONS_DIR, paste0(cell, "_logcts_observed.bed")),
    out_file = file.path(OUT_BASE, "consensus_peak_preds",
                         paste0(cell, "_cons_peaks_obs_pred.txt.gz"))
  )
}

# --- DA peaks ---
cat("Extracting differentially accessible (DA) peak predictions...\n")
for (cell in CELLS) {
  extract_peak_preds(
    h5_file  = file.path(DA_DIR, paste0(cell, "_chrombpnet_nobias_predictions.h5")),
    bed_file = file.path(DA_DIR, paste0(cell, "_logcts_observed.bed")),
    out_file = file.path(OUT_BASE, "da_peak_preds",
                         paste0(cell, "_da_peaks_obs_pred.txt.gz"))
  )
}

# --- gkmSVM-matched chrombpnet predictions (endo only; for ext_data_fig3_partIII) ---
cat("Extracting gkmSVM-matched chromBPNet logcounts...\n")
for (pn in c("pos", "neg")) {
  h5_file  <- file.path(GKM_DIR, paste0("endo_", pn, "_chrombpnet_nobias_predictions.h5"))
  out_file <- file.path(OUT_BASE, "gkm_svm",
                        paste0("endo_", pn, "_chrombpnet_logcounts.txt.gz"))
  cat("  ", basename(h5_file), "\n")
  logcounts <- h5read(h5_file, "predictions/logcounts")
  write.table(data.frame(pred = logcounts), gzfile(out_file),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("\nDone.\n")
