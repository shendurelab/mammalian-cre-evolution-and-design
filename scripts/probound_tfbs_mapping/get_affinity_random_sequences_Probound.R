### Regenerate the null-distribution affinity quantiles (q99, q99p9, q99p99) in
### data/ParEndo_markerTF_cuttoffs.txt.
###
### For each ParEndo marker TF, scores a large set of random DNA sequences (at
### the TF's ProBound model size) with ProBound, then takes the 99th / 99.9th /
### 99.99th percentile of the affinity distribution as that TF's background
### cutoffs. TF_name/TFID/model_size/max_aff are carried over unchanged from
### the committed table: TFID/model_size come from the ProBound Central motif
### models, and max_aff is an empirical calibration constant (max observed
### affinity across curated maximal-tile DMS CRE sequences) that this script
### does not re-derive.
###
### Usage: Rscript scripts/probound_tfbs_mapping/get_affinity_random_sequences_Probound.R
### Requires: Java 8+ and a built ProBound jar (see scripts/probound_tfbs_mapping/README.md),
### exposed via $PROBOUND_JAR.

suppressMessages(library(dplyr))

RNG_SEED <- 20240414
N_RANDOM <- 1e6 # sequences per model size; override for a quick smoke test
BASES <- c("A", "C", "G", "T")

### resolve repo-root-relative paths regardless of the working directory the
### script is launched from (mirrors run_probound.sh's REPO_ROOT handling)
.script_dir <- (function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
    if (length(file_arg) == 0) {
        return(getwd()) # interactive session fallback
    }
    dirname(normalizePath(file_arg))
})()
REPO_ROOT <- normalizePath(file.path(.script_dir, "..", ".."))

TF_MODEL_TABLE <- file.path(REPO_ROOT, "data", "ParEndo_markerTF_cuttoffs.txt")
OUT_TABLE <- TF_MODEL_TABLE

probound_jar <- Sys.getenv("PROBOUND_JAR")
if (!nzchar(probound_jar)) {
    stop("Set $PROBOUND_JAR to the ProBound jar path (see scripts/probound_tfbs_mapping/README.md)")
}
if (!file.exists(probound_jar)) {
    stop("PROBOUND_JAR does not point to an existing file: ", probound_jar)
}
if (nzchar(Sys.which("java")) == FALSE) {
    stop("java not found on PATH")
}

df_TF_models <- read.table(TF_MODEL_TABLE, header = TRUE, sep = "\t")
rownames(df_TF_models) <- df_TF_models$TF_name

### score N_RANDOM random sequences of a given length with ProBound's consensus
### binding-mode score, returning a numeric vector of affinities
score_random_sequences <- function(model_size, tfid, n = N_RANDOM) {
    set.seed(RNG_SEED + model_size) # reproducible per model size, independent of TF order
    random_seqs <- replicate(n, paste0(sample(BASES, model_size, replace = TRUE), collapse = ""))

    seq_file <- tempfile()
    writeLines(random_seqs, seq_file)
    on.exit(unlink(seq_file), add = TRUE)

    cmd <- sprintf(
        "java -cp %s proBoundTools/App -c 'loadMotifCentralModel(%d).buildConsensusModel().addNScoring().inputTXT(%s).bindingModeScores(/dev/stdout)'",
        shQuote(probound_jar), tfid, shQuote(seq_file)
    )
    affinity_out <- system(cmd, intern = TRUE)
    df_affinity <- read.table(text = affinity_out, col.names = c("substring", "affinity"))
    df_affinity$affinity
}

message("Scoring ", N_RANDOM, " random sequences per TF model against ProBound...")
df_quantiles <- lapply(df_TF_models$TF_name, function(TF_oi) {
    message("  ", TF_oi)
    affinity <- score_random_sequences(df_TF_models[TF_oi, "model_size"], df_TF_models[TF_oi, "TFID"])
    data.frame(
        TF_name = TF_oi,
        q99p99 = quantile(affinity, 0.9999, names = FALSE),
        q99p9 = quantile(affinity, 0.999, names = FALSE),
        q99 = quantile(affinity, 0.99, names = FALSE)
    )
}) %>% bind_rows()

df_out <- df_TF_models %>%
    dplyr::select(TF_name, TFID, model_size, max_aff) %>%
    dplyr::left_join(df_quantiles, by = "TF_name") %>%
    dplyr::select(TF_name, TFID, model_size, q99p99, q99p9, q99, max_aff)

write.table(df_out, OUT_TABLE, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote ", nrow(df_out), " TF cutoffs to ", OUT_TABLE)
