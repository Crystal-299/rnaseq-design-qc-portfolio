args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript scripts/audit.R COUNTS.tsv SAMPLES.tsv OUTPUT_DIR", call. = FALSE)
counts_path <- args[[1]]
samples_path <- args[[2]]
out_dir <- args[[3]]
if (!file.exists(counts_path) || !file.exists(samples_path)) stop("Input file is missing", call. = FALSE)
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
counts <- read_tsv(counts_path)
samples <- read_tsv(samples_path)
required <- c("sample_id", "group", "stage", "replicate")
if (ncol(counts) < 2L || names(counts)[1L] != "gene_id") stop("Counts need gene_id and sample columns", call. = FALSE)
if (!all(required %in% names(samples))) stop("Metadata is missing required columns", call. = FALSE)
if (nrow(counts) == 0L || nrow(samples) == 0L) stop("Input tables must be nonempty", call. = FALSE)
if (anyNA(counts$gene_id) || anyDuplicated(counts$gene_id)) stop("Gene IDs must be unique and nonempty", call. = FALSE)
for (field in required) if (anyNA(samples[[field]]) || any(!nzchar(trimws(as.character(samples[[field]]))))) stop(paste("Missing", field), call. = FALSE)
if (anyDuplicated(samples$sample_id) || anyDuplicated(names(counts))) stop("Duplicate sample ID or count column", call. = FALSE)
sample_ids <- names(counts)[-1L]
if (!setequal(sample_ids, samples$sample_id) || length(sample_ids) != nrow(samples)) stop("Count and metadata sample IDs do not match exactly", call. = FALSE)
samples <- samples[match(sample_ids, samples$sample_id), required, drop = FALSE]
replicate_key <- paste(samples$group, samples$stage, samples$replicate, sep = "\r")
if (anyDuplicated(replicate_key)) stop("Duplicate replicate within a group and stage", call. = FALSE)
raw <- counts[, -1L, drop = FALSE]
values <- suppressWarnings(matrix(as.numeric(as.matrix(raw)), nrow = nrow(raw), dimnames = list(NULL, sample_ids)))
if (anyNA(values) || any(!is.finite(values)) || any(values < 0) || any(values != round(values))) stop("Counts must be finite, nonnegative integers", call. = FALSE)
lib <- colSums(values)
if (any(lib == 0)) stop("At least one sample has an empty library", call. = FALSE)
qc <- data.frame(sample_id = sample_ids, group = samples$group, stage = samples$stage,
                 replicate = samples$replicate, library_size = lib,
                 detected_genes = colSums(values > 0), zero_fraction = colMeans(values == 0),
                 check.names = FALSE)
design <- aggregate(rep(1L, nrow(samples)), by = list(stage = samples$stage, group = samples$group), FUN = sum)
names(design)[3L] <- "n_samples"
design <- design[order(design$stage, design$group), , drop = FALSE]
stages <- unique(as.character(samples$stage))
eligible <- stages[vapply(stages, function(s) {
  rows <- design[design$stage == s, , drop = FALSE]
  nrow(rows) >= 2L && all(rows$n_samples >= 2L)
}, logical(1))]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(qc, file.path(out_dir, "sample_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(design, file.path(out_dir, "design.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
lines <- c(sprintf("%d samples", length(sample_ids)), sprintf("%d genes", nrow(counts)),
           sprintf("eligible replicated stages: %s", if (length(eligible)) paste(eligible, collapse = ", ") else "none"),
           "Eligibility is structural only; no statistical model was fitted.")
writeLines(lines, file.path(out_dir, "summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
