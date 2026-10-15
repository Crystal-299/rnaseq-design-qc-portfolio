# RNA-seq design and count QC

A small, dependency-free R demonstration of the checks that should precede differential expression: count integrity, sample matching, library size, detected genes, and replicate availability by group and stage.

## Provenance

This repository is a **new portfolio implementation**, informed by an earlier, actually run 30-sample RNA-seq analysis. The source project contains R scripts, dated pipeline logs with completed steps, QC tables, and a scientific summary. Its real counts, sample identifiers, plots, and unpublished results are intentionally absent here. This toy does not reproduce that study or perform DESeq2 inference.

## Requirements

- R 4.2 or newer; base R only.

## Run the synthetic example

From the repository root:

```sh
Rscript scripts/audit.R examples/counts.tsv examples/samples.tsv demo-output
```

Expected files: `demo-output/sample_qc.tsv`, `demo-output/design.tsv`, and `demo-output/summary.txt`. The summary should report `8 samples`, `5 genes`, and both stages as eligible for a replicated two-group comparison. Remove `demo-output` before another run if you want a clean example; the script creates or updates only that output directory.

## Input contract

`counts.tsv` is tab-separated with `gene_id` followed by sample IDs; counts must be finite, nonnegative integers and each library must contain at least one count. `samples.tsv` has exactly one row per count column and columns `sample_id`, `group`, `stage`, `replicate`. IDs must be unique; groups/stages/replicates must be nonempty. Replicate eligibility means **at least two samples per group per stage**; it says nothing about dispersion, batch effects, power, or whether a contrast is biologically valid.

The script rejects missing/extra samples, duplicate gene IDs, bad counts, empty libraries, and missing metadata. Inspect the printed error and correct the input; there is no automatic sample exclusion.

## Interpretation

The example demonstrates an audit gate, not a complete RNA-seq pipeline. In the original work, sample exclusions, formal interaction models, and biological interpretation depended on project-specific evidence. See [LIMITATIONS.md](LIMITATIONS.md).
