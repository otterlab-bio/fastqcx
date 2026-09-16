# CLAUDE.md

This file provides guidance for working with fastqcx.

## Project Overview

fastqcx is a Rust FASTQ quality-control library and CLI. The `fastqcx` binary emits an HTML report to stdout unless `--no-html` is supplied and optionally writes `<summary-directory>/fastqc_data.txt` for MultiQC and qctb consumers.

## Build and Quality Gates

```bash
cargo fmt --all -- --check
cargo check --all-targets --all-features --locked
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --all-targets --all-features --locked
cargo build --all-targets --all-features --locked
```

## Running the Tool

```bash
# HTML report to stdout
fastqcx --fastq input.fastq > report.html

# Workflow/qctb output only
fastqcx --fastq input.fastq --summary summary_directory --no-html

# Supported k-mer lengths are 1 through 7
fastqcx --fastq input.fastq --kmer 7 --summary summary_directory --no-html
```

Arguments:

- `-q/--fastq`: required FASTQ input; FASTA, empty reads, malformed records, truncated streams, and invalid Phred+33 qualities fail closed.
- `-k/--kmer`: k-mer length in the inclusive range 1–7; default 5.
- `-s/--summary`: directory receiving `fastqc_data.txt`.
- `--no-html`: skip HTML output.

## Typed Analysis Contract

`src/types.rs` defines:

- `KmerLength`: validated 1–7 k-mer bound.
- `AnalysisLimits`: resource bounds for tracked positions, maximum read length, and distinct read lengths.
- `AnalysisOptions`: validated options passed to `analyze_with_options`.
- `FastqAnalysisError`: fail-closed input and resource errors.
- `FastqAnalysisResult`: whole-file totals plus bounded per-position histograms.

Default limits:

- At most 10,000 positions are retained for per-position base and quality charts.
- Reads longer than 10,000,000 bases are rejected.
- More than 100,000 distinct read lengths are rejected.
- Whole-file base, quality, Q20/Q30, GC, read-length, and SeqKit totals remain complete when per-position charts are truncated.

## Processing Architecture

- `analyze()` validates the k-mer length and uses default limits.
- `analyze_with_options()` is the single FASTQ parsing and accumulation pipeline.
- `process()` renders HTML and summary output from `FastqAnalysisResult`; it does not parse the FASTQ independently.
- Needletail parse errors are returned immediately. Partial reports are never published after malformed or truncated input.

## Report Compatibility

- HTML generation performs no network requests. The generated document references pinned CDN URLs, so report generation is deterministic and offline-safe; browser-side charts require those external assets when the report is opened.
- `fastqc_data.txt` starts with standard `##FastQC` and Basic Statistics fields understood by MultiQC.
- The custom `>>Seqkit Statistics` module is consumed by qctb and follows SeqKit 2.13 semantics:
  - Q20/Q30 are base-level percentages.
  - AvgQual is derived from mean error probability, not arithmetic mean Phred.
  - SeqKit GC uses total sequence length as denominator; FastQC Basic Statistics `%GC` excludes ambiguous bases.
  - `N50_num`, length quartiles, `sum_gap`, and `sum_n` match SeqKit output behavior.
- Per-base sequence content, N content, GC, and quality positions are 1-based and percentage fields use the 0–100 scale expected by FastQC/MultiQC.

## Tests

- Unit tests in `src/process.rs` cover plain and gzip truncation, malformed FASTQ, FASTA rejection, CRLF, empty reads, IUPAC bases, low quality, variable lengths, k-mer bounds, resource limits, and SeqKit statistics.
- `tests/lib.rs` runs the real Cargo-built `fastqcx` binary and verifies CLI exit behavior, offline HTML generation, and summary contracts.
- `.github/workflows/rust.yml` pins and runs real Babraham FastQC 0.12.1, SeqKit 2.13.0, and MultiQC 1.35 compatibility smoke tests.

## Dependencies

Runtime dependencies are deliberately local-only:

- `needletail`: FASTA/Q parsing and gzip support.
- `clap`: CLI parsing.
- `tera`: report templates.
- `serde_json`: Vega-Lite data.
- `rustc-hash`: bounded analysis maps.
- `itertools`: deterministic sorting helpers.

`flate2` is a dev dependency used to generate truncated gzip fixtures in tests.
