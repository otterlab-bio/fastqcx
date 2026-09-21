# fastqcx

[![CI](https://github.com/otterlab-bio/fastqcx/actions/workflows/rust.yml/badge.svg)](https://github.com/otterlab-bio/fastqcx/actions/workflows/rust.yml)

**A Rust FASTQ quality-control reporter with self-contained HTML and FastQC-compatible summary output.**

`fastqcx` reads `.fastq`, `.fastq.gz`, `.fq`, and `.fq.gz` files, reports read length and sequence-quality statistics, and can emit `fastqc_data.txt` for MultiQC and downstream QC aggregation.

## Proof

One file, two records, written through the `--summary` contract that MultiQC and `qctb` consume:

```text
$ fastqcx --fastq mini.fastq --summary qc/mini --no-html
$ head -12 qc/mini/fastqc_data.txt
##FastQC	0.3.4
>>Basic Statistics	pass
#Measure	Value	
Filename	mini.fastq
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	2
Filtered Sequences	0	
Sequence length	32
%GC	50
>>END_MODULE
>>Per base sequence quality	warn
```

`%GC` is 50 for a 50% GC file, the encoding is recognised rather than assumed, and the per-base
module carries its own pass/warn verdict. The module labels are FastQC compatibility surface, not a
claim to reproduce every FastQC implementation detail.

## Included metrics

- Read length distribution
- Per-base sequence quality
- Per-base sequence content
- Sequence duplication levels using FastQC's bounded 50-base identity contract
- FastQC-style overrepresented sequence identities under the same bounded identity contract
- K-mer content
- GC content
- Sequence quality score summaries

## Limits and input validation

Whole-file totals (reads, bases, GC, quality summaries) are always complete.
Per-position charts are bounded by a tracking limit (default: the first
10,000 positions per read; `AnalysisLimits` in the library API): positions
beyond the bound are excluded from per-base charts only, and the HTML report
flags when the bound was hit. The bound is a library-level knob and is not
currently exposed as a CLI flag.

Inputs are validated fail-closed:

- Quality bytes must be within the supported Phred+33 range (`!` through `~`). The parser cannot infer an encoding when byte ranges overlap.
- Empty FASTQ files and empty sequences are rejected.
- Reads longer than 10,000,000 bases or more than 100,000 distinct read lengths are rejected.
- `-k/--kmer` accepts 1–7.

## Install

Build from source. `fastqcx` is not published to crates.io, so `cargo install fastqcx` fails with
`crate fastqcx does not exist`:

```bash
git clone https://github.com/otterlab-bio/fastqcx.git
cd fastqcx
cargo build --release
```

## Quick start

Write an HTML report to standard output:

```bash
fastqcx -q sample.fastq.gz > sample.fastqcx.html
```

Write the FastQC-compatible summary directory as well:

```bash
fastqcx \
  --fastq sample.fastq.gz \
  --summary qc/sample \
  > qc/sample.html
```

Generate compatibility output only:

```bash
fastqcx --fastq sample.fastq.gz --summary qc/sample --no-html
```

For multiple files, keep each summary directory separate:

```bash
for fastq_file in *.fastq.gz; do
  sample_name="${fastq_file%.fastq.gz}"
  fastqcx --fastq "$fastq_file" --summary "qc/$sample_name" > "reports/$sample_name.html"
done
```

## CLI contract

| Option | Meaning |
| --- | --- |
| `-q, --fastq` | Required FASTQ input path. Compressed input is supported. |
| `-k, --kmer` | K-mer length from 1 to 7; default is 5. |
| `-s, --summary` | Directory for FastQC-compatible `fastqc_data.txt`. |
| `--no-html` | Suppress HTML output on standard output. |

## Output contract

- HTML is written to standard output unless `--no-html` is supplied.
- `--summary` creates a FastQC-compatible summary directory for MultiQC and `qctb` consumers.
- Compatibility labels such as FastQC and SeqKit remain external protocol names; `fastqcx` does not claim to reproduce every FastQC implementation detail.

## Compatibility evidence

The GitHub Actions compatibility job runs `fastqcx`, FastQC 0.12.1, and SeqKit 2.13.0 on the
same pinned RNA-PDX `SRR30880970` R1 fixture. A deterministic duplicated-read positive control
additionally requires a non-empty overrepresented-sequence result from both implementations. The
retained evidence compares:

- 11 SeqKit integer fields exactly and 5 decimal fields within declared field-specific tolerances;
- every grouped FastQC per-base quality value (mean, median, quartiles, and 10th/90th percentiles)
  after applying FastQC's position bins;
- all 16 FastQC sequence-duplication bins and the total deduplicated percentage; and
- overrepresented sequence identities and counts exactly, with percentages within 0.01 percentage
  points, including the non-empty positive control.

Module verdicts must also agree. The resulting `otter.fastqcx-parity/v2` JSON report and all three
raw outputs are uploaded as the `fastqcx-rna-pdx-compatibility-evidence` workflow artifact. This
contract is deliberately narrower than complete FastQC equivalence: adapter content, per-tile
quality, and other unlisted FastQC modules remain outside the validated surface.

## Development

```bash
cargo fmt --all -- --check
cargo test
cargo build --release
```

## License and repository

MIT · [otterlab-bio/fastqcx](https://github.com/otterlab-bio/fastqcx)
