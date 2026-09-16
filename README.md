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
- K-mer content
- GC content
- Sequence quality score summaries

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

## Development

```bash
cargo fmt --all -- --check
cargo test
cargo build --release
```

## License and repository

MIT · [otterlab-bio/fastqcx](https://github.com/otterlab-bio/fastqcx)
