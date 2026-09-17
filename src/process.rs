use crate::types::{
    AnalysisLimits, AnalysisOptions, FastqAnalysisError, FastqAnalysisResult, KmerLength,
    QualityTotals, SequenceLengthSummary,
};
use chrono::{DateTime, Local};
use itertools::Itertools;
use needletail::{parse_fastx_file, Sequence};
use rustc_hash::FxHashMap as HashMap;
use serde_json::json;
use serde_json::Value;
use std::convert::TryFrom;
use std::error::Error;
use std::ffi::OsStr;
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::Path;
use tera::{Context, Tera};

const BASES: [char; 5] = ['A', 'C', 'G', 'T', 'N'];
const A: usize = 0;
const C: usize = 1;
const G: usize = 2;
const T: usize = 3;
const N: usize = 4;

fn base_index(base: u8) -> usize {
    match base {
        b'A' | b'a' => A,
        b'C' | b'c' => C,
        b'G' | b'g' => G,
        b'T' | b't' => T,
        _ => N,
    }
}

fn gc_percentage(sequence: &[u8]) -> usize {
    let gc_bases = sequence
        .iter()
        .filter(|&&base| matches!(base, b'G' | b'g' | b'C' | b'c'))
        .count();
    let called_bases = sequence
        .iter()
        .filter(|&&base| matches!(base, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'))
        .count();

    gc_bases
        .saturating_mul(100)
        .checked_div(called_bases)
        .unwrap_or(0)
}

fn phred33_score(encoded_quality: u8) -> io::Result<usize> {
    let score = encoded_quality.checked_sub(33).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("FASTQ quality byte {encoded_quality} is below Phred+33"),
        )
    })?;
    if score > 93 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("FASTQ quality byte {encoded_quality} exceeds supported Phred+33 range"),
        ));
    }
    Ok(score as usize)
}

fn quality_percentiles(histogram: &[usize]) -> [f32; 5] {
    let observation_count = histogram.iter().sum::<usize>();
    if observation_count == 0 {
        return [0.0; 5];
    }

    let percentile_fractions = [0.10, 0.25, 0.50, 0.75, 0.90];
    percentile_fractions.map(|fraction| {
        let percentile_index = (fraction * (observation_count - 1) as f64).floor() as usize;
        histogram_value_at_index(histogram, percentile_index) as f32
    })
}

fn histogram_value_at_index(histogram: &[usize], target_index: usize) -> usize {
    let mut cumulative_count = 0_usize;
    for (value, count) in histogram.iter().enumerate() {
        cumulative_count = cumulative_count.saturating_add(*count);
        if target_index < cumulative_count {
            return value;
        }
    }
    histogram.len().saturating_sub(1)
}

fn sequence_length_summary(
    read_lengths: &HashMap<usize, u64>,
) -> Result<SequenceLengthSummary, FastqAnalysisError> {
    let sorted_lengths: Vec<(usize, u64)> = read_lengths
        .iter()
        .map(|(&length, &count)| (length, count))
        .sorted_by_key(|(length, _)| *length)
        .collect();
    let sequence_count = sorted_lengths.iter().try_fold(0_u64, |total, (_, count)| {
        total
            .checked_add(*count)
            .ok_or(FastqAnalysisError::CountOverflow {
                resource: "sequence length count",
            })
    })?;
    let sum_len = sorted_lengths
        .iter()
        .try_fold(0_u64, |total, (length, count)| {
            let length_total =
                (*length as u64)
                    .checked_mul(*count)
                    .ok_or(FastqAnalysisError::CountOverflow {
                        resource: "sequence length sum",
                    })?;
            total
                .checked_add(length_total)
                .ok_or(FastqAnalysisError::CountOverflow {
                    resource: "sequence length sum",
                })
        })?;

    let lower_half_count = sequence_count.div_ceil(2);
    let upper_half_start = sequence_count / 2;

    let q1 = median_length_window(&sorted_lengths, 0, lower_half_count);
    let median = median_length_window(&sorted_lengths, 0, sequence_count);
    let q3 = median_length_window(&sorted_lengths, upper_half_start, lower_half_count);

    let n50_threshold = sum_len as f64 / 2.0;
    let mut cumulative_bases = 0_u64;
    let mut n50 = 0_usize;
    let mut n50_num = 0_u64;
    for (length, count) in sorted_lengths.iter().rev() {
        let length_total =
            (*length as u64)
                .checked_mul(*count)
                .ok_or(FastqAnalysisError::CountOverflow {
                    resource: "N50 length sum",
                })?;
        cumulative_bases = cumulative_bases.checked_add(length_total).ok_or(
            FastqAnalysisError::CountOverflow {
                resource: "N50 length sum",
            },
        )?;
        n50_num = n50_num
            .checked_add(1)
            .ok_or(FastqAnalysisError::CountOverflow {
                resource: "N50 length-bin count",
            })?;
        if cumulative_bases as f64 >= n50_threshold {
            n50 = *length;
            break;
        }
    }

    Ok(SequenceLengthSummary {
        sum_len,
        min_len: sorted_lengths.first().map_or(0, |(length, _)| *length),
        max_len: sorted_lengths.last().map_or(0, |(length, _)| *length),
        q1,
        median,
        q3,
        n50,
        n50_num,
    })
}

fn median_length_window(
    sorted_lengths: &[(usize, u64)],
    start_index: u64,
    window_count: u64,
) -> f64 {
    if window_count == 0 {
        return 0.0;
    }
    let middle_right = start_index + window_count / 2;
    if window_count.is_multiple_of(2) {
        let middle_left = middle_right - 1;
        (length_at_index(sorted_lengths, middle_left) as f64
            + length_at_index(sorted_lengths, middle_right) as f64)
            / 2.0
    } else {
        length_at_index(sorted_lengths, middle_right) as f64
    }
}

fn length_at_index(sorted_lengths: &[(usize, u64)], target_index: u64) -> usize {
    let mut cumulative_count = 0_u64;
    for (length, count) in sorted_lengths {
        cumulative_count = cumulative_count.saturating_add(*count);
        if target_index < cumulative_count {
            return *length;
        }
    }
    sorted_lengths.last().map_or(0, |(length, _)| *length)
}

/// Analyze a FASTQ file with the default resource limits.
pub fn analyze<P: AsRef<Path> + AsRef<OsStr>>(
    filename: P,
    k: u8,
) -> Result<FastqAnalysisResult, FastqAnalysisError> {
    let kmer_length = KmerLength::try_from(k)?;
    analyze_with_options(
        filename,
        AnalysisOptions::new(kmer_length, AnalysisLimits::default()),
    )
}

/// Analyze a FASTQ file with explicit, validated resource limits.
pub fn analyze_with_options<P: AsRef<Path> + AsRef<OsStr>>(
    filename: P,
    options: AnalysisOptions,
) -> Result<FastqAnalysisResult, FastqAnalysisError> {
    let input_path: &Path = filename.as_ref();
    let filename_string = input_path.to_string_lossy().into_owned();
    let mut reader = parse_fastx_file(&filename).map_err(|error| FastqAnalysisError::Input {
        path: filename_string.clone(),
        message: error.to_string(),
    })?;

    let mut mean_read_qualities = HashMap::default();
    let mut base_count = HashMap::default();
    let mut base_quality_count = HashMap::default();
    let mut base_totals = [0_u64; 5];
    let mut exact_n_bases = 0_u64;
    let mut gap_bases = 0_u64;
    let mut read_lengths = HashMap::default();
    let mut kmers = HashMap::default();
    let mut gc_content = vec![0_usize; 101];
    let mut quality_totals = QualityTotals::default();
    let mut total_length_bases = 0_u64;
    let mut read_count = 0_u64;
    let mut positions_truncated = false;

    while let Some(record) = reader.next() {
        let record_number = read_count
            .checked_add(1)
            .ok_or(FastqAnalysisError::CountOverflow { resource: "read" })?;
        let seqrec = record.map_err(|error| FastqAnalysisError::MalformedRecord {
            record_number,
            message: error.to_string(),
        })?;
        let sequence = seqrec.seq();
        let sequence_length = sequence.len();

        if sequence_length == 0 {
            return Err(FastqAnalysisError::EmptyRead { record_number });
        }
        if sequence_length > options.limits.max_read_length() {
            return Err(FastqAnalysisError::ResourceLimitExceeded {
                resource: "read length",
                limit: options.limits.max_read_length(),
                observed: sequence_length,
            });
        }
        if !read_lengths.contains_key(&sequence_length)
            && read_lengths.len() >= options.limits.max_distinct_read_lengths()
        {
            return Err(FastqAnalysisError::ResourceLimitExceeded {
                resource: "distinct read length",
                limit: options.limits.max_distinct_read_lengths(),
                observed: read_lengths.len() + 1,
            });
        }

        let qualities = seqrec
            .qual()
            .ok_or(FastqAnalysisError::MissingQualities { record_number })?;
        if qualities.len() != sequence_length {
            return Err(FastqAnalysisError::QualityLengthMismatch {
                record_number,
                sequence_length,
                quality_length: qualities.len(),
            });
        }

        let count_read_length = read_lengths.entry(sequence_length).or_insert(0_u64);
        *count_read_length =
            count_read_length
                .checked_add(1)
                .ok_or(FastqAnalysisError::CountOverflow {
                    resource: "read length",
                })?;
        total_length_bases = total_length_bases
            .checked_add(sequence_length as u64)
            .ok_or(FastqAnalysisError::CountOverflow {
                resource: "total sequence length",
            })?;
        gc_content[gc_percentage(&sequence)] = gc_content[gc_percentage(&sequence)]
            .checked_add(1)
            .ok_or(FastqAnalysisError::CountOverflow {
                resource: "GC histogram",
            })?;

        let mut read_quality_sum = 0_u64;
        for (position, &encoded_quality) in qualities.iter().enumerate() {
            let score =
                phred33_score(encoded_quality).map_err(|_| FastqAnalysisError::InvalidQuality {
                    record_number,
                    position: position + 1,
                    encoded_quality,
                })?;
            read_quality_sum = read_quality_sum.checked_add(score as u64).ok_or(
                FastqAnalysisError::CountOverflow {
                    resource: "read quality",
                },
            )?;
            quality_totals.total_bases = quality_totals.total_bases.checked_add(1).ok_or(
                FastqAnalysisError::CountOverflow {
                    resource: "quality base",
                },
            )?;
            quality_totals.score_sum = quality_totals.score_sum.checked_add(score as u64).ok_or(
                FastqAnalysisError::CountOverflow {
                    resource: "quality score",
                },
            )?;
            quality_totals.error_probability_sum += 10_f64.powf(-(score as f64) / 10.0);
            if score >= 20 {
                quality_totals.q20_bases = quality_totals.q20_bases.checked_add(1).ok_or(
                    FastqAnalysisError::CountOverflow {
                        resource: "Q20 base",
                    },
                )?;
            }
            if score >= 30 {
                quality_totals.q30_bases = quality_totals.q30_bases.checked_add(1).ok_or(
                    FastqAnalysisError::CountOverflow {
                        resource: "Q30 base",
                    },
                )?;
            }

            if position < options.limits.max_tracked_positions() {
                let histogram = base_quality_count
                    .entry(position)
                    .or_insert_with(|| vec![0_usize; 94]);
                histogram[score] =
                    histogram[score]
                        .checked_add(1)
                        .ok_or(FastqAnalysisError::CountOverflow {
                            resource: "per-position quality",
                        })?;
            }
        }

        let mean_read_quality = read_quality_sum / sequence_length as u64;
        let mean_quality_count = mean_read_qualities
            .entry(mean_read_quality)
            .or_insert(0_u64);
        *mean_quality_count =
            mean_quality_count
                .checked_add(1)
                .ok_or(FastqAnalysisError::CountOverflow {
                    resource: "mean read quality",
                })?;

        for (position, &base) in sequence.iter().enumerate() {
            if matches!(base, b'N' | b'n') {
                exact_n_bases = exact_n_bases
                    .checked_add(1)
                    .ok_or(FastqAnalysisError::CountOverflow { resource: "N base" })?;
            }
            if matches!(base, b'-' | b' ' | b'.') {
                gap_bases = gap_bases
                    .checked_add(1)
                    .ok_or(FastqAnalysisError::CountOverflow {
                        resource: "gap base",
                    })?;
            }
            let index = base_index(base);
            base_totals[index] = base_totals[index]
                .checked_add(1)
                .ok_or(FastqAnalysisError::CountOverflow { resource: "base" })?;
            if position < options.limits.max_tracked_positions() {
                let position_count = base_count.entry(position).or_insert_with(|| vec![0_u64; 5]);
                position_count[index] = position_count[index].checked_add(1).ok_or(
                    FastqAnalysisError::CountOverflow {
                        resource: "per-position base",
                    },
                )?;
            }
        }
        positions_truncated |= sequence_length > options.limits.max_tracked_positions();

        let normalized_sequence = seqrec.normalize(false);
        let reverse_complement = normalized_sequence.reverse_complement();
        for (_, kmer, _) in
            normalized_sequence.canonical_kmers(options.kmer_length.get(), &reverse_complement)
        {
            let count = kmers.entry(kmer.to_owned()).or_insert(0_u64);
            *count = count
                .checked_add(1)
                .ok_or(FastqAnalysisError::CountOverflow { resource: "k-mer" })?;
        }

        read_count = record_number;
    }

    if read_count == 0 {
        return Err(FastqAnalysisError::EmptyInput);
    }

    let called_bases = base_totals[A]
        .checked_add(base_totals[C])
        .and_then(|count| count.checked_add(base_totals[G]))
        .and_then(|count| count.checked_add(base_totals[T]))
        .ok_or(FastqAnalysisError::CountOverflow {
            resource: "called base",
        })?;
    let avg_gc = if called_bases == 0 {
        0.0
    } else {
        (base_totals[G] + base_totals[C]) as f64 / called_bases as f64 * 100.0
    };

    Ok(FastqAnalysisResult {
        filename: filename_string,
        read_count,
        avg_read_length: total_length_bases as f64 / read_count as f64,
        avg_gc,
        read_lengths,
        base_quality_count,
        base_count,
        base_totals,
        exact_n_bases,
        gap_bases,
        gc_distribution: gc_content,
        kmer_counts: kmers,
        mean_read_qualities,
        quality_totals,
        positions_truncated,
        kmer_length: options.kmer_length,
    })
}

/// CLI 入口函数
pub fn process<P: AsRef<Path> + AsRef<OsStr>>(
    filename: P,
    k: u8,
    summary: Option<P>,
    skip_html: bool,
) -> Result<(), Box<dyn Error>> {
    let analysis = analyze(&filename, k)?;
    let read_count = analysis.read_count;
    let avg_gc = analysis.avg_gc;
    let avg_read_length = analysis.avg_read_length;
    let read_lengths = &analysis.read_lengths;
    let base_quality_count = &analysis.base_quality_count;
    let base_count = &analysis.base_count;
    let gc_content = &analysis.gc_distribution;
    let kmers = &analysis.kmer_counts;
    let mean_read_qualities = &analysis.mean_read_qualities;
    let quality_totals = analysis.quality_totals;
    let base_totals = analysis.base_totals;
    let exact_n_bases = analysis.exact_n_bases;
    let gap_bases = analysis.gap_bases;
    let positions_truncated = analysis.positions_truncated;
    let k = analysis.kmer_length.get();

    // Data for base per position
    let mut n_warn = "pass";
    let mut base_count_data = Vec::new();
    let mut base_count_percentage = Vec::new();
    let mut gc_content_per_base = Vec::new();
    let mut base_warning = "pass";
    for (position, bases) in base_count.iter().sorted_by_key(|(position, _)| *position) {
        let tmp_sum = bases.iter().sum::<u64>();
        let called_base_count = tmp_sum - bases.get(N).unwrap();
        let tmp_gc = if called_base_count == 0 {
            0.0
        } else {
            (bases.get(G).unwrap() + bases.get(C).unwrap()) as f64 / called_base_count as f64
        };
        let display_position = position + 1;
        gc_content_per_base.push(json!({
            "pos": display_position,
            "pct": tmp_gc * 100.0,
        }));
        let percentages = bases
            .iter()
            .enumerate()
            .map(|(base_index, &count)| {
                let base = BASES[base_index];
                let denominator = if base == 'N' {
                    tmp_sum
                } else {
                    called_base_count
                };
                let percentage = if denominator == 0 {
                    0.0
                } else {
                    count as f64 / denominator as f64
                };
                if base == 'N' && percentage >= 0.20 {
                    n_warn = "fail"
                } else if base == 'N' && percentage >= 0.05 && n_warn != "fail" {
                    n_warn = "warn"
                };
                (base, percentage * 100.0)
            })
            .collect::<HashMap<char, f64>>();
        base_count_percentage.push(json!({
            "pos": display_position,
            "base_map": percentages,
        }));
        let gc_diff = i64::abs(*bases.get(G).unwrap() as i64 - *bases.get(C).unwrap() as i64)
            as f64
            / tmp_sum as f64;
        let at_diff = i64::abs(*bases.get(A).unwrap() as i64 - *bases.get(T).unwrap() as i64)
            as f64
            / tmp_sum as f64;
        if gc_diff >= 0.20 || at_diff >= 0.20 {
            base_warning = "fail"
        } else if (gc_diff >= 0.10 || at_diff >= 0.10) && base_warning != "fail" {
            base_warning = "warn"
        }
        for (base, &count) in bases.iter().enumerate() {
            base_count_data.push(json!({"pos": position + 1, "count": count, "base": BASES[base]}));
        }
    }

    let mut bpp_specs: Value =
        serde_json::from_str(include_str!("report/base_per_pos_specs.json"))?;
    bpp_specs["data"]["values"] = json!(base_count_data);

    // Data for read lengths
    let read_length_warn = if read_lengths.len() > 1 {
        "warn"
    } else {
        "pass"
    };
    let mut read_length_data = Vec::new();
    for (length, count) in read_lengths.iter().sorted_by_key(|(length, _)| *length) {
        read_length_data.push(json!({"length": length, "count": count}));
    }

    let mut rle_specs: Value =
        serde_json::from_str(include_str!("report/read_lengths_specs.json"))?;
    rle_specs["data"]["values"] = json!(read_length_data);

    // Data for mean read qualities
    let mut mean_read_quality_data = Vec::new();
    for (quality, count) in mean_read_qualities
        .iter()
        .sorted_by_key(|(quality, _)| *quality)
    {
        mean_read_quality_data.push(json!({"score": quality, "count": count}))
    }

    let mut sqc_specs: Value =
        serde_json::from_str(include_str!("report/sequence_quality_score_specs.json"))?;
    sqc_specs["data"]["values"] = json!(mean_read_quality_data);

    // Data for kmer quantities
    let mut kmer_data = Vec::new();
    for (kmer, count) in kmers {
        kmer_data.push(json!({"k_mer": std::str::from_utf8(kmer).unwrap(), "count": count}))
    }

    let mut overly_represented = Vec::new();
    let mut overly_represented_warn = "pass";
    let total_kmers = kmers.values().sum::<u64>();
    for (km, occ) in kmers
        .iter()
        .sorted_by(|(_, a), (_, b)| Ord::cmp(&b, &a))
        .take(5)
    {
        let percentage = *occ as f64 / total_kmers as f64;
        if percentage >= 0.01 {
            overly_represented_warn = "fail";
            overly_represented.push(json!({"k_mer": std::str::from_utf8(km).unwrap(), "count": occ, "pct": percentage, "or": "Yes"}));
        } else if percentage >= 0.002 {
            if overly_represented_warn != "fail" {
                overly_represented_warn = "warn";
            }
            overly_represented.push(json!({"k_mer": std::str::from_utf8(km).unwrap(), "count": occ, "pct": percentage, "or": "No"}));
        };
    }

    let mut counter_specs: Value = serde_json::from_str(include_str!("report/counter_specs.json"))?;
    counter_specs["data"]["values"] = json!(kmer_data);

    // Data for GC content
    let mut gc_data = Vec::new();
    // TODO(lhepler): I had to add this filter to get the same histogram in the output report
    // but I think it may be more correct to leave it out. Arguable, though
    for (perc, &count) in gc_content
        .iter()
        .enumerate()
        .filter(|(_, &count)| count > 0)
    {
        gc_data.push(json!({"gc_pct": perc, "count": count, "type": "gc"}))
    }

    let mut gc_specs: Value = serde_json::from_str(include_str!("report/gc_content_specs.json"))?;
    gc_specs["data"]["values"] = json!(gc_data);

    // Data for base quality per position
    let mut base_quality_warn = "pass";
    let mut base_per_pos_data = Vec::new();
    for (position, qualities) in base_quality_count
        .iter()
        .sorted_by_key(|(position, _)| *position)
    {
        let (sum, len) = qualities
            .iter()
            .enumerate()
            .fold((0_usize, 0_usize), |(s, l), (q, c)| (s + q * c, l + c));
        let avg = sum as f64 / len as f64;
        let values = quality_percentiles(qualities);
        if values.get(2).unwrap() <= &20_f32 {
            base_quality_warn = "fail"
        } else if values.get(2).unwrap() <= &25_f32 && base_quality_warn != "fail" {
            base_quality_warn = "warn"
        }
        base_per_pos_data.push(json!({
        "pos": position + 1,
        "average": avg,
        "lower": values.first().unwrap(),
        "q1": values.get(1).unwrap(),
        "median":values.get(2).unwrap(),
        "q3": values.get(3).unwrap(),
        "upper": values.get(4).unwrap(),
        }));
    }

    let mut qpp_specs: Value =
        serde_json::from_str(include_str!("report/quality_per_pos_specs.json"))?;
    qpp_specs["data"]["values"] = json!(base_per_pos_data);

    let plots = json!({
        "k-mer quantities": {"short": "count", "specs": counter_specs.to_string()},
        "gc content": {"short": "gc", "specs": gc_specs.to_string()},
        "base sequence quality": {"short": "base", "specs": qpp_specs.to_string()},
        "sequence quality score": {"short": "qual", "specs": sqc_specs.to_string()},
        "base sequence content": {"short": "cont", "specs": bpp_specs.to_string()},
        "read lengths": {"short": "rlen", "specs": rle_specs.to_string()},
    });

    let file = Path::new(&filename)
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| analysis.filename.clone());
    let meta = json!({
        "file name": {"name": "file name", "value": file},
        "canonical": {"name": "canonical", "value": "True"},
        "k": {"name": "k", "value": k},
        "total reads": {"name": "total reads", "value": read_count},
        "average GC content": {"name": "average GC content", "value": avg_gc},
        "average read length": {"name": "average read length", "value": avg_read_length},
        "position statistics truncated": {"name": "position statistics truncated", "value": positions_truncated},
    });

    let mut templates = Tera::default();
    let mut context = Context::new();

    // Always set these variables needed by summary template
    let local: DateTime<Local> = Local::now();
    context.insert("time", &local.format("%a %b %e %T %Y").to_string());
    context.insert("version", &env!("CARGO_PKG_VERSION"));
    context.insert("positions_truncated", &positions_truncated);

    if !skip_html {
        templates.add_raw_template("report.html.tera", include_str!("report/report.html.tera"))?;
        context.insert("plots", &plots);
        context.insert("meta", &meta);
        let html = templates.render("report.html.tera", &context)?;
        io::stdout().write_all(html.as_bytes())?;
    }

    if let Some(path) = summary {
        let output_path = Path::new(&path);
        std::fs::create_dir_all(output_path)?;

        let length_summary = sequence_length_summary(read_lengths)?;
        let sum_n = exact_n_bases;
        let seqkit_gc = if length_summary.sum_len == 0 {
            0.0
        } else {
            (base_totals[G] + base_totals[C]) as f64 / length_summary.sum_len as f64 * 100.0
        };
        let q20_pct = quality_totals.q20_percent();
        let q30_pct = quality_totals.q30_percent();
        let avg_qual = quality_totals.average_quality();

        templates.add_raw_template(
            "fastqc_summary.txt.tera",
            include_str!("report/fastqc_summary.txt.tera"),
        )?;
        context.insert("filename", &file);
        context.insert("reads", &read_count);
        context.insert("avg_read_length", &avg_read_length);
        context.insert("avg_gc", &avg_gc);
        context.insert("bpp_data", &base_per_pos_data);
        context.insert("mean_read_quality_data", &mean_read_quality_data);
        context.insert("base_count", &base_count_percentage);
        context.insert("base_warning", &base_warning);
        context.insert("read_lengths", &read_length_data);
        context.insert("read_length_warn", &read_length_warn);
        context.insert("gc_data", &gc_data);
        context.insert("gc_per_base", &gc_content_per_base);
        context.insert("overly_represented", &overly_represented);
        context.insert("overly_represented_warn", &overly_represented_warn);
        context.insert("n_warn", &n_warn);
        context.insert("base_quality_warn", &base_quality_warn);

        // 添加 seqkit 风格的统计指标
        context.insert("seqkit_format", &"FASTQ".to_string());
        context.insert("seqkit_type", &"DNA".to_string());
        context.insert("sum_len", &length_summary.sum_len);
        context.insert("min_len", &length_summary.min_len);
        context.insert("max_len", &length_summary.max_len);
        context.insert("q1", &length_summary.q1);
        context.insert("q2", &length_summary.median);
        context.insert("q3", &length_summary.q3);
        context.insert("sum_gap", &gap_bases);
        context.insert("n50", &length_summary.n50);
        context.insert("n50_num", &length_summary.n50_num);
        context.insert("q20_pct", &q20_pct);
        context.insert("q30_pct", &q30_pct);
        context.insert("avg_qual", &avg_qual);
        context.insert("seqkit_gc", &seqkit_gc);
        context.insert("sum_n", &sum_n);

        let txt = templates.render("fastqc_summary.txt.tera", &context)?;
        let mut file = File::create(output_path.join("fastqc_data.txt"))?;
        file.write_all(txt.as_bytes())?;
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use super::{
        analyze, analyze_with_options, gc_percentage, phred33_score, process, quality_percentiles,
        sequence_length_summary,
    };
    use crate::types::{
        AnalysisLimits, AnalysisOptions, FastqAnalysisError, KmerLength, MAX_KMER_LENGTH,
    };
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use rustc_hash::FxHashMap as HashMap;
    use std::convert::TryFrom;
    use std::fs;
    use std::io::Write;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::time::{SystemTime, UNIX_EPOCH};

    static NEXT_TEMPORARY_PATH: AtomicU64 = AtomicU64::new(0);

    fn temporary_path(extension: &str) -> std::path::PathBuf {
        let unique_suffix = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let sequence = NEXT_TEMPORARY_PATH.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "fastqcx-{}-{unique_suffix}-{sequence}.{extension}",
            std::process::id()
        ))
    }

    fn temporary_fastq(contents: &str) -> std::path::PathBuf {
        let path = temporary_path("fastq");
        fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn test_gc_percentage_excludes_n_bases() {
        assert_eq!(gc_percentage(b"GCNN"), 100);
        assert_eq!(gc_percentage(b"ATNN"), 0);
        assert_eq!(gc_percentage(b"NNNN"), 0);
    }

    #[test]
    fn test_analyze_all_n_fastq() {
        let path = temporary_fastq("@all_n\nNNNN\n+\nIIII\n");
        let result = analyze(&path, 2).unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(result.read_count, 1);
        assert_eq!(result.avg_gc, 0.0);
        assert_eq!(result.gc_distribution[0], 1);
    }

    #[test]
    fn test_analyze_empty_fastq_returns_error() {
        let path = temporary_fastq("");
        assert!(analyze(&path, 2).is_err());
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_analyze_rejects_malformed_record_after_valid_data() {
        let path = temporary_fastq("@valid\nAC\n+\nII\n@broken\nAC\n+\nI\n");
        let error = analyze(&path, 2).unwrap_err();
        fs::remove_file(path).unwrap();

        assert!(matches!(
            error,
            FastqAnalysisError::MalformedRecord {
                record_number: 2,
                ..
            }
        ));
    }

    #[test]
    fn test_analyze_rejects_truncated_gzip() {
        let path = temporary_path("fastq.gz");
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder
            .write_all(b"@valid\nACGT\n+\nIIII\n@broken\nACGT\n+\nIIII\n")
            .unwrap();
        let mut compressed = encoder.finish().unwrap();
        compressed.truncate(compressed.len() / 2);
        fs::write(&path, compressed).unwrap();

        assert!(analyze(&path, 2).is_err());
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_analyze_rejects_fasta_input() {
        let path = temporary_fastq(">not_fastq\nACGT\n");
        let error = analyze(&path, 2).unwrap_err();
        fs::remove_file(path).unwrap();

        assert!(matches!(
            error,
            FastqAnalysisError::MissingQualities { record_number: 1 }
        ));
    }

    #[test]
    fn test_analyze_rejects_empty_read() {
        let path = temporary_fastq("@empty\n\n+\n\n");
        assert!(analyze(&path, 2).is_err());
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_analyze_rejects_out_of_range_kmer_lengths() {
        let path = temporary_fastq("@read\nACGT\n+\nIIII\n");
        assert!(matches!(
            analyze(&path, 0),
            Err(FastqAnalysisError::InvalidKmerLength { requested: 0 })
        ));
        assert!(matches!(
            analyze(&path, MAX_KMER_LENGTH + 1),
            Err(FastqAnalysisError::InvalidKmerLength { .. })
        ));
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn test_position_limit_bounds_histograms_but_preserves_global_totals() {
        let sequence = "ACGTACGTACGT";
        let qualities = "!!!!IIIIIIII";
        let path = temporary_fastq(&format!("@long\n{sequence}\n+\n{qualities}\n"));
        let limits = AnalysisLimits::new(4, 100, 10).unwrap();
        let options = AnalysisOptions::new(KmerLength::try_from(2).unwrap(), limits);
        let result = analyze_with_options(&path, options).unwrap();
        fs::remove_file(path).unwrap();

        assert!(result.positions_truncated);
        assert_eq!(result.base_count.len(), 4);
        assert_eq!(result.base_quality_count.len(), 4);
        assert_eq!(result.base_totals.iter().sum::<u64>(), 12);
        assert_eq!(result.quality_totals.total_bases, 12);
        assert_eq!(result.quality_totals.q30_percent(), 8.0 / 12.0 * 100.0);
        assert_eq!(
            result.quality_totals.arithmetic_average_quality(),
            80.0 / 3.0
        );
    }

    #[test]
    fn default_position_budget_is_machine_enforced() {
        let sequence = "A".repeat(super::AnalysisLimits::default().max_tracked_positions() + 1);
        let qualities = "I".repeat(sequence.len());
        let path = temporary_fastq(&format!("@bounded\n{sequence}\n+\n{qualities}\n"));
        let result = analyze(&path, MAX_KMER_LENGTH).unwrap();
        fs::remove_file(path).unwrap();

        assert!(result.positions_truncated);
        assert_eq!(
            result.base_count.len(),
            super::AnalysisLimits::default().max_tracked_positions()
        );
        assert_eq!(
            result.base_quality_count.len(),
            super::AnalysisLimits::default().max_tracked_positions()
        );
        assert!(result.kmer_counts.len() <= 4_usize.pow(MAX_KMER_LENGTH as u32));
    }

    #[test]
    fn test_read_length_limit_fails_closed() {
        let path = temporary_fastq("@long\nACGTA\n+\nIIIII\n");
        let limits = AnalysisLimits::new(4, 4, 10).unwrap();
        let options = AnalysisOptions::new(KmerLength::try_from(2).unwrap(), limits);
        let error = analyze_with_options(&path, options).unwrap_err();
        fs::remove_file(path).unwrap();

        assert!(matches!(
            error,
            FastqAnalysisError::ResourceLimitExceeded {
                resource: "read length",
                limit: 4,
                observed: 5
            }
        ));
    }

    #[test]
    fn test_distinct_read_length_limit_fails_closed() {
        let path = temporary_fastq("@one\nA\n+\nI\n@two\nAC\n+\nII\n@three\nACG\n+\nIII\n");
        let limits = AnalysisLimits::new(10, 10, 2).unwrap();
        let options = AnalysisOptions::new(KmerLength::try_from(1).unwrap(), limits);
        let error = analyze_with_options(&path, options).unwrap_err();
        fs::remove_file(path).unwrap();

        assert!(matches!(
            error,
            FastqAnalysisError::ResourceLimitExceeded {
                resource: "distinct read length",
                limit: 2,
                observed: 3
            }
        ));
    }

    #[test]
    fn test_global_gc_is_weighted_by_called_bases() {
        let path = temporary_fastq("@short\nG\n+\nI\n@long\nAAAAAAAAA\n+\nIIIIIIIII\n");
        let result = analyze(&path, 1).unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(result.avg_gc, 10.0);
    }

    #[test]
    fn test_length_summary_does_not_expand_per_read_values() {
        let mut lengths = HashMap::default();
        lengths.insert(100, 1_000_000_000);
        lengths.insert(200, 1);
        let summary = sequence_length_summary(&lengths).unwrap();

        assert_eq!(summary.sum_len, 100_000_000_200);
        assert_eq!(summary.min_len, 100);
        assert_eq!(summary.max_len, 200);
        assert_eq!(summary.n50, 100);
    }

    #[test]
    fn test_crlf_fastq_is_accepted() {
        let path = temporary_fastq("@read\r\nACGT\r\n+\r\nIIII\r\n");
        let result = analyze(&path, 2).unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(result.read_count, 1);
        assert_eq!(result.quality_totals.total_bases, 4);
    }

    #[test]
    fn test_phred33_boundaries() {
        assert_eq!(phred33_score(b'!').unwrap(), 0);
        assert_eq!(phred33_score(b'~').unwrap(), 93);
        assert!(phred33_score(32).is_err());
        assert!(phred33_score(127).is_err());
    }
    #[test]
    fn test_quality_percentiles_cover_fastqc_columns() {
        let histogram = vec![1; 100];
        assert_eq!(
            quality_percentiles(&histogram),
            [9.0, 24.0, 49.0, 74.0, 89.0]
        );
    }

    #[test]
    fn test_quality_percentiles_handle_repeated_scores() {
        let mut histogram = [0_usize; 76];
        histogram[25] = 75;
        histogram[75] = 25;
        assert_eq!(
            quality_percentiles(&histogram),
            [25.0, 25.0, 25.0, 25.0, 75.0]
        );
    }

    #[test]
    fn test_quality_percentiles_empty_histogram_is_zeroed() {
        assert_eq!(quality_percentiles(&[]), [0.0; 5]);
        assert_eq!(quality_percentiles(&[0; 10]), [0.0; 5]);
    }

    #[test]
    fn test_analyze_treats_iupac_bases_as_n() {
        let path = temporary_fastq("@iupac\nACGTRY\n+\nIIIIII\n");
        let result = analyze(&path, 2).unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(result.base_count[&4][4], 1);
        assert_eq!(result.base_count[&5][4], 1);
        assert_eq!(result.avg_gc, 50.0);
    }

    #[test]
    fn test_summary_quality_statistics_are_base_weighted() {
        let path = temporary_fastq("@read1\nAC\n+\n5?\n@read2\nGT\n+\n+I\n");
        let output_dir = path.with_extension("summary");
        process(&path, 2, Some(&output_dir), true).unwrap();

        let summary = fs::read_to_string(output_dir.join("fastqc_data.txt")).unwrap();
        let filename = path.file_name().unwrap().to_string_lossy();
        let statistics_line = summary
            .lines()
            .find(|line| line.starts_with(filename.as_ref()))
            .unwrap()
            .trim_end_matches(">>END_MODULE");
        let fields: Vec<&str> = statistics_line.split('\t').collect();

        assert_eq!(fields[14].parse::<f64>().unwrap(), 75.0);
        assert_eq!(fields[15].parse::<f64>().unwrap(), 50.0);
        let expected_avg_quality = -10.0_f64
            * ((10.0_f64.powf(-20.0 / 10.0)
                + 10.0_f64.powf(-30.0 / 10.0)
                + 10.0_f64.powf(-10.0 / 10.0)
                + 10.0_f64.powf(-40.0 / 10.0))
                / 4.0)
                .log10();
        assert!((fields[16].parse::<f64>().unwrap() - expected_avg_quality).abs() < 0.01);

        fs::remove_dir_all(output_dir).unwrap();
        fs::remove_file(path).unwrap();
    }
}
