use rustc_hash::FxHashMap;
use std::convert::TryFrom;
use std::error::Error;
use std::fmt;

pub const MIN_KMER_LENGTH: u8 = 1;
pub const MAX_KMER_LENGTH: u8 = 7;
pub const DEFAULT_MAX_TRACKED_POSITIONS: usize = 10_000;
pub const DEFAULT_MAX_READ_LENGTH: usize = 10_000_000;
pub const DEFAULT_MAX_DISTINCT_READ_LENGTHS: usize = 100_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KmerLength(u8);

impl KmerLength {
    pub const fn get(self) -> u8 {
        self.0
    }
}

impl TryFrom<u8> for KmerLength {
    type Error = FastqAnalysisError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        if (MIN_KMER_LENGTH..=MAX_KMER_LENGTH).contains(&value) {
            Ok(Self(value))
        } else {
            Err(FastqAnalysisError::InvalidKmerLength { requested: value })
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AnalysisLimits {
    max_tracked_positions: usize,
    max_read_length: usize,
    max_distinct_read_lengths: usize,
}

impl AnalysisLimits {
    pub fn new(
        max_tracked_positions: usize,
        max_read_length: usize,
        max_distinct_read_lengths: usize,
    ) -> Result<Self, FastqAnalysisError> {
        for (name, value) in [
            ("max_tracked_positions", max_tracked_positions),
            ("max_read_length", max_read_length),
            ("max_distinct_read_lengths", max_distinct_read_lengths),
        ] {
            if value == 0 {
                return Err(FastqAnalysisError::InvalidLimit { name, value });
            }
        }

        Ok(Self {
            max_tracked_positions,
            max_read_length,
            max_distinct_read_lengths,
        })
    }

    pub const fn max_tracked_positions(self) -> usize {
        self.max_tracked_positions
    }

    pub const fn max_read_length(self) -> usize {
        self.max_read_length
    }

    pub const fn max_distinct_read_lengths(self) -> usize {
        self.max_distinct_read_lengths
    }
}

impl Default for AnalysisLimits {
    fn default() -> Self {
        Self {
            max_tracked_positions: DEFAULT_MAX_TRACKED_POSITIONS,
            max_read_length: DEFAULT_MAX_READ_LENGTH,
            max_distinct_read_lengths: DEFAULT_MAX_DISTINCT_READ_LENGTHS,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AnalysisOptions {
    pub kmer_length: KmerLength,
    pub limits: AnalysisLimits,
}

impl AnalysisOptions {
    pub const fn new(kmer_length: KmerLength, limits: AnalysisLimits) -> Self {
        Self {
            kmer_length,
            limits,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct QualityTotals {
    pub total_bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    pub score_sum: u64,
    pub error_probability_sum: f64,
}

impl QualityTotals {
    pub fn q20_percent(self) -> f64 {
        percentage(self.q20_bases, self.total_bases)
    }

    pub fn q30_percent(self) -> f64 {
        percentage(self.q30_bases, self.total_bases)
    }

    pub fn average_quality(self) -> f64 {
        if self.total_bases == 0 || self.error_probability_sum == 0.0 {
            0.0
        } else {
            -10.0 * (self.error_probability_sum / self.total_bases as f64).log10()
        }
    }

    pub fn arithmetic_average_quality(self) -> f64 {
        if self.total_bases == 0 {
            0.0
        } else {
            self.score_sum as f64 / self.total_bases as f64
        }
    }
}

fn percentage(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64 * 100.0
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SequenceLengthSummary {
    pub sum_len: u64,
    pub min_len: usize,
    pub max_len: usize,
    pub q1: f64,
    pub median: f64,
    pub q3: f64,
    pub n50: usize,
    pub n50_num: u64,
}

#[derive(Debug, Clone)]
pub struct FastqAnalysisResult {
    pub filename: String,
    pub read_count: u64,
    pub avg_read_length: f64,
    pub avg_gc: f64,
    pub read_lengths: FxHashMap<usize, u64>,
    pub base_quality_count: FxHashMap<usize, Vec<usize>>,
    pub base_count: FxHashMap<usize, Vec<u64>>,
    pub base_totals: [u64; 5],
    pub exact_n_bases: u64,
    pub gap_bases: u64,
    pub gc_distribution: Vec<usize>,
    pub kmer_counts: FxHashMap<Vec<u8>, u64>,
    pub mean_read_qualities: FxHashMap<u64, u64>,
    pub quality_totals: QualityTotals,
    pub positions_truncated: bool,
    pub kmer_length: KmerLength,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FastqAnalysisError {
    InvalidKmerLength {
        requested: u8,
    },
    InvalidLimit {
        name: &'static str,
        value: usize,
    },
    Input {
        path: String,
        message: String,
    },
    MalformedRecord {
        record_number: u64,
        message: String,
    },
    EmptyRead {
        record_number: u64,
    },
    MissingQualities {
        record_number: u64,
    },
    QualityLengthMismatch {
        record_number: u64,
        sequence_length: usize,
        quality_length: usize,
    },
    InvalidQuality {
        record_number: u64,
        position: usize,
        encoded_quality: u8,
    },
    ResourceLimitExceeded {
        resource: &'static str,
        limit: usize,
        observed: usize,
    },
    CountOverflow {
        resource: &'static str,
    },
    EmptyInput,
}

impl fmt::Display for FastqAnalysisError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKmerLength { requested } => write!(
                formatter,
                "invalid k-mer length {requested}; expected {MIN_KMER_LENGTH}..={MAX_KMER_LENGTH}"
            ),
            Self::InvalidLimit { name, value } => {
                write!(formatter, "invalid analysis limit {name}={value}; expected a positive value")
            }
            Self::Input { path, message } => {
                write!(formatter, "failed to open FASTQ {path}: {message}")
            }
            Self::MalformedRecord {
                record_number,
                message,
            } => write!(
                formatter,
                "malformed FASTQ record {record_number}: {message}"
            ),
            Self::EmptyRead { record_number } => {
                write!(formatter, "FASTQ record {record_number} has an empty sequence")
            }
            Self::MissingQualities { record_number } => write!(
                formatter,
                "record {record_number} has no quality string; FASTQ input is required"
            ),
            Self::QualityLengthMismatch {
                record_number,
                sequence_length,
                quality_length,
            } => write!(
                formatter,
                "FASTQ record {record_number} has sequence length {sequence_length} but quality length {quality_length}"
            ),
            Self::InvalidQuality {
                record_number,
                position,
                encoded_quality,
            } => write!(
                formatter,
                "FASTQ record {record_number} has invalid Phred+33 byte {encoded_quality} at position {position}"
            ),
            Self::ResourceLimitExceeded {
                resource,
                limit,
                observed,
            } => write!(
                formatter,
                "FASTQ {resource} limit exceeded: observed {observed}, limit {limit}"
            ),
            Self::CountOverflow { resource } => {
                write!(formatter, "FASTQ {resource} counter overflowed")
            }
            Self::EmptyInput => write!(formatter, "FASTQ contains no reads"),
        }
    }
}

impl Error for FastqAnalysisError {}

#[cfg(test)]
mod tests {
    use super::{AnalysisLimits, FastqAnalysisError, KmerLength, MAX_KMER_LENGTH};
    use std::convert::TryFrom;

    #[test]
    fn kmer_length_accepts_only_bounded_values() {
        assert_eq!(KmerLength::try_from(1).unwrap().get(), 1);
        assert_eq!(
            KmerLength::try_from(MAX_KMER_LENGTH).unwrap().get(),
            MAX_KMER_LENGTH
        );
        assert!(matches!(
            KmerLength::try_from(0),
            Err(FastqAnalysisError::InvalidKmerLength { requested: 0 })
        ));
        assert!(KmerLength::try_from(MAX_KMER_LENGTH + 1).is_err());
    }

    #[test]
    fn analysis_limits_must_be_positive() {
        assert!(AnalysisLimits::new(0, 1, 1).is_err());
        assert!(AnalysisLimits::new(1, 0, 1).is_err());
        assert!(AnalysisLimits::new(1, 1, 0).is_err());
        assert!(AnalysisLimits::new(1, 1, 1).is_ok());
    }
}
