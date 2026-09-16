//! fastqcx - A fast quality control tool for FASTQ files
//!
//! This crate provides both a library API and a command-line tool for
//! analyzing FASTQ sequencing files.
//!
//! # Library Usage
//!
//! ```no_run
//! use fastqcx::analyze;
//!
//! // Analyze a FASTQ file
//! let result = analyze("sample.fastq", 5).unwrap();
//! println!("Total reads: {}", result.read_count);
//! println!("Average GC: {}%", result.avg_gc);
//! ```

pub mod process;
pub mod types;

pub use process::{analyze, analyze_with_options};
pub use types::{
    AnalysisLimits, AnalysisOptions, FastqAnalysisError, FastqAnalysisResult, KmerLength,
    QualityTotals, SequenceLengthSummary,
};
