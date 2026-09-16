use clap::{Arg, Command};
use fastqcx::KmerLength;
use std::convert::TryFrom;
use std::error::Error;

fn parse_kmer_length(value: &str) -> Result<KmerLength, String> {
    let parsed = value
        .parse::<u8>()
        .map_err(|error| format!("invalid k-mer length {value:?}: {error}"))?;
    KmerLength::try_from(parsed).map_err(|error| error.to_string())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("fastqcx")
        .about("A FASTQ quality control tool inspired by FastQC")
        .version(env!("CARGO_PKG_VERSION"))
        .author("Felix W. <fxwiegand@wgdnet.de>")
        .arg(
            Arg::new("fastq")
                .short('q')
                .long("fastq")
                .value_name("FILE")
                .help("The input FASTQ file to use.")
                .required(true)
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            Arg::new("k")
                .short('k')
                .long("kmer")
                .value_name("K")
                .help("K-mer length for counting (supported range: 1-7).")
                .default_value("5")
                .value_parser(parse_kmer_length),
        )
        .arg(
            Arg::new("summary")
                .short('s')
                .long("summary")
                .value_name("DIRECTORY")
                .required(false)
                .help("Create DIRECTORY/fastqc_data.txt for qctb and MultiQC consumers.")
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            Arg::new("no_html")
                .long("no-html")
                .required(false)
                .default_value("false")
                .help("Skip HTML report output to stdout.")
                .action(clap::ArgAction::SetTrue),
        )
        .get_matches();

    let fastq_file = matches
        .get_one::<String>("fastq")
        .expect("required by clap");
    let kmer_length = matches
        .get_one::<KmerLength>("k")
        .copied()
        .expect("defaulted by clap");
    let summary = matches.get_one::<String>("summary");
    let skip_html = matches.get_one::<bool>("no_html").copied().unwrap_or(false);

    fastqcx::process::process(fastq_file, kmer_length.get(), summary, skip_html)
}
