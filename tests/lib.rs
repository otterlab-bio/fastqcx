use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn temporary_path(name: &str) -> PathBuf {
    let unique_suffix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock should be after UNIX epoch")
        .as_nanos();
    std::env::temp_dir().join(format!("fastqcx-integration-{unique_suffix}-{name}"))
}

fn run_fastqcx(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fastqcx"))
        .args(arguments)
        .env("HTTPS_PROXY", "http://127.0.0.1:9")
        .env("HTTP_PROXY", "http://127.0.0.1:9")
        .env("NO_PROXY", "")
        .output()
        .expect("fastqcx should start")
}

fn extract_module_positions(summary: &str, module_name: &str) -> Vec<usize> {
    let module_header = format!(">>{module_name}\t");
    let mut in_module = false;
    let mut positions = Vec::new();

    for line in summary.lines() {
        if line.starts_with(&module_header) {
            in_module = true;
            continue;
        }
        if !in_module {
            continue;
        }
        if line.starts_with(">>END_MODULE") {
            break;
        }
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let position = line
            .split('\t')
            .next()
            .expect("module row should contain a position")
            .parse::<usize>()
            .expect("module position should be an unsigned integer");
        positions.push(position);
    }

    positions
}

#[test]
fn html_report_generation_does_not_fetch_remote_assets() {
    let output = run_fastqcx(&["--fastq", "tests/resources/example.fastq"]);

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let html = String::from_utf8(output.stdout).expect("HTML should be UTF-8");
    assert!(html.starts_with("<!DOCTYPE html>"));
    assert!(html.contains("src=\"https://cdn.jsdelivr.net/npm/vega@5.17.0\""));
    assert!(html.contains("fastqcx report"));
    assert!(
        html.len() < 500_000,
        "HTML unexpectedly contains fetched CDN payloads"
    );
}

#[test]
fn summary_preserves_fastqc_and_qctb_contracts() {
    let summary_directory = temporary_path("summary");
    let output = run_fastqcx(&[
        "--fastq",
        "tests/resources/example.fastq",
        "--summary",
        summary_directory
            .to_str()
            .expect("temporary path should be UTF-8"),
        "--no-html",
    ]);

    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let summary = fs::read_to_string(summary_directory.join("fastqc_data.txt"))
        .expect("summary should be published");
    assert!(summary.starts_with("##FastQC\t"));
    assert!(summary.contains(">>Basic Statistics\t"));
    assert!(summary.contains(">>Seqkit Statistics\tpass"));
    assert!(summary.contains("Q20(%)\tQ30(%)\tAvgQual"));

    assert!(summary.contains("1\t30.135\t33\t31\t34\t26\t34\n"));
    assert!(summary.contains(
        "1\t25.668449197860966\t22.459893048128343\t24.598930481283425\t27.27272727272727\n"
    ));

    let expected_positions: Vec<usize> = (1..=101).collect();
    for module_name in [
        "Per base sequence quality",
        "Per base sequence content",
        "Per base GC content",
        "Per base N content",
    ] {
        assert_eq!(
            extract_module_positions(&summary, module_name),
            expected_positions,
            "{} positions should be numerically ordered",
            module_name
        );
    }

    fs::remove_dir_all(summary_directory).expect("temporary summary should be removable");
}

#[test]
fn malformed_fastq_returns_nonzero_without_publishing_summary() {
    let input_path = temporary_path("malformed.fastq");
    let summary_directory = temporary_path("malformed-summary");
    fs::write(&input_path, "@valid\nAC\n+\nII\n@broken\nAC\n+\nI\n")
        .expect("malformed fixture should be written");

    let output = run_fastqcx(&[
        "--fastq",
        input_path.to_str().expect("temporary path should be UTF-8"),
        "--summary",
        summary_directory
            .to_str()
            .expect("temporary path should be UTF-8"),
        "--no-html",
    ]);

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("FASTQ record 2"),
        "unexpected stderr: {}",
        stderr
    );
    assert!(!summary_directory.exists());

    fs::remove_file(input_path).expect("temporary input should be removable");
}

#[test]
fn cli_rejects_unbounded_kmer_length() {
    let output = run_fastqcx(&[
        "--fastq",
        "tests/resources/example.fastq",
        "--kmer",
        "8",
        "--no-html",
    ]);

    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("expected 1..=7"));
}

#[test]
fn summary_reports_fastqc_style_duplication_and_overrepresented_sequences() {
    let fixture = temporary_path("duplication.fastq");
    let summary_directory = temporary_path("duplication-summary");
    let shared_prefix = "A".repeat(50);
    let records = [
        format!("{shared_prefix}{}", "C".repeat(30)),
        format!("{shared_prefix}{}", "G".repeat(30)),
        "C".repeat(80),
        "G".repeat(80),
    ];
    let mut fastq = String::new();
    for (index, sequence) in records.iter().enumerate() {
        fastq.push_str(&format!(
            "@read-{index}\n{sequence}\n+\n{}\n",
            "I".repeat(sequence.len())
        ));
    }
    fs::write(&fixture, fastq).expect("fixture should be written");

    let output = run_fastqcx(&[
        "--fastq",
        fixture.to_str().expect("fixture path should be UTF-8"),
        "--summary",
        summary_directory
            .to_str()
            .expect("summary path should be UTF-8"),
        "--no-html",
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );

    let summary = fs::read_to_string(summary_directory.join("fastqc_data.txt"))
        .expect("summary should be published");
    assert!(
        summary.contains(">>Sequence Duplication Levels\t"),
        "summary must expose FastQC's duplication-level module"
    );
    assert!(
        summary.contains("#Total Deduplicated Percentage\t75"),
        "three distinct sequence identities across four reads yield 75% deduplicated reads"
    );
    assert!(
        summary.contains("1\t50"),
        "two singleton reads account for half of all reads"
    );
    assert!(
        summary.contains("2\t50"),
        "the two long reads sharing their first 50 bases account for half of all reads"
    );
    assert!(
        summary.contains(">>Overrepresented sequences\t"),
        "summary must expose FastQC's overrepresented-sequences module"
    );
    assert!(
        summary.contains(&format!("{shared_prefix}\t2\t50")),
        "FastQC's 50-base identity rule must group long reads before reporting them"
    );

    fs::remove_file(fixture).expect("fixture should be removable");
    fs::remove_dir_all(summary_directory).expect("summary should be removable");
}
