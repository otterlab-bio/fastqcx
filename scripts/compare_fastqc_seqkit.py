#!/usr/bin/env python3
"""Compare fastqcx output with pinned FastQC and SeqKit reference outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


FASTQC_NUMERIC_TOLERANCE = 0.01


def parse_modules(text: str) -> dict[str, dict[str, object]]:
    modules: dict[str, dict[str, object]] = {}
    current: dict[str, object] | None = None
    for line in text.splitlines():
        ends_module = line.endswith(">>END_MODULE")
        if ends_module:
            line = line.removesuffix(">>END_MODULE")
        if line.startswith(">>"):
            name, status = line[2:].split("\t", 1)
            current = {"status": status, "comments": [], "rows": []}
            modules[name] = current
        elif current is not None and line:
            key = "comments" if line.startswith("#") else "rows"
            current[key].append(line.split("\t"))
        if ends_module:
            current = None
    return modules


def parse_statistic(module: dict[str, object], name: str) -> str:
    for section in ("comments", "rows"):
        for row in module[section]:
            if row[0].lstrip("#") == name:
                return row[1]
    raise AssertionError(f"missing {name!r} in module")


def parse_base_range(label: str) -> range:
    bounds = [int(value) for value in label.split("-", 1)]
    return range(bounds[0], bounds[-1] + 1)


def compare_per_base_quality(
    fastqcx_modules: dict[str, dict[str, object]],
    fastqc_modules: dict[str, dict[str, object]],
) -> dict[str, object]:
    name = "Per base sequence quality"
    actual = fastqcx_modules[name]
    expected = fastqc_modules[name]
    assert actual["status"] == expected["status"], (
        name,
        actual["status"],
        expected["status"],
    )

    actual_rows = {int(row[0]): [float(value) for value in row[1:]] for row in actual["rows"]}
    columns = [
        "mean",
        "median",
        "lower_quartile",
        "upper_quartile",
        "percentile_10",
        "percentile_90",
    ]
    maximum_difference = 0.0
    values_compared = 0
    for expected_row in expected["rows"]:
        positions = list(parse_base_range(expected_row[0]))
        grouped_actual = [
            sum(actual_rows[position][column] for position in positions) / len(positions)
            for column in range(len(columns))
        ]
        for column, actual_value, expected_value in zip(
            columns, grouped_actual, map(float, expected_row[1:])
        ):
            difference = abs(actual_value - expected_value)
            assert difference <= FASTQC_NUMERIC_TOLERANCE, (
                name,
                expected_row[0],
                column,
                actual_value,
                expected_value,
                difference,
            )
            maximum_difference = max(maximum_difference, difference)
            values_compared += 1

    return {
        "status": actual["status"],
        "position_groups_compared": len(expected["rows"]),
        "values_compared": values_compared,
        "tolerance": FASTQC_NUMERIC_TOLERANCE,
        "maximum_difference": round(maximum_difference, 6),
    }


def compare_duplication(
    fastqcx_modules: dict[str, dict[str, object]],
    fastqc_modules: dict[str, dict[str, object]],
) -> dict[str, object]:
    name = "Sequence Duplication Levels"
    actual = fastqcx_modules[name]
    expected = fastqc_modules[name]
    assert actual["status"] == expected["status"], (
        name,
        actual["status"],
        expected["status"],
    )

    actual_deduplicated = float(parse_statistic(actual, "Total Deduplicated Percentage"))
    expected_deduplicated = float(parse_statistic(expected, "Total Deduplicated Percentage"))
    deduplicated_difference = abs(actual_deduplicated - expected_deduplicated)
    assert deduplicated_difference <= FASTQC_NUMERIC_TOLERANCE, (
        name,
        "Total Deduplicated Percentage",
        actual_deduplicated,
        expected_deduplicated,
        deduplicated_difference,
    )

    actual_levels = {row[0]: float(row[1]) for row in actual["rows"]}
    expected_levels = {row[0]: float(row[1]) for row in expected["rows"]}
    assert actual_levels.keys() == expected_levels.keys(), (
        actual_levels.keys(),
        expected_levels.keys(),
    )
    maximum_difference = deduplicated_difference
    for level, expected_percentage in expected_levels.items():
        difference = abs(actual_levels[level] - expected_percentage)
        assert difference <= FASTQC_NUMERIC_TOLERANCE, (
            name,
            level,
            actual_levels[level],
            expected_percentage,
            difference,
        )
        maximum_difference = max(maximum_difference, difference)

    return {
        "status": actual["status"],
        "deduplicated_percentage": actual_deduplicated,
        "levels_compared": len(actual_levels),
        "tolerance": FASTQC_NUMERIC_TOLERANCE,
        "maximum_difference": round(maximum_difference, 6),
    }


def compare_overrepresented_sequences(
    fastqcx_modules: dict[str, dict[str, object]],
    fastqc_modules: dict[str, dict[str, object]],
) -> dict[str, object]:
    name = "Overrepresented sequences"
    actual = fastqcx_modules[name]
    expected = fastqc_modules[name]
    assert actual["status"] == expected["status"], (
        name,
        actual["status"],
        expected["status"],
    )

    actual_sequences = {
        row[0]: {"count": int(row[1]), "percentage": float(row[2])}
        for row in actual["rows"]
    }
    expected_sequences = {
        row[0]: {"count": int(row[1]), "percentage": float(row[2])}
        for row in expected["rows"]
    }
    assert actual_sequences.keys() == expected_sequences.keys(), (
        actual_sequences.keys(),
        expected_sequences.keys(),
    )
    maximum_difference = 0.0
    for sequence, expected_values in expected_sequences.items():
        actual_values = actual_sequences[sequence]
        assert actual_values["count"] == expected_values["count"], (
            name,
            sequence,
            actual_values["count"],
            expected_values["count"],
        )
        difference = abs(actual_values["percentage"] - expected_values["percentage"])
        assert difference <= FASTQC_NUMERIC_TOLERANCE, (
            name,
            sequence,
            actual_values["percentage"],
            expected_values["percentage"],
            difference,
        )
        maximum_difference = max(maximum_difference, difference)

    return {
        "status": actual["status"],
        "sequences_compared": len(actual_sequences),
        "tolerance": FASTQC_NUMERIC_TOLERANCE,
        "maximum_difference": round(maximum_difference, 6),
    }


def compare_seqkit(fastqcx_modules: dict[str, dict[str, object]], seqkit_text: str) -> dict[str, object]:
    module = fastqcx_modules["Seqkit Statistics"]
    header = module["comments"][0]
    values = module["rows"][0]
    fastqcx_stats = dict(zip([field.lstrip("#") for field in header], values))

    seqkit_lines = seqkit_text.splitlines()
    seqkit_stats = dict(zip(seqkit_lines[0].split("\t"), seqkit_lines[1].split("\t")))
    assert fastqcx_stats["format"] == seqkit_stats["format"] == "FASTQ"
    assert fastqcx_stats["type"] == seqkit_stats["type"] == "DNA"

    integer_fields = [
        "num_seqs",
        "sum_len",
        "min_len",
        "max_len",
        "Q1",
        "Q2",
        "Q3",
        "sum_gap",
        "N50",
        "N50_num",
        "sum_n",
    ]
    integer_statistics = {}
    for field in integer_fields:
        actual = int(float(fastqcx_stats[field]))
        expected = int(float(seqkit_stats[field]))
        assert actual == expected, (field, actual, expected)
        integer_statistics[field] = actual

    decimal_tolerances = {
        "avg_len": 0.05,
        "Q20(%)": 0.5,
        "Q30(%)": 0.5,
        "AvgQual": 0.01,
        "GC(%)": 0.01,
    }
    decimal_statistics = {}
    largest_difference = 0.0
    for field, tolerance in decimal_tolerances.items():
        actual = float(fastqcx_stats[field])
        expected = float(seqkit_stats[field])
        difference = abs(actual - expected)
        assert difference <= tolerance, (field, actual, expected, difference)
        largest_difference = max(largest_difference, difference)
        decimal_statistics[field] = {
            "fastqcx": actual,
            "seqkit": expected,
            "tolerance": tolerance,
            "difference": round(difference, 6),
        }

    return {
        "integer_fields_compared": len(integer_fields),
        "integer_fields_matched": len(integer_fields),
        "integer_statistics": integer_statistics,
        "decimal_fields_compared": len(decimal_tolerances),
        "decimal_fields_within_tolerance": len(decimal_tolerances),
        "decimal_statistics": decimal_statistics,
        "largest_decimal_deviation": round(largest_difference, 6),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastqcx-summary", type=Path, required=True)
    parser.add_argument("--fastqc-summary", type=Path, required=True)
    parser.add_argument("--seqkit-stats", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    fastqcx_text = arguments.fastqcx_summary.read_text()
    fastqc_text = arguments.fastqc_summary.read_text()
    fastqcx_modules = parse_modules(fastqcx_text)
    fastqc_modules = parse_modules(fastqc_text)

    seqkit = compare_seqkit(fastqcx_modules, arguments.seqkit_stats.read_text())
    fastqc_total = int(parse_statistic(fastqc_modules["Basic Statistics"], "Total Sequences"))
    fastqc_gc = round(float(parse_statistic(fastqc_modules["Basic Statistics"], "%GC")))
    fastqcx_total = int(parse_statistic(fastqcx_modules["Basic Statistics"], "Total Sequences"))
    fastqcx_gc = round(float(parse_statistic(fastqcx_modules["Basic Statistics"], "%GC")))
    assert fastqcx_total == fastqc_total == seqkit["integer_statistics"]["num_seqs"]
    assert fastqcx_gc == fastqc_gc

    fastqc_module_results = {
        "per_base_sequence_quality": compare_per_base_quality(fastqcx_modules, fastqc_modules),
        "sequence_duplication_levels": compare_duplication(fastqcx_modules, fastqc_modules),
        "overrepresented_sequences": compare_overrepresented_sequences(
            fastqcx_modules, fastqc_modules
        ),
    }
    report = {
        "schema_version": "otter.fastqcx-parity/v2",
        **seqkit,
        "official_fastqc_total_sequences": fastqc_total,
        "official_fastqc_gc_percent": fastqc_gc,
        "fastqc_modules": fastqc_module_results,
        "fastqc_per_base_values_compared": fastqc_module_results[
            "per_base_sequence_quality"
        ]["values_compared"],
        "fastqc_per_base_maximum_difference": fastqc_module_results[
            "per_base_sequence_quality"
        ]["maximum_difference"],
        "fastqc_duplication_levels_compared": fastqc_module_results[
            "sequence_duplication_levels"
        ]["levels_compared"],
        "fastqc_duplication_maximum_difference": fastqc_module_results[
            "sequence_duplication_levels"
        ]["maximum_difference"],
        "fastqc_overrepresented_sequences_compared": fastqc_module_results[
            "overrepresented_sequences"
        ]["sequences_compared"],
        "fastqc_overrepresented_maximum_difference": fastqc_module_results[
            "overrepresented_sequences"
        ]["maximum_difference"],
        "equal": True,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(
        "parity report: "
        f"SeqKit {report['integer_fields_matched']}/{report['integer_fields_compared']} exact integers, "
        f"FastQC per-base {report['fastqc_modules']['per_base_sequence_quality']['values_compared']} values, "
        f"duplication {report['fastqc_modules']['sequence_duplication_levels']['levels_compared']} levels, "
        f"overrepresented {report['fastqc_modules']['overrepresented_sequences']['sequences_compared']} sequences"
    )


if __name__ == "__main__":
    main()
