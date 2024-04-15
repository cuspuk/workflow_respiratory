import json
import sys
from dataclasses import dataclass


@dataclass
class SegmentCoverage:
    segment: str
    length: int
    coverage_at_depth: dict[int, float]


def process_json(json_file: str, depth: int):
    with open(json_file) as f:
        segments = [SegmentCoverage(**x) for x in json.load(f)]
    genome_length = 0
    total = 0
    for segment in segments:
        total += segment.coverage_at_depth[depth] * segment.length
        genome_length += segment.length
    return total / genome_length


def evaluate_mapping_quality(
    jsons: list[str],
    reference_names: list[str],
    min_genome_fraction_with_10x_coverage: float,
    output_for_passed: str,
    output_for_failed: str,
):
    passed_refs: dict[str, float] = {}
    failed_refs: dict[str, float] = {}
    for ref, json_file in zip(reference_names, jsons):
        avg = process_json(json_file, 10)
        if avg >= min_genome_fraction_with_10x_coverage:
            passed_refs[ref] = avg
        else:
            failed_refs[ref] = avg

    with open(output_for_passed, "w") as out_file:
        for passed_ref, coverage in passed_refs.items():
            out_file.write(f"{passed_ref}\t{coverage}\n")

    with open(output_for_failed, "w") as out_file:
        for failed_ref, coverage in failed_refs.items():
            out_file.write(f"{failed_ref}\t{coverage}\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_mapping_quality(
        snakemake.input.jsons,
        snakemake.params.reference_names,
        snakemake.params.min_genome_fraction_with_10x_coverage,
        snakemake.output.passed_refs,
        snakemake.output.failed_refs,
    )
