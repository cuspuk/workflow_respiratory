import json
import os
import sys
from dataclasses import asdict, dataclass

MIN_DEPTH = 10


@dataclass
class SegmentCoverage:
    segment: str
    length: int
    coverage_at_depth: dict[int, float]


@dataclass
class Reference:
    name: str
    segments: list[SegmentCoverage]

    def average_depth(self, depth: int):
        genome_length = 0
        total = 0
        for segment in self.segments:
            total += segment.coverage_at_depth[depth] * segment.length
            genome_length += segment.length
        return total / genome_length


def coverage_with_at_least(depths: list[int], threshold: int) -> float:
    return len([d for d in depths if d >= threshold]) / len(depths)


def load_depths(filename: str) -> list[SegmentCoverage]:
    rows: dict[str, list[int]] = {}
    with open(filename) as f:
        for line in f.readlines():
            segment, _, depth = line.strip().split("\t")
            if segment not in rows:
                rows[segment] = []
            rows[segment].append(int(depth))

    segments: list[SegmentCoverage] = []
    for segment, depths in rows.items():
        coverages = {j: coverage_with_at_least(depths, j) for j in range(1, 100)}
        segments.append(SegmentCoverage(segment, len(depths), coverages))

    return segments


def evaluate_mapping_quality(
    depths_files: list[str],
    reference_names: list[str],
    threshold: float,
    tsv_out: str,
    json_out: str,
):
    refs: list[Reference] = []
    for ref, depths_file in zip(reference_names, depths_files):
        refs.append(Reference(name=ref, segments=load_depths(depths_file)))

    os.makedirs(os.path.dirname(os.path.abspath(tsv_out)), exist_ok=True)
    with open(tsv_out, "w") as out_file:
        for ref in refs:
            avg = ref.average_depth(MIN_DEPTH)
            passed = "PASS" if avg >= threshold else "FAIL"
            out_file.write(f"{passed}\t{ref.name}\t{avg}\n")

    os.makedirs(os.path.dirname(os.path.abspath(json_out)), exist_ok=True)
    with open(json_out, "w") as out_file:
        refs_dict = [asdict(x) for x in refs]
        json.dump(refs_dict, out_file, indent=4)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_mapping_quality(
        snakemake.input.depths,
        snakemake.params.reference_names,
        snakemake.params.threshold,
        snakemake.output.tsv,
        snakemake.output.json,
    )
