import json
import os
import sys
from dataclasses import asdict, dataclass, field


@dataclass
class SegmentCoverage:
    segment: str
    length: int
    depths: list[int] = field(default_factory=list)
    coverage_at_depth: dict[int, float] = field(default_factory=dict)

    def __post_init__(self):
        self.coverage_min_depth = {j: self.coverage_with_at_least(j) for j in range(1, 100)}

    def coverage_with_at_least(self, threshold: int) -> float:
        return len([d for d in self.depths if d >= threshold]) / self.length


def load_depths(filename: str) -> list[SegmentCoverage]:
    rows: dict[str, list[int]] = {}
    with open(filename) as f:
        for line in f.readlines():
            segment, _, depth = line.strip().split("\t")
            if segment not in rows:
                rows[segment] = []
            rows[segment].append(int(depth))

    return [SegmentCoverage(segment, len(depths), depths) for segment, depths in rows.items()]


def produce_depths(filename: str, output: str):
    depths = load_depths(filename)
    depths_dct = [asdict(x, dict_factory=lambda x: {k: v for (k, v) in x if k != "depths"}) for x in depths]
    json_object = json.dumps(depths_dct, indent=4)

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, "w") as f:
        f.write(json_object)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    produce_depths(
        filename=snakemake.input.txt,
        output=snakemake.output.json,
    )
