import os
import sys

sys.stderr = open(snakemake.log[0], "w")


def _parse_reference_name(filepath: str) -> str:
    # results/mapping/{sample}/deduplicated/{reference}.count
    return os.path.basename(filepath).removesuffix(".count")


def _count_is_zero(filepath: str) -> int:
    with open(filepath, "r") as f:
        line = f.readline().strip()
        return int(line) == 0


def parse_read_counts_in_bams(*, filepaths: list[str], output_file: str):
    not_empty = [_parse_reference_name(filepath) for filepath in filepaths if not _count_is_zero(filepath)]

    with open(output_file, "w") as f:
        f.write("\n".join(not_empty))


if __name__ == "__main__":
    parse_read_counts_in_bams(filepaths=snakemake.input.read_counts, output_file=snakemake.output[0])
