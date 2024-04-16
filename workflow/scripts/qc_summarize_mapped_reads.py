import os
import sys

sys.stderr = open(snakemake.log[0], "w")


def parse_ref_name(filepath: str) -> str:
    return os.path.basename(filepath).removesuffix(".count")


def parse_count(filepath: str) -> int:
    with open(filepath, "r") as f:
        line = f.readline().strip()
        return int(line)


def summarize_mapped_reads(*, filepaths: list[str], output_file: str):
    ref_counts = {parse_ref_name(filepath): parse_count(filepath) for filepath in filepaths}

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)

    with open(output_file, "w") as f:
        for name, count in ref_counts.items():
            passed = "PASS" if count > 0 else "FAIL"
            f.write(f"{passed}\t{name}\t{count}\n")


if __name__ == "__main__":
    summarize_mapped_reads(filepaths=snakemake.input.read_counts, output_file=snakemake.output[0])
