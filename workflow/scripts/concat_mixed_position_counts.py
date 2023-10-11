import os
import sys


def concat_counts(counts: list[str], references: list[str]):
    results_dict: dict[str, int] = {}
    if not counts:
        print("There are no mixed position files to concat", file=sys.stderr)
        return results_dict

    print(f"There are {len(counts)} mixed_positions counts to copy", file=sys.stderr)
    for reference, count_file in zip(references, counts):
        with open(count_file, "r") as f:
            results_dict[reference] = int(f.readline().strip())
    return results_dict


def run_concat_counts(counts: list[str], references: list[str], out_file: str):
    if os.path.exists(out_file):
        os.remove(out_file)

    results = concat_counts(counts, references)

    with open(out_file, "w") as f:
        for ref, count in results.items():
            f.write(f"{ref}:{count}\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    run_concat_counts(
        snakemake.input.mixed_positions_counts,
        snakemake.params.reference_names,
        snakemake.output.summary,
    )
