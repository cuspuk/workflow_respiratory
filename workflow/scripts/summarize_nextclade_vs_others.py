import sys


def summarize_results(others_csv: str, nextclade_files: list[str], passed_references: list[str], out_summary: str):
    lines = []
    with open(input.others, "r") as f:
        lines = f.readlines()
    with open(out_summary, "w") as f:
        for line in lines:
            f.write(line + "\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    summarize_results(
        others_csv=snakemake.input.others,
        nextclade_files=snakemake.input.nextclade,
        passed_references=snakemake.input.passed_references,
        out_summary=snakemake.output,
    )
