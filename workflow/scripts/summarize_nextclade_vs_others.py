import sys


def summarize_results(others_csv: str, nextclade_tsv: list[str], nextclade_refs: str, out_summary: str):
    nextclade_references = []
    with open(nextclade_refs, "r") as f:
        nextclade_references = [line.strip().split()[0] for line in f.readlines()]

    other_references = []
    with open(others_csv, "r") as f:
        other_references = [line.strip() for line in f.readlines()]

    with open(out_summary, "w") as f:
        f.write("category\treference_name\tnextclade_dir\n")
        for ref in nextclade_references:
            for tsv in nextclade_tsv:
                if ref in tsv:
                    f.write("nextclade\t{}\t{}\n".format(ref, tsv))
                    break
        for ref in other_references:
            f.write(f"other\t{ref}\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    summarize_results(
        others_csv=snakemake.input.others,
        nextclade_tsv=snakemake.input.nextclade_tsv,
        nextclade_refs=snakemake.input.nextclade_refs,
        out_summary=snakemake.output[0],
    )
