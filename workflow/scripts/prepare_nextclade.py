import sys


def summarize_results(
    nextclade_out: str, others_out: str, metadata: dict[str, list[tuple[str, str, str]]], references: list[str]
):
    nextclades: list[tuple[str, str, str, str]] = []
    others: list[str] = []

    for reference in references:
        if reference in metadata:
            for mapping in metadata[reference]:
                nextclades.append((reference, *mapping))
        else:
            others.append(reference)

    with open(nextclade_out, "w") as f:
        for reference in nextclades:
            f.write(f"{reference[0]}\t{reference[1]}\t{reference[2]}\t{reference[3]}\n")

    with open(others_out, "w") as f:
        for reference in others:
            f.write(f"{reference}\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    summarize_results(
        nextclade_out=snakemake.output.nextclade,
        others_out=snakemake.output.others,
        metadata=snakemake.params.metadata,
        references=snakemake.params.references,
    )
