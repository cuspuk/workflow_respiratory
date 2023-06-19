import sys


class InvalidMetadataFile(Exception):
    """Raised when the metadata file cannot be parsed into 4 values"""


def summarize_results(references_file: str, nextclade_out: str, others_out: str, metadata: dict[str, tuple[str, str]]):
    nextclades: dict[str, tuple[str, str]] = {}
    others: list[str] = []

    references: list[str] = [line.strip() for line in open(references_file, "r").readlines()]
    print(f"Found {len(references)} references", file=sys.stderr)

    for reference in references:
        if metadata[reference][1]:
            nextclades[reference] = metadata[reference]
        else:
            others.append(reference)

    with open(nextclade_out, "w") as f:
        for reference in nextclades:
            f.write(f"{reference}\t{nextclades[reference][0]}\t{nextclades[reference][1]}\n")

    with open(others_out, "w") as f:
        for reference in others:
            f.write(f"{reference}\n")


def load_metadata(metadata_file: str):
    mapping: dict[str, tuple[str, str]] = {}
    with open(metadata_file, "r") as f:
        for line in f.readlines():
            try:
                name, _, nextclade, tag = line.split(",")
                mapping[name] = (nextclade, tag)
            except ValueError:
                raise InvalidMetadataFile("Metadata table {} does not have 4 columns".format(metadata_file))
    return mapping


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    metadata = load_metadata(snakemake.input.metadata)
    summarize_results(
        references_file=snakemake.input.references,
        nextclade_out=snakemake.output.nextclade,
        others_out=snakemake.output.others,
        metadata=metadata,
    )
