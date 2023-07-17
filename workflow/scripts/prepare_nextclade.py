import sys


class InvalidMetadataFile(Exception):
    """Raised when the metadata file cannot be parsed into 4 values"""


def summarize_results(
    references_file: str, nextclade_out: str, others_out: str, metadata: dict[str, list[tuple[str, str, str]]]
):
    nextclades: list[tuple[str, str, str, str]] = []
    others: list[str] = []

    references: list[str] = [line.strip() for line in open(references_file, "r").readlines()]
    print(f"Found {len(references)} references", file=sys.stderr)

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


def load_metadata(metadata_file: str):
    mapping: dict[str, list[tuple[str, str, str]]] = {}
    with open(metadata_file, "r") as f:
        for line in f.readlines():
            try:
                name, segment, nextclade, tag = line.strip().split(",")
                tag = tag if tag else "default"
                if name not in mapping:
                    mapping[name] = []
                mapping[name].append((segment, nextclade, tag))
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
