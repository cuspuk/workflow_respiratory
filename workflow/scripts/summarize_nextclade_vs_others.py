import json
import os
import sys
from typing import Literal, TypedDict


class ConsensusResult(TypedDict):
    category: Literal["nextclade", "other"]
    reference_name: str
    virus: str
    relative_path_to_nextclade_dir: str | None


class InvalidMetadataFile(Exception):
    """Raised when the metadata file cannot be parsed into 4 values"""


def load_metadata(metadata_file: str):
    mapping: dict[str, str] = {}
    with open(metadata_file, "r") as f:
        for line in f.readlines():
            try:
                name, virus, _, _ = line.strip().split(",")
                mapping[name] = virus
            except ValueError:
                raise InvalidMetadataFile("Metadata table {} does not have 4 columns".format(metadata_file))
    return mapping


def summarize_results(
    others_csv: str, nextclade_tsv: list[str], nextclade_refs: str, out_summary_json: str, mapping: dict[str, str]
):
    nextclade_references = []
    with open(nextclade_refs, "r") as f:
        nextclade_references = [line.strip().split()[0] for line in f.readlines()]

    other_references = []
    with open(others_csv, "r") as f:
        other_references = [line.strip() for line in f.readlines()]

    results: list[ConsensusResult] = []
    for ref in nextclade_references:
        for tsv in nextclade_tsv:
            if ref in tsv:
                results.append(
                    {
                        "category": "nextclade",
                        "reference_name": ref,
                        "virus": mapping[ref],
                        "relative_path_to_nextclade_dir": os.path.relpath(os.path.dirname(tsv), out_summary_json),
                    }
                )
                break
    for ref in other_references:
        results.append(
            {"category": "other", "reference_name": ref, "virus": mapping[ref], "relative_path_to_nextclade_dir": None}
        )

    with open(out_summary_json, "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    metadata = load_metadata(snakemake.input.metadata)
    summarize_results(
        others_csv=snakemake.input.others,
        nextclade_tsv=snakemake.input.nextclade_tsv,
        nextclade_refs=snakemake.input.nextclade_refs,
        out_summary_json=snakemake.output[0],
        mapping=metadata,
    )
