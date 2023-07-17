import json
import os
import sys
from typing import Literal, TypedDict


class NextcladeResult(TypedDict):
    segment: str
    relative_path_to_nextclade_dir: str


class ConsensusResult(TypedDict):
    category: Literal["nextclade", "other"]
    reference_name: str
    nextclade_results: list[NextcladeResult]


def summarize_results(others_csv: str, nextclade_tsv: list[str], nextclade_refs_file: str, out_summary_json: str):
    requested_tuples: list[tuple[str, str]] = []
    with open(nextclade_refs_file, "r") as f:
        for line in f.readlines():
            name, segment, _, _ = line.strip().split()
            requested_tuples.append((name, segment))

    other_references = []
    with open(others_csv, "r") as f:
        other_references = [line.strip() for line in f.readlines()]

    results: list[ConsensusResult] = []

    for ref, segment in requested_tuples:
        results.append({"category": "nextclade", "reference_name": ref, "nextclade_results": []})

    for tsv in nextclade_tsv:
        result_ref = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(tsv))))
        segment_nextclade = os.path.basename(os.path.dirname(tsv)).split("__")[0]

        for result in results:
            if result["reference_name"] == result_ref:
                result["nextclade_results"].append(
                    {
                        "segment": segment_nextclade,
                        "relative_path_to_nextclade_dir": os.path.relpath(
                            os.path.dirname(tsv), os.path.dirname(out_summary_json)
                        ),
                    }
                )
            break

    for ref in other_references:
        results.append({"category": "other", "reference_name": ref, "nextclade_results": []})

    with open(out_summary_json, "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    summarize_results(
        others_csv=snakemake.input.others,
        nextclade_tsv=snakemake.input.nextclade_tsv,
        nextclade_refs_file=snakemake.input.nextclade_refs,
        out_summary_json=snakemake.output[0],
    )
