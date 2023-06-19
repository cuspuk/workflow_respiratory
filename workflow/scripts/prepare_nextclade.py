import os
import sys


def summarize_results(references_file: str, nextclade_out: str, others_out: str):
    mapping = {
        "sars_cov_2": ("sars-cov-2", "default"),  # ???
        "rsv_a_2017": ("rsv_a", "default"),
        "rsv_b_2019": ("rsv_b", "default"),
        "yamagata_2013": ("flu_yam_ha", "default"),
        "victoria_2021": ("flu_vic_ha", "default"),  # flu_vic_na ??
        "h1n1_2019": ("flu_h1n1pdm_ha", "default"),  # ci NA?
        "h3n2_2021": ("flu_h3n2_ha", "default"),  # a NA ?
        "monkeypox": ("hMPXV", "default"),  # ci ine?
    }

    nextclades: dict[str, tuple[str, str]] = {}
    others: list[str] = []

    references: list[str] = [line.strip() for line in open(references_file, "r").readlines()]
    print(f"Found {len(references)} references", file=sys.stderr)
    for reference in references:
        if reference in mapping:
            nextclades[reference] = mapping[reference]
        else:
            others.append(reference)

    with open(nextclade_out, "w") as f:
        for reference in nextclades:
            f.write(f"{reference}\t{nextclades[reference][0]}\t{nextclades[reference][1]}\n")

    with open(others_out, "w") as f:
        for reference in others:
            f.write(f"{reference}\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    summarize_results(snakemake.input[0], snakemake.output.nextclade, snakemake.output.others)
