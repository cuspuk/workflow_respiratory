import sys

import pandas as pd
import panel as pn

sys.stderr = open(snakemake.log[0], "w")


def nextclade_to_html(tsv_input_path: str, output_html: str):
    print(f"Converting {tsv_input_path} into html", file=sys.stderr)
    nextclade_tsv = pd.read_csv(tsv_input_path, delimiter="\t")

    filters = {
        "seqName": {"type": "input", "func": "like", "placeholder": "Enter seqName"},
    }

    df = pn.widgets.Tabulator(nextclade_tsv, header_filters=filters, disabled=True, theme="modern")
    df.save(output_html)


if __name__ == "__main__":
    nextclade_to_html(tsv_input_path=snakemake.input[0], output_html=snakemake.output[0])
