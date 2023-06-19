import pandas as pd
import panel as pn


def ivar_variants_to_html(tsv_input_path: str, output_html: str, reference_name: str):
    ivar_df = pd.read_csv(tsv_input_path, delimiter="\t")

    tabulator_formatters = {"PASS": {"type": "tickCross"}}

    filters = {
        "REGION": {"type": "input", "func": "like", "placeholder": "Enter region"},
        "ALT_FREQ": {"type": "number", "func": ">=", "placeholder": "Enter minimum"},
        "ALT_DP": {"type": "number", "func": ">=", "placeholder": "Enter minimum"},
        "TOTAL_DP": {"type": "number", "func": ">=", "placeholder": "Enter minimum"},
    }

    df = pn.widgets.Tabulator(
        ivar_df, header_filters=filters, disabled=True, theme="modern", formatters=tabulator_formatters
    )
    df.save(output_html, title=reference_name)


if __name__ == "__main__":
    ivar_variants_to_html(
        tsv_input_path=snakemake.input[0], output_html=snakemake.output[0], reference_name=snakemake.wildcards.reference
    )
