__author__ = "Palo Misenko"

import os
import sys

import numpy as np
import pandas as pd


def assign_pass(coverage: float) -> str:
    if coverage > 1 or coverage < 0:
        raise ValueError(f'Coverage "{coverage}" not in 0 - 1 interval!')
    if coverage >= 0.9:
        return "PASS"
    elif coverage < 0.9 and coverage >= 0.5:
        return "WARN"
    return "FAIL"


def group_coverage_filter(group) -> list[int]:
    if len(group) == 1:
        return [group.index[0]]
    group_max = group["coverage"].max()

    # Calculate differences between each value and maximum in group
    cov_diffs = group["coverage"].apply(lambda x: group_max - x)
    cov = group["coverage"].index

    # Leave just indexes where coverage differs less than 0.3 from max in group
    results = map(lambda x: x[0], filter(lambda x: x[1] < 0.5, zip(cov, cov_diffs)))
    return list(results)


def merge_nextclade_tsvs(files: list[str], output_tsv: str):
    parent_dir = os.path.dirname(output_tsv)
    if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)

    if not files:
        print("No files to process. Leaving empty output file", file=sys.stderr)
        with open(output_tsv, "w") as _:
            pass
        return

    print(f"Processing {len(files)} files", file=sys.stderr)
    all_attributes = set()
    dataframes = []
    for file in files:
        if os.path.getsize(file) == 0:
            print(f"File {file} is empty. Continuing...", file=sys.stderr)
            continue

        df = pd.read_csv(file, sep="\t")
        if "type" not in df:
            df["type"] = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(file))))

        df["QC"] = df["coverage"].apply(lambda x: assign_pass(x))

        if "index" in df.columns:
            df.drop(columns=["index"], inplace=True)
        # Only accept samples with specified clade and coverage
        # Samples without clade or coverage are wrongly mapped sequences
        df = df[pd.notnull(df["clade"]) & pd.notnull(df["coverage"])]
        dataframes.append(df)
        all_attributes = all_attributes.union(set(df.columns))

    if len(dataframes) == 0:
        print("All files have been empty. Leaving empty output file", file=sys.stderr)
        with open(output_tsv, "w") as _:
            pass
        return

    # Normalize tables to same columns
    for df in dataframes:
        for attribute in all_attributes:
            if attribute not in df.columns:
                # Insert empty column to dataframe
                df.insert(loc=4, column=attribute, value=["" for _ in range(df.shape[0])])

    # Use ordering of biggest table
    index = np.asarray(map(len, dataframes)).argmax()
    columns_order = list(dataframes[index].columns)

    # Reorder columns of tables and concat
    df_final = pd.concat(map(lambda x: x[columns_order], dataframes))

    # Visual updates of table and final reorder
    for column in ["RBD", "short_clade", "glycosylation"]:
        if column in df_final.columns:
            df_final.insert(6, column, df_final.pop(column))
    i = 0
    for column in ["seqName", "type", "clade", "G_clade", "Nextclade_pango", "QC", "coverage"]:
        if column in df_final.columns:
            df_final.insert(i, column, df_final.pop(column))
            i += 1

    df_final.sort_values(["seqName", "coverage"], inplace=True)
    df_final.reset_index(inplace=True, drop=True)

    df_final.to_csv(output_tsv, index=False, sep="\t")

    # Delete samples where coverage is too small.
    # Group by 'seqName' and determine which rows to keep for each group
    # idx = df_final.groupby(['seqName']).apply(group_coverage_filter).explode()
    # df_final.iloc[idx].to_csv(output_tsv, index=False, sep='\t')


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    merge_nextclade_tsvs(
        snakemake.input.nextclade_tsvs,
        snakemake.output.merged_tsv,
    )
