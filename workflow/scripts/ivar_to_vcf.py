# Inspired by https://github.com/jts/ncov-tools/blob/master/workflow/scripts/ivar_variants_to_vcf.py

import os
import re
import sys

sys.stderr = open(snakemake.log[0], "w")


def _prepare_header(filename: str):
    header = (
        "##fileformat=VCFv4.2\n"
        "##source=iVar\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
        '##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">\n'
        '##FILTER=<ID=FAIL,Description="Result of p-value > 0.05">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">\n'
        '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">\n'
        '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">\n'
        '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">\n'
        '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Depth of alternate base on reverse reads">\n'
        '##FORMAT=<ID=ALT_QUAL,Number=1,Type=String,Description="Mean quality of alternate base">\n'
        '##FORMAT=<ID=ALT_FREQ,Number=1,Type=String,Description="Frequency of alternate base">\n'
    )
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + filename + "\n"
    return header


def _prepare_output(output_vcf_path: str):
    out_dir = os.path.dirname(output_vcf_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


def ivar_variants_to_vcf(
    *, ivar_tsv_path: str, all_vcf_path: str, filtered_vcf_path: str | None, min_allele_freq: int = 0
):
    filename = os.path.splitext(ivar_tsv_path)[0]
    header = _prepare_header(filename)

    var_list: list[tuple[str, str, str, str]] = []

    all_lines: list[str] = []
    filt_lines: list[str] = []

    fout = open(all_vcf_path, "w")
    fout.write(header)

    with open(ivar_tsv_path, "r") as f:
        for line in f:
            if re.match("REGION", line):
                continue

            tab_values: list[str] = re.split("\t", line)
            chrom, pos, ref, alt = tab_values[0], tab_values[1], tab_values[2], tab_values[3]
            passed = "PASS" if tab_values[13] == "TRUE" else "FAIL"

            if alt[0] == "+":
                alt = ref + alt[1:]
            elif alt[0] == "-":
                ref += alt[1:]
                alt = tab_values[2]

            info = "DP=" + tab_values[11]
            FORMAT = "GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ"
            SAMPLE = ":".join(["1"] + tab_values[4:11])
            out_line = "\t".join([chrom, pos, ".", ref, alt, ".", passed, info, FORMAT, SAMPLE])

            if (chrom, pos, ref, alt) not in var_list:
                var_list.append((chrom, pos, ref, alt))
                all_lines.append(out_line)
                if filtered_vcf_path and passed == "PASS" and float(tab_values[10]) >= min_allele_freq:
                    filt_lines.append(out_line)

    print(f"Found {len(all_lines)} lines", file=sys.stderr)
    print(f"After filtering there are {len(filt_lines)} lines", file=sys.stderr)

    _prepare_output(all_vcf_path)
    with open(all_vcf_path, "w") as f:
        f.write(header)
        for line in all_lines:
            f.write(line)

    if filtered_vcf_path:
        _prepare_output(all_vcf_path)
        with open(filtered_vcf_path, "w") as f:
            f.write(header)
            for line in filt_lines:
                f.write(line)


if __name__ == "__main__":
    all_variants = snakemake.output.all
    filtered = snakemake.output.get("filtered", None)
    ivar_variants_to_vcf(
        ivar_tsv_path=snakemake.input[0],
        all_vcf_path=all_variants,
        filtered_vcf_path=filtered,
        min_allele_freq=0,
    )
