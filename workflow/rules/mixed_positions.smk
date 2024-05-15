rule ivar__get_variants:
    input:
        bam="results/mapping/{sample}/deduplicated/{reference}.bam",
        bai="results/mapping/{sample}/deduplicated/{reference}.bam.bai",
        ref=infer_reference_fasta,
    output:
        tsv="results/variants/{sample}/{reference}/all.tsv",
    params:
        samtools_params=parse_samtools_params_for_variants(),
        ivar_params=parse_ivar_params_for_variants(),
    log:
        "logs/variants/{sample}/ivar/{reference}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/ivar/variants"


rule ivar__variants_to_vcf:
    input:
        "results/variants/{sample}/{reference}/all.tsv",
    output:
        all=report(
            "results/variants/{sample}/{reference}/all.vcf",
            category="{sample}",
            labels={
                "Type": "Variants - vcf",
                "Reference": "{reference}",
            },
        ),
        filtered="results/variants/{sample}/{reference}/passed_only.vcf",
    log:
        "logs/ivar/variants_to_vcf/{sample}/{reference}.log",
    localrule: True
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/ivar/vcf_convert"


rule report__ivar_variants_to_html:
    input:
        "results/variants/{sample}/{reference}/all.tsv",
    output:
        report(
            "results/variants/{sample}/{reference}/all.html",
            category="{sample}",
            labels={
                "Type": "Variants - all",
                "Reference": "{reference}",
            },
        ),
    log:
        "logs/report/{sample}/ivar_variants/{reference}.log",
    localrule: True
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/ivar/html_convert"


rule custom__compute_mixed_positions:
    input:
        "results/variants/{sample}/{reference}/all.tsv",
    output:
        mixed_positions="results/variants/{sample}/{reference}/mixed_positions.tsv",
        readcount="results/variants/{sample}/{reference}/mixed_positions_count.tsv",
    params:
        alt_depth=config["mixed_positions_params"]["filtering"]["min_alt_depth"],
        min_alt_freq=config["mixed_positions_params"]["filtering"]["min_alt_freq"],
        max_alt_freq=config["mixed_positions_params"]["filtering"]["max_alt_freq"],
        total_depth=config["mixed_positions_params"]["filtering"]["min_total_depth"],
    log:
        "logs/variants/{sample}/variants/{reference}.log",
    localrule: True
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/ivar/mixed_positions"


rule report__mixed_positions_to_html:
    input:
        "results/variants/{sample}/{reference}/mixed_positions.tsv",
    output:
        report(
            "results/variants/{sample}/{reference}/mixed_positions.html",
            caption="../report/mixed_positions.rst",
            category="{sample}",
            labels={
                "Type": "Variants - mixed positions",
                "Reference": "{reference}",
            },
        ),
    log:
        "logs/report/{sample}/mixed_positions/{reference}.log",
    localrule: True
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/ivar/html_convert"


rule custom__concat_mixed_position_counts:
    input:
        mixed_positions_counts=infer_mixed_positions_for_passed_references_only,
        reports=infer_variant_reports_for_passed_references_only,
    output:
        summary=report(
            "results/variants/{sample}/mixed_positions_summary.txt",
            category="{sample}",
            labels={
                "Type": "Mixed positions summary",
                "Reference": "-",
            },
        ),
    params:
        reference_names=lambda wildcards, input: [
            os.path.basename(os.path.dirname(filename)) for filename in input.mixed_positions_counts
        ],
    log:
        "logs/variants/{sample}/concat_mixed_position_counts.log",
    localrule: True
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/concat_mixed_position_counts.py"
