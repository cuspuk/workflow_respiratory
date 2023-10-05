rule ivar__get_variants:
    input:
        bam="results/mapping/{sample}/deduplicated/{reference}.bam",
        bai="results/mapping/{sample}/deduplicated/{reference}.bam.bai",
        ref=get_reference_fasta,
    output:
        "results/variants/{sample}/{reference}/all.tsv",
    params:
        samtools_params=parse_samtools_params_for_variants(),
        ivar_params=parse_ivar_params_for_variants(),
    log:
        "logs/variants/{sample}/ivar/{reference}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.6.0/wrappers/ivar/variants"


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
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/ivar_to_vcf.py"


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
    conda:
        "../envs/python_panel.yaml"
    script:
        "../scripts/tsv_html.py"


rule custom__compute_mixed_positions:
    input:
        "results/variants/{sample}/{reference}/all.tsv",
    output:
        mixed_positions="results/variants/{sample}/{reference}/mixed_positions.tsv",
        readcount=temp("results/variants/{sample}/{reference}/mixed_positions_count.tsv"),
    params:
        alt_depth=config["mixed_positions_params"]["filtering"]["min_alt_depth"],
        min_alt_freq=config["mixed_positions_params"]["filtering"]["min_alt_freq"],
        max_alt_freq=config["mixed_positions_params"]["filtering"]["max_alt_freq"],
        total_depth=config["mixed_positions_params"]["filtering"]["min_total_depth"],
    log:
        "logs/variants/{sample}/variants/{reference}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/compute_mixed_positions.py"


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
    conda:
        "../envs/python_panel.yaml"
    script:
        "../scripts/tsv_html.py"


rule custom__concatenate_mixed_positions:
    input:
        mixed_positions=get_mixed_positions_for_passed_references_only,
        reports=get_variant_reports_for_passed_references_only,
    output:
        report(
            "results/variants/{sample}/mixed_positions_summary.txt",
            category="{sample}",
            labels={
                "Type": "Mixed positions summary",
                "Reference": "-",
            },
        ),
    log:
        "logs/variants/{sample}/summary.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input.mixed_positions} > {output} 2>&1"
