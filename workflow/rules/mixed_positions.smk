rule ivar__get_variants:
    input:
        bam="results/mapping/{reference}/deduplicated/{sample}.bam",
        bai="results/mapping/{reference}/deduplicated/{sample}.bam.bai",
        ref=get_reference_fasta,
    output:
        "results/variants/{sample}/{reference}/all.tsv",
    params:
        out_prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        samtools_params=parse_samtools_params_for_variants(),
        ivar_params=parse_ivar_params_for_variants(),
    log:
        "logs/variants/{sample}/ivar/{reference}.log",
    conda:
        "../envs/ivar.yaml"
    shell:
        "(samtools mpileup -aa --no-BAQ {params.samtools_params} {input.bam}"
        " |"
        " ivar variants -p {params.out_prefix} -r {input.ref} {params.ivar_params}) 1> {log} 2>&1"


rule report__ivar_variants_to_html:
    input:
        "results/variants/{sample}/{reference}/all.tsv",
    output:
        report(
            "results/variants/{sample}/{reference}/all.html",
            category="Mixed positions",
            labels={
                "reference": "{reference}",
                "step": "All variants",
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
            category="Mixed positions",
            caption="../report/mixed_positions.rst",
            labels={
                "reference": "{reference}",
                "step": "Only mixed positions",
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
            labels={
                "Name": "Summary of mixed positions",
            },
        ),
    log:
        "logs/variants/{sample}/summary.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input.mixed_positions} > {output} 2>&1"
