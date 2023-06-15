rule ivar__get_variants:
    input:
        bam="results/mapping/{reference}/deduplicated/{sample}.bam",
        bai="results/mapping/{reference}/deduplicated/{sample}.bam.bai",
        ref=get_reference_fasta,
    output:
        "results/variants/{sample}/{reference}/original.tsv",
    params:
        out_prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        samtools_params=parse_samtools_params_for_variants(),
        ivar_params=parse_ivar_params_for_variants(),
    log:
        "logs/variants/{sample}/ivar/{reference}.log",
    conda:
        "../envs/ivar.yaml"
    shell:
        "(samtools mpileup --no-BAQ {params.samtools_params} {input.bam}"
        " |"
        " ivar variants -p {params.out_prefix} -r {input.ref} {params.ivar_params}) 1> {log} 2>&1"


rule custom__compute_mixed_positions:
    input:
        "results/variants/{sample}/{reference}.tsv",
    output:
        mixed="results/variants/{sample}/mixed_positions/{reference}.tsv",
        count=temp("results/variants/{sample}/mixed_positions/{reference}_count.tsv"),
    params:
        alt_depth=2,
        alt_freq=0.1,
        total_depth=10,
    log:
        "logs/variants/{sample}/variants/{reference}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/compute_mixed_positions.py"


rule custom__concatenate_mixed_positions:
    input:
        get_mixed_positions_for_passed_references_only(),
    output:
        "results/variants/{sample}/mixed_positions_summary.tsv",
    log:
        "logs/variants/{sample}/summary.log",
    shell:
        "cat {input} > {output} 2>&1"
