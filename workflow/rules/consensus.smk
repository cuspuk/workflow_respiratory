checkpoint index_passed_references:
    input:
        reference=get_reference_fasta,
    output:
        os.path.join(config["reference_panel_dirpath"], "{reference}.fa.fai"),
    log:
        "logs/checkpoints/reference_segments/{reference}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/samtools/faidx"


rule ivar__create_consensus_per_segment:
    input:
        bam="results/mapping/{sample}/deduplicated/{reference}.bam",
        bai="results/mapping/{sample}/deduplicated/{reference}.bam.bai",
    output:
        consensus=temp("results/consensus/{sample}/{reference}/{segment}.fa"),
    params:
        out_prefix=lambda wildcards, output: os.path.splitext(output.consensus)[0],
        samtools_params=parse_samtools_params(),
        ivar_params=parse_ivar_params(),
    log:
        "logs/ivar/create_consensus_per_segment/{sample}/{reference}/{segment}.log",
    conda:
        "../envs/ivar.yaml"
    shell:
        "("
        " samtools mpileup --no-BAQ --region {wildcards.segment} {params.samtools_params} {input.bam}"
        " |"
        " ivar consensus -p {params.out_prefix} -i {wildcards.sample}_{wildcards.segment} {params.ivar_params}"
        ") 1> {log} 2>&1"


rule concat__consensus_from_segments:
    input:
        consensuses=get_consensus_per_reference_segment,
    output:
        report(
            "results/consensus/{sample}/{reference}.fa",
            category="Consensus",
            labels={
                "Name": "{reference}",
            },
        ),
    log:
        "logs/concat/consensus_from_segments/{sample}/{reference}.log",
    conda:
        "../envs/ivar.yaml"
    shell:
        "cat {input.consensuses} 1> {output} 2> {log}"
