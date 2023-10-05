checkpoint index_passed_references:
    input:
        reference=get_reference_fasta,
    output:
        protected("{reference_dir}/{reference}.fa.fai"),
    log:
        "{reference_dir}/logs/samtools__prepare_fai_index/{reference}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/samtools/faidx"


rule ivar__create_consensus_per_segment:
    input:
        bam="results/mapping/{sample}/deduplicated/{reference}.bam",
        bai="results/mapping/{sample}/deduplicated/{reference}.bam.bai",
    output:
        consensus="results/consensus/{sample}/{reference}/{segment}.fa",
    params:
        samtools_params=parse_samtools_params_with_region,
        ivar_params=parse_ivar_params(),
    log:
        "logs/ivar/create_consensus_per_segment/{sample}/{reference}/{segment}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.6.0/wrappers/ivar/consensus"


rule concat__consensus_from_segments:
    input:
        consensuses=get_consensus_per_reference_segment,
    output:
        report(
            "results/consensus/{sample}/{reference}.fa",
            category="{sample}",
            labels={
                "Type": "Consensus",
                "Reference": "{reference}",
            },
        ),
    log:
        "logs/concat/consensus_from_segments/{sample}/{reference}.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input.consensuses} 1> {output} 2> {log}"


rule concat__consensuses_for_references:
    input:
        consensuses=get_consensuses_to_merge_for_reference,
    output:
        report(
            "results/_aggregation/{reference}.fa",
            category="_aggregation",
            labels={
                "Type": "Consensus for {reference}",
            },
        ),
    log:
        "logs/concat/consensuses_for_references/{reference}.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "cat {input.consensuses} 1> {output} 2> {log}"


rule aggregate__all_consensuses:
    input:
        consensuses=get_all_aggregated_consensuses,
    output:
        touch("results/checkpoints/aggregated_all_consensuses.txt"),
    log:
        "logs/aggregate/all_consensuses.log",
