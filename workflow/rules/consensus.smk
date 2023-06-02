checkpoint mapping_quality_evaluation:
    input:
        get_all_qualimap_dirs,
    output:
        passed_refs="results/mapping/passed_references/{sample}.txt",
    params:
        min_mean_coverage=config["consensus_params"]["reference_criteria"]["min_mean_coverage"],
    log:
        "logs/checkpoints/mapping_quality_evaluation/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/mapping_quality_evaluation.py"


rule ivar__create_consensus_per_segment:
    input:
        bam="results/mapping/{reference}/deduplicated/{sample}.bam",
        bai="results/mapping/{reference}/deduplicated/{sample}.bam.bai",
    output:
        consensus=temp("results/consensus/{sample}/{reference}/{segment}.fa"),
    params:
        out_prefix=lambda wildcards, output: os.path.splitext(output.idx[0])[0],
        samtools_params=parse_samtools_params(),
        ivar_params=parse_ivar_params(),
    log:
        "logs/consensus/{sample}/{reference}/{segment}.log",
    conda:
        "../envs/ivar.yaml"
    shell:
        "(samtools mpileup --region {wildcards.segment} {params.samtools_params} {input.bam}"
        " | ivar consensus -p {params.out_prefix} -i {wildcards.sample}_{wildcards.segment} {params.ivar_params}) > {log}"


checkpoint faidx_reference_segments:
    input:
        reference=get_reference_fasta,
    output:
        f"{config['reference_panel_dirpath']}/{{reference}}.fa.fai",
    log:
        "logs/checkpoints/reference_segments/{reference}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/samtools/faidx"


rule ivar__aggregate_consensus_from_segments:
    input:
        get_consensus_per_reference_segment,
    output:
        "results/consensus/{sample}/{reference}.fa",
    log:
        "logs/consensus/{sample}/{reference}.log",
    conda:
        "../envs/ivar.yaml"
    shell:
        "cat {input.consensus_lst} > {output} 2> {log}"


# rule ivar__create_consensus_fasta_segment:
#     input:
#         bam="results/mapping/{reference}/deduplicated/{sample}.bam",
#         bai="results/mapping/{reference}/deduplicated/{sample}.bam.bai",
#         fai   = get_reference_faidx,
#     output:
#         consensus = 'results/consensus/{sample}/{reference}.fa'
#     params:
#         out_prefix = 'results/consensus/{sample}/{reference}'
#         samtools_params = parse_samtools_params(),
#         ivar_params = parse_ivar_params(),
#     log:
#         out = 'consensus/{reference}-{panel}/log/{sample}.fa.out',
#         err = 'consensus/{reference}-{panel}/log/{sample}.fa.err'
#     conda:
#         '../envs/ivar.yaml'
#     shell:
#         "rm {log}; "
#         "for REF in `cut -f 1 {input.fai}`; do"
#         " samtools mpileup --region $REF {params.samtools_params} {input.bam} 2>> {log}"
#         " |"
#         " ivar consensus -p {params.out_prefix}_$REF -i {wildcards.sample}_$REF {params.ivar} 1>> {log} 2>>&1"
#         " OUT_FILE='{params.out_prefix}_$REF'.fa;"
#         " QUAL_FILE='{params.out_prefix}_$REF'.qual.txt;"
#         " cat $OUT_FILE >> {output.fasta} 2>> {log.err};"
#         " rm $OUT_FILE; rm $QUAL_FILE;"
#         " done 2>> {log.err};"


rule aggregate_consensus:
    input:
        get_consensus_for_passed_references_only,
    output:
        "results/summary/{sample}/aggr_result.txt",
    log:
        "logs/test/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        "touch {output} > {log} 2>&1"
