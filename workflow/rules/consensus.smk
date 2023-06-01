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


rule ivar__aggregate_consensus_from_segments:
    input:
        faidx=get_reference_faidx,
        consensus_lst=get_consensus_per_reference_segment,
    output:
        "results/consensus/{sample}/{reference}.fa",
    log:
        "logs/consensus/{sample}/{reference}.log",
    conda:
        "../envs/ivar.yaml"
    shell:
        "cat {input.consensus_lst} > {output} 2> {log}"


# def parse_samtools_params():
#     count_orphans = '--count-orphans' if config["ivar_consensus_params"]['count_orphans'] else '', # , False) else '',
#     max_depth = config["ivar_consensus_params"]['max_read_depth'], #, 1000),
#     min_mapping_quality = config["ivar_consensus_params"]['min_mapping_quality'], # 0, 0),
#     min_base_quality = config["ivar_consensus_params"]['min_base_quality'], # 13
#     no_base_alignment_quality = '--no-BAQ' if config["ivar_consensus_params"]['no_base_alignment_quality'] else '', #, False) else '',
#     ivar_quality_threshold = config["ivar_consensus_params"]['consensus_base_quality_threshold'], # 20),
#     ivar_frequency_threshold = config["ivar_consensus_params"]['consensus_frequency_threshold'], # 0),
#     ivar_min_depth = config["ivar_consensus_params"]['min_consensus_depth',] #  10),

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
