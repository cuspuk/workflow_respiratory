rule bwa__build_index:
    input:
        "{reference_panel_dir}/{reference}.fa",
    output:
        idx=multiext(
            "{reference_panel_dir}/bwa_index/{reference}/{reference}",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output.idx[0])[0],
        approach="bwtsw",
    log:
        "{reference_panel_dir}/bwa_index/logs/{reference}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/bwa/index"


rule custom__infer_and_store_read_group:
    input:
        "results/reads/original/{sample}_R1.fastq.gz",
    output:
        read_group="results/reads/original/read_group/{sample}.txt",
    log:
        "logs/custom/infer_and_store_read_group/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/save_read_group.py"


rule bwa__map_reads_to_reference:
    input:
        reads=[
            "results/reads/%s/{sample}_R1.fastq.gz" % MAPPING_INPUT_STEP,
            "results/reads/%s/{sample}_R2.fastq.gz" % MAPPING_INPUT_STEP,
        ],
        index=get_bwa_index,
        read_group="results/reads/original/read_group/{sample}.txt",
    output:
        bam=temp("results/mapping/{reference}/mapped/{sample}.bam"),
    log:
        "logs/bwa/mapping/{reference}/{sample}.log",
    benchmark:
        "benchmarks/bwa/map_reads_to_reference/{reference}/{sample}.benchmark"
    threads: config["threads"]
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/bwa/map"


rule samtools__sort_mapped_reads:
    input:
        ref=get_reference_fasta,
        bam="results/mapping/{reference}/mapped/{sample}.bam",
    output:
        bam=temp("results/mapping/{reference}/sorted/{sample}.bam"),
    log:
        "logs/samtools/sorting/{reference}/{sample}.log",
    benchmark:
        "benchmarks/samtools/sort_mapped_reads/{reference}/{sample}.benchmark"
    threads: config["threads"]
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/samtools/sort"


rule samtools__bam_index:
    input:
        bam="results/mapping/{reference}/{step}/{sample}.bam",
    output:
        bai="results/mapping/{reference}/{step}/{sample}.bam.bai",
    benchmark:
        "benchmarks/samtools/indexing/{step}/{reference}/{sample}.benchmark"
    threads: config["threads"]
    log:
        "logs/samtools/bam_index/{step}/{reference}/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/samtools/index"


rule picard__mark_duplicates:
    input:
        bam="results/mapping/{reference}/sorted/{sample}.bam",
        bai="results/mapping/{reference}/sorted/{sample}.bam.bai",
    output:
        bam="results/mapping/{reference}/deduplicated/{sample}.bam",
        stat="results/mapping/{reference}/deduplicated/{sample}.stats",
    log:
        "logs/picard/mark_duplicates/{reference}/{sample}.log",
    benchmark:
        "benchmarks/picard/mark_duplicates/{reference}/{sample}.benchmark"
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/picard/markduplicates"


rule qualimap__mapping_quality_report:
    input:
        bam="results/mapping/{reference}/{step}/{sample}.bam",
        bai="results/mapping/{reference}/{step}/{sample}.bam.bai",
    output:
        report_dir=report(
            directory("results/mapping/{reference}/{step}/bamqc/{sample}"),
            category="BamQC",
            subcategory="{step}",
            htmlindex="qualimapReport.html",
        ),
    params:
        extra=[
            "--paint-chromosome-limits",
            "-outformat PDF:HTML",
        ],
    resources:
        mem_mb=10000,
    log:
        "logs/qualimap/mapping_quality_report/{reference}/{step}/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/qualimap/bamqc"


checkpoint mapping_quality_evaluation:
    input:
        get_all_qualimap_dirs,
    output:
        passed_refs="results/mapping/passed_references/{sample}.txt",
    params:
        min_mean_coverage=50,
    log:
        "logs/checkpoints/mapping_quality_evaluation/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/mapping_quality_evaluation.py"


rule fake__create_consensus:
    input:
        bam="results/mapping/{reference}/deduplicated/{sample}.bam",
        bai="results/mapping/{reference}/deduplicated/{sample}.bam.bai",
        fasta=get_reference_fasta,
    output:
        consensus="results/consensus/{sample}/{reference}.fa",
    log:
        "logs/consensus/{sample}/{reference}.log",
    conda:
        "../envs/python.yaml"
    shell:
        "cat {wildcards.reference} > {output} 2>&1"


# rule ivar__create_consensus_fasta:
#     input:
#         bam="results/mapping/{reference}/deduplicated/{sample}.bam",
#         bai="results/mapping/{reference}/deduplicated/{sample}.bam.bai",
#         fasta   = get_reference_fasta,
#         fai   = get_reference_faidx,
#     output:
#         consensus = 'results/consensus/{sample}/{reference}.fa'
#     params:
#         out_prefix = 'results/consensus/{sample}/{reference}'
#         # count_orphans = '--count-orphans' if method_config.get('count_orphans', False) else '',
#         max_depth = method_config.get('max_depth', 8000),
#         min_mapping_quality = method_config.get('min_mapping_quality', 0),
#         min_base_quality = method_config.get('min_base_quality', 13),
#         no_base_alignment_quality = '--no-BAQ' if method_config.get('no_base_alignment_quality', False) else '',
#         ivar_quality_threshold = method_config.get('ivar_quality_threshold', 20),
#         ivar_frequency_threshold = method_config.get('ivar_frequency_threshold', 0),
#         ivar_min_depth = method_config.get('ivar_min_depth', 10),
#     log:
#         out = 'consensus/{reference}-{panel}/log/{sample}.fa.out',
#         err = 'consensus/{reference}-{panel}/log/{sample}.fa.err'
#     conda:
#         '../envs/ivar_consensus.yaml'
#     shell:
#         '''
#         for REF in `cut -f 1 {input.fai}`; do \
#             samtools mpileup \
#                 --region $REF \
#                 {params.no_base_alignment_quality} \
#                 -aa \
#                 {params.count_orphans} \
#                 --max-depth {params.max_depth} \
#                 --min-MQ {params.min_mapping_quality} \
#                 --min-BQ {params.min_base_quality} \
#                 {input.bam} \
#                 2> {log.err} \
#             | ivar consensus \
#                 -p {params.out_prefix}_$REF \
#                 -i {wildcards.sample}_$REF \
#                 -q {params.ivar_quality_threshold} \
#                 -t {params.ivar_frequency_threshold} \
#                 -m {params.ivar_min_depth} \
#                 1> {log.out} \
#                 2>> {log.err}; \
#             OUT_FILE="{params.out_prefix}_$REF".fa; \
#             QUAL_FILE="{params.out_prefix}_$REF".qual.txt; \
#             cat $OUT_FILE \
#                 >> {output.fasta} \
#                 2>> {log.err};
#             rm $OUT_FILE; \
#             rm $QUAL_FILE; \
#         done \
#             2>> {log.err};
#         '''


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
