rule bwa__build_index:
    input:
        "{reference_panel_dir}/{reference}.fa",
    output:
        idx=protected(
            multiext(
                "{reference_panel_dir}/bwa_index/{reference}/{reference}",
                ".amb",
                ".ann",
                ".bwt",
                ".pac",
                ".sa",
            )
        ),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output.idx[0])[0],
        approach="bwtsw",
    log:
        "{reference_panel_dir}/bwa_index/logs/{reference}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/bwa/index"


rule custom__infer_read_group:
    input:
        infer_fastq_for_mapping,
    output:
        read_group="results/reads/.read_groups/{sample}.txt",
    params:
        sample_id=lambda wildcards: wildcards.sample,
    log:
        "logs/custom/infer_and_store_read_group/{sample}.log",
    localrule: True
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.14.3/wrappers/custom/read_group"


rule bwa__map_reads_to_reference:
    input:
        reads=infer_fastq_for_mapping,
        index=infer_bwa_index,
        read_group="results/reads/.read_groups/{sample}.txt",
    output:
        bam=temp("results/mapping/{sample}/mapped/{reference}.bam"),
    params:
        filter="-F 4",  # discard unmapped reads
    threads: min(config["threads"]["mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping,
    log:
        "logs/bwa/map_reads_to_reference/{sample}/{reference}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/bwa/map"


rule samtools__bam_index:
    input:
        bam="results/mapping/{sample}/mapped/{reference}.bam",
    output:
        bai="results/mapping/{sample}/mapped/{reference}.bam.bai",
    threads: min(config["threads"]["bam_index"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_bam_index,
    log:
        "logs/samtools/bam_index/mapped/{sample}/{reference}.log",
    wrapper:
        "v4.7.0/bio/samtools/index"


rule picard__mark_duplicates:
    input:
        bams="results/mapping/{sample}/mapped/{reference}.bam",
        bai="results/mapping/{sample}/mapped/{reference}.bam.bai",
    output:
        bam="results/mapping/{sample}/deduplicated/{reference}.bam",
        idx="results/mapping/{sample}/deduplicated/{reference}.bam.bai",
        metrics=temp("results/mapping/{sample}/deduplicated/{reference}.stats"),
    resources:
        mem_mb=get_mem_mb_for_picard,
    log:
        "logs/picard/mark_duplicates/{sample}/{reference}.log",
    wrapper:
        "v4.0.0/bio/picard/markduplicates"


rule samtools__view_number_of_reads:
    input:
        "results/mapping/{sample}/deduplicated/{reference}.bam",
    output:
        temp("results/mapping/{sample}/deduplicated/{reference}.count"),
    log:
        "logs/samtools/view_number_of_reads/{reference}/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -c {input} > {output} 2> {log}"


checkpoint checkpoint_get_all_nonempty_bams:
    input:
        read_counts=expand("results/mapping/{{sample}}/deduplicated/{reference}.count", reference=REFERENCES),
    output:
        report(
            "results/checkpoints/mapped_reads_per_reference/{sample}.tsv",
            category="{sample}",
            labels={
                "Type": "BAMs with read counts",
                "Reference": "-",
            },
        ),
    log:
        "logs/checkpoints/mapped_reads_per_reference/{sample}.log",
    conda:
        "../envs/python.yaml"
    localrule: True
    script:
        "../scripts/qc_summarize_mapped_reads.py"


rule qualimap__mapping_quality_report:
    input:
        bam="results/mapping/{sample}/{step}/{reference}.bam",
        bai="results/mapping/{sample}/{step}/{reference}.bam.bai",
    output:
        report_dir=report(
            directory("results/mapping/{sample}/{step}/bamqc/{reference}"),
            category="{sample}",
            labels={
                "Type": "Qualimap - {step}",
                "Reference": "{reference}",
            },
            htmlindex="qualimapReport.html",
        ),
    params:
        extra=[
            "--paint-chromosome-limits",
            "-outformat PDF:HTML",
        ],
    resources:
        mem_mb=get_mem_mb_for_qualimap,
    log:
        "logs/qualimap/mapping_quality_report/{sample}/{step}/{reference}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/qualimap/bamqc"


rule samtools__depth:
    input:
        bams=["results/mapping/{sample}/{step}/{reference}.bam"],
        bai="results/mapping/{sample}/{step}/{reference}.bam.bai",
    output:
        temp("results/mapping/{sample}/{step}/depths/{reference}.txt"),
    log:
        "logs/samtools/depth/{sample}/{step}/{reference}.log",
    params:
        extra="-a",
    wrapper:
        "v4.0.0/bio/samtools/depth"


checkpoint checkpoint_mapping_evaluation:
    input:
        depths=infer_depths_for_nonempty_bams,
        qualimaps=infer_qualimaps_for_nonempty_bams,
    output:
        tsv=report(
            "results/checkpoints/mapping_evaluation/{sample}.tsv",
            category="{sample}",
            labels={
                "Type": "List of passed references",
                "Reference": "-",
            },
        ),
        json="results/checkpoints/mapping_evaluation/{sample}.json",
    params:
        reference_names=lambda wildcards, input: [
            os.path.basename(os.path.splitext(filename)[0]) for filename in input.depths
        ],
        threshold=config["consensus_params"]["reference_criteria"]["min_genome_fraction_with_10x_coverage"],
    localrule: True
    log:
        "logs/checkpoints/mapping_evaluation/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/qc_mapping_evaluation.py"


rule copy__passed_bams:
    input:
        bams_and_idxes=infer_passed_bams_and_bais,
    output:
        output_dir=directory("results/checkpoints/passed_deduplicated_bams/{sample}"),
    log:
        "logs/copy/passed_bams/{sample}.log",
    conda:
        "../envs/python.yaml"
    localrule: True
    script:
        "../scripts/copy_passed_bams.py"
