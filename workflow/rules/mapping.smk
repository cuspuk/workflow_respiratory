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
        "https://github.com/xsitarcik/wrappers/raw/v1.12.2/wrappers/bwa/index"


rule custom__infer_and_store_read_group:
    input:
        get_one_fastq_file,
    output:
        read_group="results/reads/original/read_group/{sample}.txt",
    params:
        sample_id=lambda wildcards: wildcards.sample,
    localrule: True
    log:
        "logs/custom/infer_and_store_read_group/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.12.2/wrappers/custom/read_group"


rule bwa__map_reads_to_reference:
    input:
        reads=[
            "results/reads/%s/{sample}_R1.fastq.gz" % MAPPING_INPUT_STEP,
            "results/reads/%s/{sample}_R2.fastq.gz" % MAPPING_INPUT_STEP,
        ],
        index=get_bwa_index,
        read_group="results/reads/original/read_group/{sample}.txt",
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
        "https://github.com/xsitarcik/wrappers/raw/v1.12.2/wrappers/bwa/map"


rule samtools__bam_index:
    input:
        bam="results/mapping/{sample}/{step}/{reference}.bam",
    output:
        bai="results/mapping/{sample}/{step}/{reference}.bam.bai",
    threads: min(config["threads"]["bam_index"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_bam_index,
    log:
        "logs/samtools/bam_index/{sample}/{reference}_{step}.log",
    wrapper:
        "v3.3.6/bio/samtools/index"


rule picard__mark_duplicates:
    input:
        bams="results/mapping/{sample}/mapped/{reference}.bam",
        bai="results/mapping/{sample}/mapped/{reference}.bam.bai",
    output:
        bam="results/mapping/{sample}/deduplicated/{reference}.bam",
        metrics=temp("results/mapping/{sample}/deduplicated/{reference}.stats"),
    params:
        extra="--VALIDATION_STRINGENCY SILENT",
    resources:
        mem_mb=get_mem_mb_for_picard,
    log:
        "logs/picard/mark_duplicates/{sample}/{reference}.log",
    wrapper:
        "v3.3.6/bio/picard/markduplicates"


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


checkpoint nonempty_bams:
    input:
        read_counts=expand("results/mapping/{{sample}}/deduplicated/{reference}.count", reference=REFERENCES),
    output:
        report(
            "results/checkpoints/nonempty_bams/{sample}.txt",
            category="{sample}",
            labels={
                "Type": "List of non empty BAMs",
                "Reference": "-",
            },
        ),
    log:
        "logs/checkpoints/nonempty_bams/{sample}.log",
    conda:
        "../envs/python.yaml"
    localrule: True
    script:
        "../scripts/parse_read_counts_in_bam.py"


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
        "https://github.com/xsitarcik/wrappers/raw/v1.12.2/wrappers/qualimap/bamqc"


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
        "v3.3.6/bio/samtools/depth"


rule custom__depth_json:
    input:
        txt="results/mapping/{sample}/{step}/depths/{reference}.txt",
    output:
        json="results/mapping/{sample}/{step}/depths/{reference}.json",
    log:
        "logs/custom/depth_json/{sample}/{step}/{reference}.log",
    localrule: True
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/depths.py"


checkpoint mapping_quality_evaluation:
    input:
        qualimaps=get_all_qualimap_dirs,
        jsons=get_all_depths_jsons,
    output:
        passed_refs=report(
            "results/checkpoints/passed_references/{sample}.txt",
            category="{sample}",
            labels={
                "Type": "List of passed references",
                "Reference": "-",
            },
        ),
    params:
        reference_names=lambda wildcards, input: [os.path.basename(filename) for filename in input.qualimaps],
        criteria=config["consensus_params"]["reference_criteria"],
    localrule: True
    log:
        "logs/checkpoints/mapping_quality_evaluation/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/mapping_quality_evaluation.py"


rule copy__passed_bams:
    input:
        bams_and_idxes=get_deduplicated_bams_with_idxes_for_sample,
    output:
        output_dir=directory("results/checkpoints/passed_deduplicated_bams/{sample}"),
    log:
        "logs/copy/passed_bams/{sample}.log",
    conda:
        "../envs/python.yaml"
    localrule: True
    script:
        "../scripts/copy_passed_bams.py"
