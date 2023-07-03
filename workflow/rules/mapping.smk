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
        bam=temp("results/mapping/{sample}/mapped/{reference}.bam"),
    params:
        filter="-F 4",  # discard unmapped reads
    threads: min(config["threads"]["mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping,
    log:
        "logs/bwa/map_reads_to_reference/{sample}/{reference}.log",
    benchmark:
        "benchmarks/bwa/map_reads_to_reference/{sample}/{reference}.benchmark"
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.7/wrappers/bwa/map"


rule samtools__bam_index:
    input:
        bam="results/mapping/{sample}/{step}/{reference}.bam",
    output:
        bai="results/mapping/{sample}/{step}/{reference}.bam.bai",
    benchmark:
        "benchmarks/samtools/bam_index/{step}/{reference}/{sample}.benchmark"
    threads: min(config["threads"]["bam_index"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_bam_index,
    log:
        "logs/samtools/bam_index/{sample}/{reference}_{step}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/samtools/index"


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
    benchmark:
        "benchmarks/picard/mark_duplicates/{sample}/{reference}.benchmark"
    wrapper:
        "v2.1.1/bio/picard/markduplicates"


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
            labels={
                "Name": "List of non empty BAMs",
            },
        ),
    log:
        "logs/checkpoints/nonempty_bams/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/parse_read_counts_in_bam.py"


rule qualimap__mapping_quality_report:
    input:
        bam="results/mapping/{sample}/{step}/{reference}.bam",
        bai="results/mapping/{sample}/{step}/{reference}.bam.bai",
    output:
        report(
            directory("results/mapping/{sample}/{step}/bamqc/{reference}"),
            category="Reports",
            labels={
                "Type": "Qualimap",
                "Name": "{reference}",
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
        "v2.1.1/bio/qualimap/bamqc"


checkpoint mapping_quality_evaluation:
    input:
        qualimaps=get_all_qualimap_dirs,
    output:
        passed_refs=report(
            "results/checkpoints/passed_references/{sample}.txt",
            labels={
                "Name": "List of passed references",
            },
        ),
    params:
        reference_names=lambda wildcards, input: [
            os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(filename))))
            for filename in input.qualimaps
        ],
        criteria=config["consensus_params"]["reference_criteria"],
    log:
        "logs/checkpoints/mapping_quality_evaluation/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/mapping_quality_evaluation.py"
