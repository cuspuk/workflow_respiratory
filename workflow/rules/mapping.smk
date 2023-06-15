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
    params:
        filter="-F 4",  # discard unmapped reads
    threads: 2
    resources:
        mem_mb=4096,
    log:
        "logs/bwa/mapping/{reference}/{sample}.log",
    benchmark:
        "benchmarks/bwa/map_reads_to_reference/{reference}/{sample}.benchmark"
    threads: config["threads"]
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.7/wrappers/bwa/map"


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
        bam="results/mapping/{reference}/mapped/{sample}.bam",
        bai="results/mapping/{reference}/mapped/{sample}.bam.bai",
    output:
        bam="results/mapping/{reference}/deduplicated/{sample}.bam",
        stat="results/mapping/{reference}/deduplicated/{sample}.stats",
    log:
        "logs/picard/mark_duplicates/{reference}/{sample}.log",
    benchmark:
        "benchmarks/picard/mark_duplicates/{reference}/{sample}.benchmark"
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/picard/markduplicates"


rule samtools__view_number_of_reads:
    input:
        "results/mapping/{reference}/deduplicated/{sample}.bam",
    output:
        temp("results/mapping/{reference}/deduplicated/{sample}_counter.txt"),
    log:
        "logs/samtools/count_reads/{reference}/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -c {input} > {output} 2> {log}"


checkpoint nonempty_bams:
    input:
        counters=expand("results/mapping/{reference}/deduplicated/{{sample}}_counter.txt", reference=REFERENCES),
    output:
        "results/checkpoints/nonempty_bams/{sample}.tsv",
    log:
        "logs/checkpoints/nonempty_bams/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/parse_read_counts_in_bam.py"


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
