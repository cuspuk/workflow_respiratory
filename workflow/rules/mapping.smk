rule bwa__build_index:
    input:
        f"{config['reference_panel_dirpath']}/{{reference}}.fa",
    output:
        idx=multiext(
            f"{config['reference_panel_dirpath']}/bwa_index/{{reference}}/{{reference}}",
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
        "logs/bwa/build_index/{reference}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/bwa/index"


rule custom__infer_and_store_read_group:
    input:
        "results/reads/%s/{sample}_R1.fastq.gz" % MAPPING_INPUT_STEP,
    output:
        read_group="results/reads/%s/read_group/{sample}.txt" % MAPPING_INPUT_STEP,
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
        index=multiext(
            f"{config['reference_panel_dirpath']}/bwa_index/{{reference}}/{{reference}}",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
        read_group="results/reads/%s/read_group/{sample}.txt" % MAPPING_INPUT_STEP,
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
        ref=f"{config['reference_panel_dirpath']}/{{reference}}.fa",
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
