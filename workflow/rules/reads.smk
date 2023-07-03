rule cutadapt__trim_reads_pe:
    input:
        r1="results/reads/original/{sample}_R1.fastq.gz",
        r2="results/reads/original/{sample}_R2.fastq.gz",
    output:
        r1=temp("results/reads/trimmed/{sample}_R1.fastq.gz"),
        r2=temp("results/reads/trimmed/{sample}_R2.fastq.gz"),
        report="results/reads/trimmed/{sample}_cutadapt.json",
    params:
        overlap=config["reads__trimming"]["overlap"],
        error_rate=config["reads__trimming"]["error_rate"],
        times=config["reads__trimming"]["times"],
        action=config["reads__trimming"]["action"],
        extra=get_cutadapt_extra_pe(),
    resources:
        mem_mb=get_mem_mb_for_trimming,
    threads: min(config["threads"]["trimming"], config["max_threads"])
    log:
        "logs/cutadapt/trim_reads_pe/{sample}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.1/wrappers/cutadapt/paired"


rule kraken__decontaminate:
    input:
        r1="results/reads/trimmed/{sample}_R1.fastq.gz",
        r2="results/reads/trimmed/{sample}_R2.fastq.gz",
        kraken_output="results/kraken/{sample}.kraken",
        kraken_report="results/kraken/{sample}.kreport2",
    output:
        r1=temp("results/reads/decontaminated/{sample}_R1.fastq"),
        r2=temp("results/reads/decontaminated/{sample}_R2.fastq"),
    params:
        taxid=" ".join(str(taxa_id) for taxa_id in config["reads__decontamination"]["exclude_taxa_ids"]),
        children="--include-children" if config["reads__decontamination"]["exclude_children"] else "",
        parents="--include-parents" if config["reads__decontamination"]["exclude_ancestors"] else "",
    threads: min(config["threads"]["kraken"], config["max_threads"])
    log:
        "logs/kraken/decontaminate/{sample}.log",
    conda:
        "../envs/krakentools.yaml"
    shell:
        "extract_kraken_reads.py -k {input.kraken_output} -r {input.kraken_report} -s {input.r1} -s2 {input.r2}"
        " -o {output.r1} -o2 {output.r2} -t {params.taxid} --exclude --fastq-output > {log} 2>&1"


rule pigz__compress_decontaminated:
    input:
        "results/reads/decontaminated/{sample}_{strand}.fastq",
    output:
        "results/reads/decontaminated/{sample}_{strand}.fastq.gz",
    threads: 1
    params:
        level=6,
    log:
        "logs/pigz/compress/{sample}_{strand}.log",
    conda:
        "../envs/pigz.yaml"
    shell:
        "pigz {input} -{params.level} -c -p {threads} > {output} 2> {log}"


rule fastqc__quality_report:
    input:
        read="results/reads/{step}/{fastq}.fastq.gz",
    output:
        html=report(
            "results/reads/{step}/fastqc/{fastq}.html",
            category="Reports",
            labels={
                "Type": "Fastqc",
                "Name": "{fastq}",
            },
        ),
        zip="results/reads/{step}/fastqc/{fastq}.zip",
        qc_data="results/reads/{step}/fastqc/{fastq}/fastqc_data.txt",
        summary_txt="results/reads/{step}/fastqc/{fastq}/summary.txt",
    threads: min(config["threads"]["fastqc"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_fastqc,
    log:
        "logs/fastqc/{step}/{fastq}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.3/wrappers/fastqc/quality"
