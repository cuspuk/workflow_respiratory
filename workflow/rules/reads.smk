rule custom__link_samples_to_workdir:
    input:
        r1=os.path.join(config["samples_dirpath"], "{sample}_R1.fastq.gz"),
        r2=os.path.join(config["samples_dirpath"], "{sample}_R2.fastq.gz"),
    output:
        r1="results/reads/original/{sample}_R1.fastq.gz",
        r2="results/reads/original/{sample}_R2.fastq.gz",
    params:
        dirpath=lambda w, input: os.path.splitext(output.r1)[0],
    log:
        "logs/custom/link_samples_to_workdir/{sample}.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "mkdir -p {params.dirpath} > {log} 2>&1; "
        "ln -s {input.r1} {output.r1} >> {log} 2>&1; "
        "ln -s {input.r2} {output.r2} >> {log} 2>&1; "


rule trimmomatic__trim_reads_pe:
    input:
        r1="results/reads/original/{sample}_R1.fastq.gz",
        r2="results/reads/original/{sample}_R2.fastq.gz",
    output:
        r1="results/reads/trimmed/{sample}_R1.fastq.gz",
        r2="results/reads/trimmed/{sample}_R2.fastq.gz",
        r1_unpaired=temp("results/reads/trimmed/{sample}_R1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/reads/trimmed/{sample}_R2.unpaired.fastq.gz"),
    params:
        trimmer=get_trimmers_from_config(),  # list
        extra="-phred33",
        compression_level="-9",
    threads: config["threads"]
    log:
        "logs/trimmomatic/{sample}.log",
    wrapper:
        "v1.31.1/bio/trimmomatic/pe"


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
    threads: config["threads"]
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
    threads: config["threads"]
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
    log:
        "logs/fastqc/{step}/{fastq}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.3/wrappers/fastqc/quality"
