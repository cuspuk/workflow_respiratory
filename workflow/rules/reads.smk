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
    output:
        r1="results/reads/decontaminated/{sample}_R1.fastq.gz",
        r2="results/reads/decontaminated/{sample}_R2.fastq.gz",
    params:
        taxid=" ".join(str(taxa_id) for taxa_id in config["reads__decontamination"]["exclude_taxa_ids"]),
        children="--include-children" if config["reads__decontamination"]["exclude_children"] else "",
        parents="--include-parents" if config["reads__decontamination"]["exclude_ancestors"] else "",
    threads: config["threads"]
    log:
        "logs/kraken/decontamination/{sample}.log",
    conda:
        "../envs/krakentools.yaml"
    shell:
        "extract_kraken_reads.py --kraken {input.kraken_output} -s1 {input.r1} -s {input.r2}"
        " -o {output.r1} -o2 {output.r2} --taxid {params.taxid} --exclude"


rule fastqc__quality_report:
    input:
        read="results/reads/{step}/{fastq}.fastq.gz",
    output:
        html=report(
            "results/reads/{step}/fastqc/{fastq}.html",
            category="FastQC",
            subcategory="{step}",
        ),
        zip="results/reads/{step}/fastqc/{fastq}.zip",
        qc_data="results/reads/{step}/fastqc/{fastq}/fastqc_data.txt",
        summary_txt="results/reads/{step}/fastqc/{fastq}/summary.txt",
    log:
        "logs/fastqc/{step}/{fastq}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.3/wrappers/fastqc/quality"
