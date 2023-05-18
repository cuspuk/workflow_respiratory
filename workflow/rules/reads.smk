rule trimmomatic__trim_reads_pe:
    input:
        r1="results/reads/original/{sample}_R1.fastq.gz",
        r2="results/reads/original/{sample}_R2.fastq.gz",
    output:
        r1="results/reads/trimmed/{sample}_R1.fastq.gz",
        r2="results/reads/trimmed/{sample}_R2.fastq.gz",
        r1_unpaired=temp("results/reads/trimmed/{sample}_R1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/reads/trimmed/{sample}_R2.unpaired.fastq.gz"),
    log:
        "logs/trimmomatic/{sample}.log",
    params:
        trimmer=get_trimmers_from_config(),  # list
        extra="-phred33",
        compression_level="-9",
    threads: config["threads"]
    wrapper:
        "v1.31.1/bio/trimmomatic/pe"


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
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/fastqc/quality"
