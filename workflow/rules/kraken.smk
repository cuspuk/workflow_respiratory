rule curl__download_kraken_db:
    output:
        directory("resources/kraken/{tag}"),
    params:
        url=lambda wildcards, output: "https://genome-idx.s3.amazonaws.com/kraken/{tag}.tar.gz".format(
            tag=os.path.basename(output[0])
        ),
    retries: 3
    log:
        "logs/curl/download_kraken_db/{tag}.log",
    conda:
        "../envs/curl.yaml"
    shell:
        "mkdir -p {output} && curl -SL {params.url} | tar zxvf - -C {output} 2> {log}"


rule kraken__analysis:
    input:
        db="resources/kraken/{tag}".format(tag=config["kraken_db"]),
        r1="results/reads/trimmed/{sample}_R1.fastq.gz",
        r2="results/reads/trimmed/{sample}_R2.fastq.gz",
    output:
        kraken_output=temp("results/kraken/{sample}.kraken"),
        report="results/kraken/{sample}.kreport2",
    threads: config["threads"]
    params:
        extra="--gzip-compressed",
    log:
        "logs/kraken/analysis/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    shell:
        "(kraken2 --db {input.db} --threads {threads} --paired {params.extra}"
        " --memory-mapping --report {output.report} {input.r1} {input.r2} 1> {output.kraken_output}) 2> {log}"


rule krona__update_taxonomy:
    output:
        directory("resources/krona/"),
    log:
        "logs/krona/update_taxonomy.log",
    conda:
        "../envs/kraken2.yaml"
    shell:
        "ktUpdateTaxonomy.sh {output} 1> {log} 2>&1"


rule krona__create_chart:
    input:
        kraken_output="results/kraken/{sample}.kreport2",
        tax_db="resources/krona/",
    output:
        report(
            "results/kraken/kronas/{sample}.html",
            category="Kronas",
        ),
    log:
        "logs/krona/create_chart/{sample}.log",
    params:
        extra="-m 3 -t 5",
    conda:
        "../envs/kraken2.yaml"
    shell:
        "ktImportTaxonomy {params.extra} -tax {input.tax_db} -o {output} {input.kraken_output} 1> {log} 2>&1"
