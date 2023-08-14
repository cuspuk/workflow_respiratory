checkpoint select_references_for_nextclade:
    input:
        references="results/checkpoints/passed_references/{sample}.txt",
        metadata=os.path.join(config["reference_panel_dirpath"], "nextclade_mapping.csv"),
    output:
        nextclade="results/checkpoints/for_nextclade/{sample}.tsv",
        others="results/checkpoints/for_others/{sample}.tsv",
    log:
        "logs/checkpoints/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/prepare_nextclade.py"


rule nextclade__download_nextclade_dataset:
    output:
        directory("resources/nextclade/{reference}/{name}__{accession}__{version}"),
    conda:
        "../envs/nextclade.yaml"
    log:
        "logs/nextclade/download/{reference}/{name}__{accession}__{version}.log",
    shell:
        "nextclade dataset get --name {wildcards.name} --reference {wildcards.accession} --tag {wildcards.version}"
        " --output-dir {output} 1> {log} 2>&1"


rule nextclade__run_nextclade:
    input:
        fa="results/consensus/{sample}/{reference}/{segment}.fa",
        nextclade_data="resources/nextclade/{reference}/{name}__{accession}__{version}",
    output:
        "results/consensus/{sample}/nextclade/{reference}/{segment}/{name}__{accession}__{version}/nextclade.tsv",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    conda:
        "../envs/nextclade.yaml"
    log:
        "logs/nextclade/run/{sample}/{reference}/{segment}/{name}__{accession}__{version}.log",
    shell:
        "nextclade run --input-dataset={input.nextclade_data} --output-all={params.outdir} {input.fa} 1> {log} 2>&1;"


rule aggregate__nextclade_results:
    input:
        nextclade_tsv=get_nextclade_results,
        nextclade_refs="results/checkpoints/for_nextclade/{sample}.tsv",
        others="results/checkpoints/for_others/{sample}.tsv",
        other_results=get_others_results,
    output:
        report(
            "results/consensus/{sample}/nextclade/reference_summary.json",
            labels={
                "Name": "Consensus summary",
            },
        ),
    conda:
        "../envs/python.yaml"
    log:
        "logs/summary/{sample}.log",
    script:
        "../scripts/summarize_nextclade_vs_others.py"
