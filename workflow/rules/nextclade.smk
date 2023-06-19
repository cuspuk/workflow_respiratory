checkpoint select_references_for_nextclade:
    input:
        references="results/mapping/passed_references/{sample}.txt",
        metadata="resources/metadata.csv",
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
        directory("resources/nextclade/{reference}/{name}_{tag}"),
    conda:
        "../envs/nextclade.yaml"
    log:
        "logs/nextclade/download/{reference}/{name}_{tag}.log",
    shell:
        "nextclade dataset get --name {wildcards.name} --reference {wildcards.tag} --output-dir {output} 1> {log} 2>&1"


rule nextclade__run_nextclade:
    input:
        fa="results/consensus/{sample}/{reference}.fa",
        nextclade_data="resources/nextclade/{reference}/{name}_{tag}",
    output:
        "results/consensus/{sample}/nextclade/{reference}/{name}_{tag}/nextclade.tsv",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    conda:
        "../envs/nextclade.yaml"
    log:
        "logs/nextclade/run/{sample}/{reference}/{name}_{tag}.log",
    shell:
        "nextclade run --input-dataset={input.nextclade_data} --output-all={params.outdir} {input.fa} 1> {log} 2>&1"


rule aggregate__nextclade_results:
    input:
        nextclade_tsv=get_nextclade_results,
        nextclade_refs="results/checkpoints/for_nextclade/{sample}.tsv",
        others="results/checkpoints/for_others/{sample}.tsv",
        other_results=get_others_results,
        metadata="resources/metadata.csv",
    output:
        report(
            "results/summary/{sample}/reference_summary.json",
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
