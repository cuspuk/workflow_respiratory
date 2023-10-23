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
        directory(os.path.join(config["reference_panel_dirpath"], "nextclade", "{name}__{accession}__{version}")),
    conda:
        "../envs/nextclade.yaml"
    log:
        os.path.join(config["reference_panel_dirpath"], "nextclade", "{name}__{accession}__{version}.log"),
    shell:
        "nextclade dataset get --name {wildcards.name} --reference {wildcards.accession} --tag {wildcards.version}"
        " --output-dir {output} 1> {log} 2>&1"


rule nextclade__run_nextclade:
    input:
        fa="results/consensus/{sample}/{reference}/{segment}.fa",
        nextclade_data=os.path.join(config["reference_panel_dirpath"], "nextclade", "{name}__{accession}__{version}"),
    output:
        "results/nextclade/{sample}/{reference}/{segment}/{name}__{accession}__{version}/nextclade.tsv",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    conda:
        "../envs/nextclade.yaml"
    log:
        "logs/nextclade/run/{sample}/{reference}/{segment}/{name}__{accession}__{version}.log",
    shell:
        "nextclade run --input-dataset={input.nextclade_data} --output-all={params.outdir} {input.fa} 1> {log} 2>&1;"


rule nextclade__merge_results_for_sample:
    input:
        nextclade_tsvs=get_nextclade_results_for_sample,
    output:
        merged_tsv="results/nextclade/{sample}/_merged/nextclade.tsv",
    conda:
        "../envs/python_pd.yaml"
    log:
        "logs/nextclade/merge_results_for_sample/{sample}.log",
    script:
        "../scripts/merge.py"


rule nextclade__merge_all_results:
    input:
        nextclade_tsvs=get_merged_nextclade_results(),
    output:
        merged_tsv="results/_aggregation/nextclade/nextclade.tsv",
    conda:
        "../envs/python_pd.yaml"
    log:
        "logs/aggregate/nextclade__merge_all_results.log",
    script:
        "../scripts/merge.py"


rule nextclade__to_html:
    input:
        "results/_aggregation/nextclade/nextclade.tsv",
    output:
        report(
            "results/_aggregation/nextclade/nextclade.html",
            category="_aggregation",
            labels={
                "Type": "Nextclade - merged all",
            },
        ),
    conda:
        "../envs/python_panel.yaml"
    log:
        "logs/aggregate/nextclade__to_html.log",
    script:
        "../scripts/nextclade_tsv_into_html.py"


rule aggregate__all_results:
    input:
        nextclade_tsv=get_nextclade_results_for_sample,
        nextclade_consensuses=get_nextclade_consensuses_for_sample,
        nextclade_refs="results/checkpoints/for_nextclade/{sample}.tsv",
        others="results/checkpoints/for_others/{sample}.tsv",
        other_results=get_others_results,
    output:
        "results/summary/{sample}.json",
    conda:
        "../envs/python.yaml"
    log:
        "logs/summary/{sample}.log",
    script:
        "../scripts/summarize_nextclade_vs_others.py"
