checkpoint select_references_for_nextclade:
    input:
        references="results/checkpoints/mapping_evaluation/{sample}.tsv",
        metadata=os.path.join(config["reference_panel_dirpath"], "nextclade_mapping.csv"),
    output:
        nextclade="results/checkpoints/for_nextclade/{sample}.tsv",
        others="results/checkpoints/for_others/{sample}.tsv",
    params:
        metadata=get_nextclade_metadata(),
        references=lambda wildcards: get_passed_references(wildcards.sample),
    log:
        "logs/checkpoints/{sample}.log",
    conda:
        "../envs/python.yaml"
    localrule: True
    script:
        "../scripts/prepare_nextclade.py"


rule nextclade__download_nextclade_dataset:
    output:
        "{prefix_dir}/nextclade/{name}/{version}/sequences.fasta",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    conda:
        "../envs/nextclade.yaml"
    wildcard_constraints:
        name="nextstrain\/.*",
        version="\d{4}-\d{2}-\d{2}--\d{2}-\d{2}-\d{2}Z",
    log:
        "{prefix_dir}/logs/{name}/download_{version}.log",
    shell:
        "nextclade dataset get --name {wildcards.name} --tag {wildcards.version}"
        " --output-dir {params.outdir} 1> {log} 2>&1"


rule nextclade__run_nextclade:
    input:
        fa="results/consensus/{sample}/{reference}/{segment}.fa",
        data=infer_relevant_nextclade_data,
    output:
        tsv="results/nextclade/{sample}/{reference}/{segment}/nextclade.tsv",
    params:
        outdir=lambda wildcards, output: os.path.dirname(os.path.realpath(output[0])),
        nextclade_dir=lambda wildcards, input: os.path.dirname(input.data),
    conda:
        "../envs/nextclade.yaml"
    log:
        "logs/nextclade/run/{sample}/{reference}/{segment}.log",
    shell:
        "nextclade run --input-dataset={params.nextclade_dir} --output-all={params.outdir} {input.fa} 1> {log} 2>&1"


rule nextclade__merge_results_for_sample:
    input:
        nextclade_tsvs=infer_nextclade_results_for_sample,
    output:
        merged_tsv="results/nextclade/{sample}/_merged/nextclade.tsv",
    params:
        names=lambda wildcards, input: [
            os.path.basename(os.path.dirname(os.path.dirname(f))) for f in input.nextclade_tsvs
        ],
    log:
        "logs/nextclade/merge_results_for_sample/{sample}.log",
    localrule: True
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.2/wrappers/nextclade/merge_tsvs"


rule nextclade__merge_all_results:
    input:
        nextclade_tsvs=get_merged_nextclade_results(),
    output:
        merged_tsv="results/_aggregation/nextclade/nextclade.tsv",
    params:
        names=[],
    log:
        "logs/aggregate/nextclade__merge_all_results.log",
    localrule: True
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.2/wrappers/nextclade/merge_tsvs"


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
    log:
        "logs/aggregate/nextclade__to_html.log",
    localrule: True
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.2/wrappers/nextclade/to_html"


rule aggregate__all_results:
    input:
        nextclade_tsv=infer_nextclade_results_for_sample,
        nextclade_consensuses=infer_nextclade_consensuses_for_sample,
        nextclade_refs="results/checkpoints/for_nextclade/{sample}.tsv",
        others="results/checkpoints/for_others/{sample}.tsv",
        other_results=infer_others_results,
    output:
        "results/summary/{sample}.json",
    conda:
        "../envs/python.yaml"
    log:
        "logs/summary/{sample}.log",
    localrule: True
    script:
        "../scripts/summarize_nextclade_vs_others.py"
