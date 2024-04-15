from snakemake.utils import validate
from snakemake.io import glob_wildcards
from functools import cache


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml")


### Layer for adapting other workflows  ###############################################################################


def get_fastq_for_mapping(wildcards):
    return reads_workflow.get_final_fastq_for_sample(wildcards.sample)


def get_sample_names():
    return reads_workflow.get_sample_names()


## GLOBAL SPACE #################################################################


def glob_references(reference_panel_dirpath: str):
    _SUFFIX = ".fa"
    _REGEX = ".*"
    location_format = os.path.join(reference_panel_dirpath, f"{{name, {_REGEX}}}{_SUFFIX}")
    return set(glob_wildcards(location_format).name)


def get_reference_dir():
    return os.path.join(config["reference_panel_dirpath"], "references")


REFERENCES = glob_references(get_reference_dir())


def get_constraints():
    return {
        "sample": "|".join(get_sample_names()),
        "reference": "|".join(REFERENCES),
    }


#### COMMON STUFF #################################################################


def get_bwa_index(wildcards):
    base_dir = os.path.join(
        config["reference_panel_dirpath"], "references", "bwa_index", wildcards.reference, wildcards.reference
    )
    return multiext(
        base_dir,
        ".amb",
        ".ann",
        ".bwt",
        ".pac",
        ".sa",
    )


def get_reference_fasta(wildcards):
    return os.path.join(config["reference_panel_dirpath"], "references", f"{wildcards.reference}.fa")


@cache
def get_references_with_non_empty_bams(sample_name: str):
    # checkpoint produces tsv of 3 values: PASS/FAIL, reference, count
    with checkpoints.checkpoint_get_all_nonempty_bams.get(sample=sample_name).output[0].open() as f:
        rows: list[tuple[str, str, int]] = [row.strip().split("\t") for row in f.readlines()]
        refs = [row[1] for row in rows if row[0] == "PASS"]
        return refs


@cache
def get_passed_references(sample_name: str):
    # checkpoint produces tsv of 3 values: PASS/FAIL, reference, average_coverage
    with checkpoints.checkpoint_mapping_evaluation.get(sample=sample_name).output.tsv.open() as f:
        rows: list[tuple[str, str, float]] = [row.strip().split("\t") for row in f.readlines()]
        refs = [row[1] for row in rows if row[0] == "PASS"]
        return refs


def infer_passed_bams_and_bais(wildcards):
    return expand(
        f"results/mapping/{wildcards.sample}/deduplicated/{{reference}}.{{ext}}",
        ext=["bam", "bam.bai"],
        reference=get_passed_references(wildcards.sample),
    )


def get_consensuses_to_merge_for_reference(wildcards):
    return [
        f"results/consensus/{sample}/{{reference}}.fa"
        for sample in get_sample_names()
        if wildcards.reference in get_passed_references(sample)
    ]


def get_all_aggregated_consensuses(wildcards):
    all_refs = [get_passed_references(sample) for sample in get_sample_names()]
    all_refs_set = set([item for sublist in all_refs for item in sublist])
    return expand("results/_aggregation/consensus/{reference}.fa", reference=all_refs_set)


def get_mixed_positions_for_passed_references_only(wildcards):
    return expand(
        f"results/variants/{wildcards.sample}/{{reference}}/mixed_positions_count.tsv",
        reference=get_passed_references(wildcards.sample),
    )


def get_variant_reports_for_passed_references_only(wildcards):
    return expand(
        f"results/variants/{wildcards.sample}/{{reference}}/{{ext}}",
        ext=["mixed_positions.html", "all.html", "all.vcf"],
        reference=get_passed_references(wildcards.sample),
    )


# def get_all_qualimap_dirs(wildcards):
#     return expand(
#         f"results/mapping/{wildcards.sample}/deduplicated/bamqc/{{reference}}",
#         reference=get_references_with_non_empty_bams(wildcards),
#     )


def get_depths_for_nonempty_bams(wildcards):
    return expand(
        f"results/mapping/{wildcards.sample}/deduplicated/depths/{{reference}}.txt",
        reference=get_references_with_non_empty_bams(wildcards),
    )


def get_consensus_per_reference_segment(wildcards):
    with checkpoints.index_passed_references.get(
        reference_dir=get_reference_dir(), reference=wildcards.reference
    ).output[0].open() as f:
        segments = [line.split()[0] for line in f.readlines()]
    return expand(f"results/consensus/{wildcards.sample}/{wildcards.reference}/{{segment}}.fa", segment=segments)


def get_nextclade_results_for_sample(wildcards):
    results = []
    with checkpoints.select_references_for_nextclade.get(sample=wildcards.sample).output.nextclade.open() as f:
        for line in f.readlines():
            ref, seg, name, version = line.split()
            results.append(f"results/nextclade/{wildcards.sample}/{ref}/{seg}/nextclade.tsv")
    return results


def infer_relevant_nextclade_data(wildcards):
    with checkpoints.select_references_for_nextclade.get(sample=wildcards.sample).output.nextclade.open() as f:
        for line in f.readlines():
            ref, seg, name, version = line.split()
            if wildcards.reference == ref and wildcards.segment == seg:
                x = os.path.realpath(config["reference_panel_dirpath"])
                return os.path.join(x, f"nextclade/{name}/{version}/sequences.fasta")


def get_nextclade_consensuses_for_sample(wildcards):
    results = []
    with checkpoints.select_references_for_nextclade.get(sample=wildcards.sample).output.nextclade.open() as f:
        for line in f.readlines():
            ref, seg, name, version = line.split()
            results.append(f"results/consensus/{wildcards.sample}/{ref}.fa")
    return results


def get_merged_nextclade_results():
    return expand("results/nextclade/{sample}/_merged/nextclade.tsv", sample=get_sample_names())


def get_others_results(wildcards):
    with checkpoints.select_references_for_nextclade.get(sample=wildcards.sample).output.others.open() as f:
        references = [line.strip() for line in f.readlines()]
    if references:
        return expand(f"results/consensus/{wildcards.sample}/{{reference}}.fa", reference=references)
    else:
        return []


#### OUTPUTS #################################################################


def get_outputs():
    sample_names = get_sample_names()
    outputs = {
        "passed_bams": expand(
            "results/checkpoints/passed_deduplicated_bams/{sample}",
            sample=sample_names,
        ),
        "nonempty_bams": expand(
            "results/checkpoints/mapped_reads/{sample}.tsv",
            sample=sample_names,
        ),
        "consensus": expand("results/summary/{sample}.json", sample=sample_names),
        "mixed_positions": expand("results/variants/{sample}/mixed_positions_summary.txt", sample=sample_names),
        "merged_nextclades": expand("results/nextclade/{sample}/_merged/nextclade.tsv", sample=sample_names),
    }
    if len(sample_names) > 1:
        outputs["aggregate_consensus"] = "results/checkpoints/aggregated_all_consensuses.txt"
        outputs["aggregate_nextclades"] = "results/_aggregation/nextclade/nextclade.html"
    return outputs


## PARAMETERS PARSING #################################################################


def parse_samtools_params():
    samtools_params = []
    if config["consensus_params"]["count_orphans"]:
        samtools_params.append("--count-orphans")

    samtools_params.append("--max-depth {value}".format(value=config["consensus_params"]["max_read_depth"]))
    samtools_params.append("--min-MQ {value}".format(value=config["consensus_params"]["min_mapping_quality"]))
    samtools_params.append("--min-BQ {value}".format(value=config["consensus_params"]["min_base_quality"]))
    return " ".join(samtools_params)


def parse_samtools_params_with_region(wildcards):
    return f"--region {wildcards.segment} {parse_samtools_params()}"


def parse_ivar_params():
    ivar_params = []

    ivar_params.append("-q {value}".format(value=config["consensus_params"]["consensus_base_quality_threshold"]))
    ivar_params.append("-t {value}".format(value=config["consensus_params"]["consensus_frequency_threshold"]))
    ivar_params.append("-m {value}".format(value=config["consensus_params"]["min_consensus_depth"]))

    return " ".join(ivar_params)


def parse_samtools_params_for_variants():
    samtools_params = []
    if config["mixed_positions_params"]["count_orphans"]:
        samtools_params.append("--count-orphans")

    samtools_params.append("--max-depth {value}".format(value=config["mixed_positions_params"]["max_read_depth"]))
    samtools_params.append("--min-MQ {value}".format(value=config["mixed_positions_params"]["min_mapping_quality"]))
    samtools_params.append("--min-BQ {value}".format(value=config["mixed_positions_params"]["min_base_quality"]))
    return " ".join(samtools_params)


def parse_ivar_params_for_variants():
    ivar_params = []

    ivar_params.append("-q {value}".format(value=config["mixed_positions_params"]["min_base_quality_threshold"]))
    ivar_params.append("-t {value}".format(value=config["mixed_positions_params"]["min_frequency_threshold"]))
    ivar_params.append("-m {value}".format(value=config["mixed_positions_params"]["min_read_depth"]))
    return " ".join(ivar_params)


### RESOURCES


def get_mem_mb_for_picard(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["picard_mem_mb"] * attempt)


def get_mem_mb_for_qualimap(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["qualimap_mem_mb"] * attempt)


def get_mem_mb_for_mapping(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping_mem_mb"] * attempt)


def get_mem_mb_for_bam_index(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["bam_index_mem_mb"] * attempt)
