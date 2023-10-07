from snakemake.utils import validate
from snakemake.io import glob_wildcards


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml")


pepfile: config["pepfile"]


validate(pep.sample_table, "../schemas/samples.schema.yaml")

## GLOBAL SPACE #################################################################

MAPPING_INPUT_STEP = "decontaminated"


def glob_references(reference_panel_dirpath: str):
    _SUFFIX = ".fa"
    _REGEX = ".*"
    location_format = os.path.join(reference_panel_dirpath, f"{{name, {_REGEX}}}{_SUFFIX}")
    return set(glob_wildcards(location_format).name)


REFERENCES = glob_references(os.path.join(config["reference_panel_dirpath"], "references"))


def get_sample_names():
    return list(pep.sample_table["sample_name"].values)


def get_one_fastq_file(wildcards):
    return pep.sample_table.loc[wildcards.sample][["fq1"]]


def get_fastq_paths(wildcards):
    return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]


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


def get_references_with_non_empty_bams(wildcards):
    with checkpoints.nonempty_bams.get(sample=wildcards.sample).output[0].open() as f:
        return f.read().splitlines()


def get_reference_fasta(wildcards):
    return os.path.join(config["reference_panel_dirpath"], "references", f"{wildcards.reference}.fa")


def get_passed_references_for_sample(sample_name: str):
    with checkpoints.mapping_quality_evaluation.get(sample=sample_name).output[0].open() as f:
        return [line.strip() for line in f.readlines()]


def get_passed_references(wildcards):
    return get_passed_references_for_sample(wildcards.sample)


def get_consensus_for_passed_references_only(wildcards):
    return expand(f"results/consensus/{wildcards.sample}/{{reference}}.fa", reference=get_passed_references(wildcards))


def get_consensuses_to_merge_for_reference(wildcards):
    return [
        f"results/consensus/{sample}/{{reference}}.fa"
        for sample in get_sample_names()
        if wildcards.reference in get_passed_references_for_sample(sample)
    ]


def get_all_aggregated_consensuses(wildcards):
    all_refs = [get_passed_references_for_sample(sample) for sample in get_sample_names()]
    all_refs_set = set([item for sublist in all_refs for item in sublist])
    return expand("results/_aggregation/{reference}.fa", reference=all_refs_set)


def get_mixed_positions_for_passed_references_only(wildcards):
    return expand(
        f"results/variants/{wildcards.sample}/{{reference}}/mixed_positions_count.tsv",
        reference=get_passed_references(wildcards),
    )


def get_variant_reports_for_passed_references_only(wildcards):
    passed_refs = get_passed_references(wildcards)
    lst1 = expand(
        f"results/variants/{wildcards.sample}/{{reference}}/mixed_positions.html",
        reference=passed_refs,
    )
    lst2 = expand(
        f"results/variants/{wildcards.sample}/{{reference}}/all.{{ext}}",
        reference=passed_refs,
        ext=["html", "vcf"],
    )
    return lst1 + lst2


def get_all_qualimap_dirs(wildcards):
    return expand(
        f"results/mapping/{wildcards.sample}/deduplicated/bamqc/{{reference}}",
        reference=get_references_with_non_empty_bams(wildcards),
    )


def get_consensus_per_reference_segment(wildcards):
    with checkpoints.index_passed_references.get(reference=wildcards.reference).output[0].open() as f:
        segments = [line.split()[0] for line in f.readlines()]
    return expand(f"results/consensus/{wildcards.sample}/{wildcards.reference}/{{segment}}.fa", segment=segments)


def get_nextclade_results(wildcards):
    results = []
    with checkpoints.select_references_for_nextclade.get(sample=wildcards.sample).output.nextclade.open() as f:
        for line in f.readlines():
            ref, seg, name, acc, version = line.split()
            results.append(f"results/nextclade/{wildcards.sample}/{ref}/{seg}/{name}__{acc}__{version}/nextclade.tsv")
    return results


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
    return {
        "fastqc_report": expand(
            "results/reads/trimmed/fastqc/{sample}_R{orientation}.html",
            sample=sample_names,
            orientation=[1, 2],
        ),
        "bams": expand(
            "results/mapping/{sample}/deduplicated/{reference}.bam",
            sample=sample_names,
            reference=REFERENCES,
        ),
        "nonempty_bams": expand(
            "results/checkpoints/nonempty_bams/{sample}.txt",
            sample=sample_names,
        ),
        "kronas": expand("results/kraken/kronas/{sample}.html", sample=sample_names),
        "consensus": expand("results/nextclade/{sample}/reference_summary.json", sample=sample_names),
        "mixed_positions": expand("results/variants/{sample}/mixed_positions_summary.txt", sample=sample_names),
        "aggregate": "results/checkpoints/aggregated_all_consensuses.txt",
    }


## PARAMETERS PARSING #################################################################


def get_cutadapt_extra() -> list[str]:
    args_lst = []
    if config["reads__trimming"].get("keep_trimmed_only", False):
        args_lst.append("--discard-untrimmed")
    if "shorten_to_length" in config["reads__trimming"]:
        args_lst.append(f"--length {config['reads__trimming']['shorten_to_length']}")
    if "cut_from_start" in config["reads__trimming"]:
        args_lst.append(f"--cut {config['reads__trimming']['cut_from_start']}")
    if "cut_from_end" in config["reads__trimming"]:
        args_lst.append(f"--cut -{config['reads__trimming']['cut_from_end']}")
    if "max_n_bases" in config["reads__trimming"]:
        args_lst.append(f"--max-n {config['reads__trimming']['max_n_bases']}")
    if "max_expected_errors" in config["reads__trimming"]:
        args_lst.append(f"--max-expected-errors {config['reads__trimming']['max_expected_errors']}")
    if param_value := config["reads__trimming"].get("anywhere_adapter", ""):
        args_lst.append(f"--anywhere file:{param_value}")
    if param_value := config["reads__trimming"].get("front_adapter", ""):
        args_lst.append(f"--front file:{param_value}")
    if param_value := config["reads__trimming"].get("regular_adapter", ""):
        args_lst.append(f"--adapter file:{param_value}")
    return args_lst


def parse_paired_cutadapt_param(pe_config, param1, param2, arg_name) -> str:
    if param1 in pe_config:
        if param2 in pe_config:
            return f"{arg_name} {pe_config[param1]}:{pe_config[param2]}"
        else:
            return f"{arg_name} {pe_config[param1]}:"
    elif param2 in pe_config:
        return f"{arg_name} :{pe_config[param2]}"
    return ""


def parse_cutadapt_comma_param(config, param1, param2, arg_name) -> str:
    if param1 in config:
        if param2 in config:
            return f"{arg_name} {config[param2]},{config[param1]}"
        else:
            return f"{arg_name} {config[param1]}"
    elif param2 in config:
        return f"{arg_name} {config[param2]},0"
    return ""


def get_cutadapt_extra_pe() -> str:
    args_lst = get_cutadapt_extra()

    cutadapt_config = config["reads__trimming"]
    if parsed_arg := parse_paired_cutadapt_param(cutadapt_config, "max_length_r1", "max_length_r2", "--maximum-length"):
        args_lst.append(parsed_arg)
    if parsed_arg := parse_paired_cutadapt_param(cutadapt_config, "min_length_r1", "min_length_r2", "--minimum-length"):
        args_lst.append(parsed_arg)
    if qual_cut_arg_r1 := parse_cutadapt_comma_param(
        cutadapt_config, "quality_cutoff_from_3_end_r1", "quality_cutoff_from_5_end_r2", "--quality-cutoff"
    ):
        args_lst.append(qual_cut_arg_r1)
    if qual_cut_arg_r2 := parse_cutadapt_comma_param(
        cutadapt_config, "quality_cutoff_from_3_end_r1", "quality_cutoff_from_5_end_r2", "-Q"
    ):
        args_lst.append(qual_cut_arg_r2)
    return " ".join(args_lst)


def get_kraken_decontamination_params():
    extra = []
    if config["reads__decontamination"]["exclude_children"]:
        extra.append("--include-children")
    if config["reads__decontamination"]["exclude_ancestors"]:
        extra.append("--include-parents")
    return " ".join(extra)


def parse_samtools_params():
    samtools_params = []
    if config["consensus_params"]["count_orphans"]:
        samtools_params.append("--count-orphans")

    samtools_params.append("--max-depth {value}".format(value=config["consensus_params"]["max_read_depth"]))
    samtools_params.append("--min-MQ {value}".format(value=config["consensus_params"]["min_mapping_quality"]))
    samtools_params.append("--min-BQ {value}".format(value=config["consensus_params"]["min_base_quality"]))
    return " ".join(samtools_params)


def parse_samtools_params_with_region(wildcards):
    return f"--region {wildcards.sample} {parse_samtools_params()}"


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


def get_mem_mb_for_trimming(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["trimming_mem_mb"] * attempt)


def get_mem_mb_for_mapping(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping_mem_mb"] * attempt)


def get_mem_mb_for_bam_index(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["bam_index_mem_mb"] * attempt)


def get_mem_mb_for_fastqc(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["fastqc_mem_mb"] * attempt)
