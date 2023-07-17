from snakemake.utils import validate
from snakemake.io import glob_wildcards


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml")


## GLOBAL SPACE #################################################################

MAPPING_INPUT_STEP = "decontaminated"


def glob_samples(regex: str):
    _PAIRED_INFO = ("_R1", "_R2")
    _BASE_LOCATION = "results/reads/original/{{sample, {regex}}}{paired_str}.fastq.gz"
    return set(glob_wildcards(_BASE_LOCATION.format(regex=regex, paired_str=_PAIRED_INFO[0])).sample)


def glob_references(reference_panel_dirpath: str):
    _SUFFIX = ".fa"
    _REGEX = ".*"
    location_format = os.path.join(reference_panel_dirpath, f"{{name, {_REGEX}}}{_SUFFIX}")
    return set(glob_wildcards(location_format).name)


SAMPLES = glob_samples(config["sample_names_regex"])
REFERENCES = glob_references(os.path.join(config["reference_panel_dirpath"], "references"))


def get_constraints():
    return {
        "sample": "|".join(SAMPLES),
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


def get_reference_faidx(wildcards):
    return os.path.join(config["reference_panel_dirpath"], "references", f"{wildcards.reference}.fa.fai")


def get_passed_references(wildcards):
    with checkpoints.mapping_quality_evaluation.get(sample=wildcards.sample).output[0].open() as f:
        return [line.strip() for line in f.readlines()]


def get_consensus_for_passed_references_only(wildcards):
    return expand(f"results/consensus/{wildcards.sample}/{{reference}}.fa", reference=get_passed_references(wildcards))


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
        f"results/variants/{wildcards.sample}/{{reference}}/all.html",
        reference=passed_refs,
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
            ref, segment, name, tag = line.split()
            results.append(
                f"results/consensus/{wildcards.sample}/nextclade/{ref}/{segment}/{name}__{tag}/nextclade.tsv"
            )
    return results


def get_others_results(wildcards):
    with checkpoints.select_references_for_nextclade.get(sample=wildcards.sample).output.others.open() as f:
        references = [line.strip() for line in f.readlines()]
    if references:
        return expand(f"results/consensus/{wildcards.sample}/{{reference}}.fa", reference=references)
    else:
        return []


#### OUTPUTS #################################################################


def get_fastqc_reports():
    return {
        "fastqc_report": expand(
            "results/reads/{step}/fastqc/{sample}_R{orientation}.html",
            sample=SAMPLES,
            orientation=[1, 2],
            step=["trimmed"],  # TODO
        )
    }


def get_bam_outputs():
    return {
        "bams": expand(
            "results/mapping/{sample}/deduplicated/{reference}.bam",
            sample=SAMPLES,
            reference=REFERENCES,
        ),
        "nonempty_bams": expand(
            "results/checkpoints/nonempty_bams/{sample}.txt",
            sample=SAMPLES,
        ),
    }


def get_krona_reports():
    return {
        "kronas": expand("results/kraken/kronas/{sample}.html", sample=SAMPLES),
    }


def get_consensus_files():
    return {"consensus": expand("results/consensus/{sample}/nextclade/reference_summary.json", sample=SAMPLES)}


def get_mixed_positions_result():
    return {"mixed_positions": expand("results/variants/{sample}/mixed_positions_summary.txt", sample=SAMPLES)}


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
    if not "paired" in config["reads__trimming"]:
        return " ".join(args_lst)

    pe_config = config["reads__trimming"]["paired"]
    if parsed_arg := parse_paired_cutadapt_param(pe_config, "max_length_r1", "max_length_r2", "--maximum-length"):
        args_lst.append(parsed_arg)
    if parsed_arg := parse_paired_cutadapt_param(pe_config, "min_length_r1", "min_length_r2", "--minimum-length"):
        args_lst.append(parsed_arg)
    if qual_cut_arg_r1 := parse_cutadapt_comma_param(
        pe_config, "quality_cutoff_from_3_end_r1", "quality_cutoff_from_5_end_r2", "--quality-cutoff"
    ):
        args_lst.append(qual_cut_arg_r1)
    if qual_cut_arg_r2 := parse_cutadapt_comma_param(
        pe_config, "quality_cutoff_from_3_end_r1", "quality_cutoff_from_5_end_r2", "-Q"
    ):
        args_lst.append(qual_cut_arg_r2)
    return " ".join(args_lst)


def parse_samtools_params():
    samtools_params = []
    if config["consensus_params"]["count_orphans"]:
        samtools_params.append("--count-orphans")

    samtools_params.append("--max-depth {value}".format(value=config["consensus_params"]["max_read_depth"]))
    samtools_params.append("--min-MQ {value}".format(value=config["consensus_params"]["min_mapping_quality"]))
    samtools_params.append("--min-BQ {value}".format(value=config["consensus_params"]["min_base_quality"]))
    return " ".join(samtools_params)


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
