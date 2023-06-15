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
    location_format = f"{reference_panel_dirpath}/{{name, {_REGEX}}}{_SUFFIX}"
    return set(glob_wildcards(location_format).name)


SAMPLES = glob_samples(config["sample_names_regex"])
REFERENCES = glob_references(config["reference_panel_dirpath"])


def get_constraints():
    return {
        "sample": "|".join(SAMPLES),
        "reference": "|".join(REFERENCES),
    }


#### COMMON STUFF #################################################################


def get_bwa_index(wildcards):
    base_dir = os.path.join(config["reference_panel_dirpath"], "bwa_index", wildcards.reference, wildcards.reference)
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
        non_empty_references = f.read().splitlines()
        return non_empty_references


def get_reference_fasta(wildcards):
    return os.path.join(config["reference_panel_dirpath"], f"{wildcards.reference}.fa")


def get_reference_faidx(wildcards):
    return os.path.join(config["reference_panel_dirpath"], f"{wildcards.reference}.fa.fai")


def get_passed_references(wildcards):  # TODO better way
    passed_refs = []

    with checkpoints.mapping_quality_evaluation.get(sample=wildcards.sample).output[0].open() as f:
        passed_refs = f.read().splitlines()

    ref_prefix = "results/mapping/"
    ref_suffix = f"deduplicated/bamqc/{wildcards.sample}/genome_results.txt"
    return [ref.removeprefix(ref_prefix).removesuffix(ref_suffix).rstrip("/") for ref in passed_refs]


def get_consensus_for_passed_references_only(wildcards):
    return expand(f"results/consensus/{wildcards.sample}/{{reference}}.fa", reference=get_passed_references(wildcards))


def get_mixed_positions_for_passed_references_only(wildcards):
    return expand(
        f"results/variants/{wildcards.sample}/{{reference}}/mixed_positions_count.tsv",
        reference=get_passed_references(wildcards),
    )


def get_all_qualimap_dirs(wildcards):
    return expand(
        f"results/mapping/{{reference}}/deduplicated/bamqc/{wildcards.sample}",
        reference=get_references_with_non_empty_bams(wildcards),
    )


def get_consensus_per_reference_segment(wildcards):
    segments = []
    with checkpoints.index_passed_references.get(reference=wildcards.reference).output[0].open() as f:
        for line in f.readlines():
            segments.append(line.split()[0])
    return expand(f"results/consensus/{wildcards.sample}/{wildcards.reference}/{{segment}}.fa", segment=segments)


def get_nextclade_results(wildcards):
    results = []
    with checkpoints.select_references_for_nextclade.get(sample=wildcards.sample).output[0].open() as f:
        for line in f.readlines():
            ref, name, tag = line.split()
            results.append(f"results/consensus/{wildcards.sample}/nextclade/{ref}/{name}_{tag}/nextclade.tsv")
    return results


#### OUTPUTS #################################################################


def get_fastqc_reports():
    return {
        "fastqc_report": expand(
            "results/reads/{step}/fastqc/{sample}_R{orientation}.html",
            sample=SAMPLES,
            orientation=[1, 2],
            step=["original"],  # TODO
        )
    }


def get_bam_outputs():
    return {
        "bams": expand(
            "results/mapping/{reference}/deduplicated/{sample}.bam",
            sample=SAMPLES,
            reference=REFERENCES,
        )
    }


def get_krona_reports():
    return {
        "kronas": expand("results/kraken/kronas/{sample}.html", sample=SAMPLES),
    }


def get_consensus_files():
    return {"consensus": expand("results/summary/{sample}/reference_summary.json", sample=SAMPLES)}


## PARAMETERS PARSING #################################################################


def get_illuminaclip(illumina_config):
    return ":".join(
        [
            "ILLUMINACLIP",
            illumina_config["path_adapters"],
            str(illumina_config["seed_mismatches"]),
            str(illumina_config["palindrome_clip_threshold"]),
            str(illumina_config["simple_clip_threshold"]),
        ]
    )


def get_trimmers_from_config() -> list[str]:
    trimmers = []
    if "illumina_clip" in config["reads__trimming"]:
        trimmers.append(get_illuminaclip(config["reads__trimming"]["illumina_clip"]))

    if "lead_trim" in config["reads__trimming"]:
        trimmers.append("LEADING:%s" % config["reads__trimming"]["lead_trim"])
    if "trail_trim" in config["reads__trimming"]:
        trimmers.append("TRAILING:%s" % config["reads__trimming"]["trail_trim"])
    if "head_crop" in config["reads__trimming"]:
        trimmers.append("HEADCROP:%s" % config["reads__trimming"]["head_crop"])
    if "crop_to_fixed_length" in config["reads__trimming"]:
        trimmers.append("CROP:%s" % config["reads__trimming"]["crop_to_fixed_length"])
    if "quality_threshold" in config["reads__trimming"]:
        trimmers.append("SLIDINGWINDOW:5:%s" % config["reads__trimming"]["quality_threshold"])
    if "min_length_filter" in config["reads__trimming"]:
        trimmers.append("MINLEN:%s" % config["reads__trimming"]["min_length_filter"])

    return trimmers


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
