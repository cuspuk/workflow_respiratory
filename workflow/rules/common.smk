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


def get_reference_fasta(wildcards):
    return os.path.join(config["reference_panel_dirpath"], f"{wildcards.reference}.fa")


def get_reference_faidx(wildcards):
    return os.path.join(config["reference_panel_dirpath"], f"{wildcards.reference}.fa.fai")


def get_consensus_for_passed_references_only(wildcards):
    passed_refs = []

    with checkpoints.mapping_quality_evaluation.get(sample=wildcards.sample).output[0].open() as f:
        passed_refs = f.read().splitlines()

    ref_prefix = "results/mapping/"
    ref_suffix = f"deduplicated/bamqc/{wildcards.sample}/genome_results.txt"
    reference_names = [ref.removeprefix(ref_prefix).removesuffix(ref_suffix).rstrip("/") for ref in passed_refs]
    return expand(f"results/consensus/{wildcards.sample}/{{reference}}.fa", reference=reference_names)


def get_all_qualimap_dirs(wildcards):
    return expand(f"results/mapping/{{reference}}/deduplicated/bamqc/{wildcards.sample}", reference=REFERENCES)


#### OUTPUTS #################################################################


def get_fastqc_reports():
    # steps = ["original", "trimmed", "decontaminated"] #TODO
    steps = ["original", "trimmed", "decontaminated"]

    return {
        "fastqc_report": expand(
            "results/reads/{step}/fastqc/{sample}_R{orientation}.html",
            sample=SAMPLES,
            orientation=[1, 2],
            step=steps,
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


def get_qualimap_reports():
    return {
        "qualimap_reports": expand(
            "results/mapping/{reference}/deduplicated/bamqc/{sample}", sample=SAMPLES, reference=REFERENCES
        )
    }


def get_krona_reports():
    return {
        "kronas": expand("results/kraken/kronas/{sample}.html", sample=SAMPLES),
    }


def get_consensus_files():
    return {"consensus": expand("results/summary/{sample}/aggr_result.txt", sample=SAMPLES)}


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
