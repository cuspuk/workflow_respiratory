rule samtools__prepare_fai_index:
    input:
        reference=f"{config['reference_panel_dirpath']}/{{reference}}.fa",
    output:
        f"{config['reference_panel_dirpath']}/{{reference}}.fa.fai",
    log:
        "logs/samtools/prepare_fai_index/{reference}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/samtools/faidx"


rule custom__create_wgs_bed_file:
    input:
        fai="resources/reference/{reference}/{reference}.fa.fai",
    output:
        bed="resources/reference/{reference}/annotation/wgs/regions.bed",
    log:
        "logs/custom/create_wgs_bed_file/{reference}.log",
    conda:
        "../envs/gawk.yaml"
    shell:
        "awk '{{print $1, 1, $2, $1}}' {input.fai} "
        " | sed 's/ /\t/g'"
        " 1> {output.bed}"
        " 2> {log};"


rule picard__prepare_dict_index:
    input:
        reference=f"{config['reference_panel_dirpath']}/{{reference}}.fa",
    output:
        seq_dict=protected(
            f"{config['reference_panel_dirpath']}/{{reference}}.dict",
        ),
    log:
        "logs/picard/prepare_dict_index/{reference}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/picard/createsequencedictionary"
