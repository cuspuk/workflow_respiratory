rule samtools__prepare_fai_index:
    input:
        reference=get_reference_fasta,
    output:
        f"{config['reference_panel_dirpath']}/{{reference}}.fa.fai",
    log:
        "logs/samtools/prepare_fai_index/{reference}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/samtools/faidx"


rule picard__prepare_dict_index:
    input:
        reference=get_reference_fasta,
    output:
        seq_dict=protected(
            f"{config['reference_panel_dirpath']}/{{reference}}.dict",
        ),
    log:
        "logs/picard/prepare_dict_index/{reference}.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/v1.5.0/wrappers/picard/createsequencedictionary"
