import os
import shutil
import sys


def copy_passed_bams(bams: list[str], outdir: str):
    if not bams:
        print("There are no bams to copy", file=sys.stderr)
        return

    print("There are %s bams to copy", file=sys.stderr)
    for bam in bams:
        print(f"Copying {bam} to {outdir}", file=sys.stderr)
        shutil.copy(bam, outdir)


def cleanup_dir(dirpath: str):
    print(f"Removing directory {dirpath}", file=sys.stderr)
    shutil.rmtree(dirpath)
    os.mkdir(dirpath)


def run_copy_passed_bams(bams: list[str], outdir: str):
    cleanup_dir(outdir)
    copy_passed_bams(bams, outdir)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    run_copy_passed_bams(
        snakemake.input.bams,
        snakemake.output.output_dir,
    )
