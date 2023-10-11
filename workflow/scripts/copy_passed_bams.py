import os
import shutil
import sys


def copy_passed_bams(bams: list[str], outdir: str):
    if not bams:
        print("There are no bams to copy", file=sys.stderr)
        return

    print(f"There are {len(bams)} bams to copy", file=sys.stderr)
    for bam in bams:
        print(f"Copying {bam} to {outdir}", file=sys.stderr)
        shutil.copy(bam, outdir)


def cleanup_dir(dirpath: str):
    print(f"Removing directory {dirpath}", file=sys.stderr)
    shutil.rmtree(dirpath)


def run_copy_passed_bams(bams: list[str], outdir: str):
    if os.path.exists(outdir):
        cleanup_dir(outdir)
    os.mkdir(outdir)
    copy_passed_bams(bams, outdir)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    run_copy_passed_bams(
        snakemake.input.bams_and_idxes,
        snakemake.output.output_dir,
    )
