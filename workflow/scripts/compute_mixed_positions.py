import csv
import sys
from typing import Any

sys.stderr = open(snakemake.log[0], "w")


class UnknownIvarHeaderFormat(Exception):
    """Raised when an unknown ivar header format is encountered."""


class MixedPositionDeterminator:
    alt_depth: int
    alt_freq: float
    total_depth: int

    def __init__(self, *, alt_depth: int, alt_freq: float, total_depth: int):
        self.alt_depth = alt_depth
        self.alt_freq = alt_freq
        self.total_depth = total_depth

    def _is_row_a_mixed_position(self, row: dict[str | Any, str | Any]):
        return (
            row["ALT_DP"] >= self.alt_depth and row["ALT_FREQ"] >= self.alt_freq and row["TOTAL_DP"] >= self.total_depth
        )

    def process_rows(self, rows: list[dict[str | Any, str | Any]]):
        return [row for row in rows if self._is_row_a_mixed_position(row)]


def load_ivar_variants(ivar_tsv: str):
    with open(ivar_tsv, "r") as ivar_file:
        ivar_reader = csv.DictReader(ivar_file, delimiter="\t")
        ivar_header: list[str] = ivar_reader.fieldnames
        ivar_rows = list(ivar_reader)

        if any(
            required_column not in ivar_header
            for required_column in ["REGION", "POS", "ALT_DP", "ALT_FREQ", "TOTAL_DP"]
        ):
            raise UnknownIvarHeaderFormat("Header: %s" % ivar_header)
    return ivar_rows, ivar_header


def compute_mixed_positions(
    *,
    ivar_tsv: str,
    out_mixed_positions_tsv: str,
    out_count_file: str,
    reference_name: str,
    mixed_position_determinator: MixedPositionDeterminator,
):
    ivar_rows, header = load_ivar_variants(ivar_tsv)
    mixed_positions = mixed_position_determinator.process_rows(ivar_rows)
    count = len(set([(i["REGION"], i["POS"]) for i in mixed_positions]))

    with open(out_mixed_positions_tsv, "w") as summary_file:
        summary_writer = csv.DictWriter(summary_file, delimiter="\t", fieldnames=header)
        summary_writer.writeheader()
        summary_writer.writerows(mixed_positions)

    with open(out_count_file, "w") as count_file:
        count_file.write(f"{reference_name}, {count}\n")


if __name__ == "__main__":
    mixed_position_determinator = MixedPositionDeterminator(
        alt_depth=int(snakemake.params["alt_depth"]),
        alt_freq=float(snakemake.params["alt_freq"]),
        total_depth=int(snakemake.params["total_depth"]),
    )
    compute_mixed_positions(
        ivar_tsv=snakemake.input[0],
        out_mixed_positions_tsv=snakemake.output.mixed_positions,
        out_count_file=snakemake.output.readcount,
        mixed_position_determinator=mixed_position_determinator,
        reference_name=snakemake.wildcards.reference,
    )
