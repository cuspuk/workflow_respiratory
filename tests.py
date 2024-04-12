import json
import os
from tempfile import TemporaryDirectory

from workflow.scripts import depths

TEST_DIR = os.path.join(".", ".tests", "test_data", "scripts")


def test_depths(filename: str, output: str):
    depths.produce_depths(filename, output)
    with open(output) as f:
        out = [depths.SegmentCoverage(**x) for x in json.load(f)]
    assert out


def run_tests(tmp_outputs_dir: str):
    test_depths(
        filename=os.path.join(TEST_DIR, "depths.txt"),
        output=os.path.join(tmp_outputs_dir, "depths.json"),
    )


if __name__ == "__main__":
    with TemporaryDirectory() as tmp_dir:
        run_tests(tmp_dir)
