import locale
import os
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass


class QualimapFloatConverter:
    def convert(self, string_value: str):
        locale.setlocale(locale.LC_ALL, "en_US.UTF-8")
        return locale.atof(string_value.removesuffix("X"))


class AbstractQualimapParser(ABC):
    @abstractmethod
    def _evaluate_parsed_value(self, value: float | int) -> bool:
        """returns true if passed the criterion"""

    @abstractmethod
    def _found_correct_line(self, line: str) -> bool:
        """returns true if the line contains the correct regex"""

    @abstractmethod
    def _parse_value_from_line(self, line: str) -> float | int:
        """parses the value from the file"""

    def parse_file(self, filename: str) -> bool:
        with open(filename, "r") as file_one:
            for line in file_one:
                if self._found_correct_line(line):
                    value = self._parse_value_from_line(line)
                    return self._evaluate_parsed_value(value)
        return False


@dataclass
class MeanCoverageParser(AbstractQualimapParser):
    _REGEX_VALUE = "     mean coverageData = "
    converter = QualimapFloatConverter()
    min_coverage: float

    def _evaluate_parsed_value(self, value: float) -> bool:
        return value > self.min_coverage

    def _found_correct_line(self, line: str) -> bool:
        if line.startswith(self._REGEX_VALUE):
            return True
        return False

    def _parse_value_from_line(self, line: str) -> float | int:
        float_str = line.removeprefix(self._REGEX_VALUE).split()[0]
        return self.converter.convert(float_str)


@dataclass
class QualimapProcessor:
    filenames: list[str]
    parsers: list[AbstractQualimapParser]

    def _evaluate_file_by_parsers(self, filename: str) -> bool:
        for parser in self.parsers:
            if not parser.parse_file(filename):
                return False
        return True

    def evaluate_files(self) -> list[str]:
        genome_results_files = [os.path.join(filename, "genome_results.txt") for filename in self.filenames]
        existing_files = [filename for filename in genome_results_files if os.path.exists(filename)]
        return [filename for filename in existing_files if self._evaluate_file_by_parsers(filename)]


def build_parsers(criteria: dict[str, float | int]) -> list[AbstractQualimapParser]:
    parsers: list[AbstractQualimapParser] = []

    for criterion_name, threshold in criteria.items():
        if criterion_name == "min_mean_coverage":
            parsers.append(MeanCoverageParser(threshold))
        else:
            raise ValueError(f"Unknown criterion: {criterion_name}")

    return parsers


def evaluate_mapping_quality(qualimap_dirs: list[str], criteria: dict[str, float | int], output_file: str):
    parsers = build_parsers(criteria)
    qualimap_processor = QualimapProcessor(qualimap_dirs, parsers)
    passed_files = qualimap_processor.evaluate_files()

    print(passed_files)
    with open(output_file, "w") as out_file:
        for passed_ref in passed_files:
            out_file.write(passed_ref)
            out_file.write("\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_mapping_quality(snakemake.input, snakemake.params, snakemake.output.passed_refs)
