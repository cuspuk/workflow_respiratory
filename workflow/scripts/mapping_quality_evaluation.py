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
    name: str

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
                    print(f"{self.name} parsed {value=}", file=sys.stderr)
                    result = self._evaluate_parsed_value(value)
                    print(f"{self.name} evaluation {result=}", file=sys.stderr)
                    return result
        return False


@dataclass
class MeanCoverageParser(AbstractQualimapParser):
    name = "MeanCoverageParser"
    _REGEX_VALUE = "     mean coverageData = "
    converter = QualimapFloatConverter()
    min_coverage: float

    def __post_init__(self):
        print(f"Built MeanCoverageParser with {self.min_coverage=}", file=sys.stderr)

    def _evaluate_parsed_value(self, value: float) -> bool:
        return value > self.min_coverage

    def _found_correct_line(self, line: str) -> bool:
        if line.startswith(self._REGEX_VALUE):
            return True
        return False

    def _parse_value_from_line(self, line: str) -> float | int:
        float_str = line.removeprefix(self._REGEX_VALUE).split()[0]
        print(f"MeanCoverageParser detected {line=} and parsed {float_str}", file=sys.stderr)
        value = self.converter.convert(float_str)
        return value


@dataclass
class GenomeFraction10xParser(AbstractQualimapParser):
    name = "GenomeFraction10xParser"
    _PREFIX = "     There is a "
    _SUFFIX = "% of reference with a coverageData >= 10X"
    converter = QualimapFloatConverter()
    min_fraction: float

    def __post_init__(self):
        print(f"Built GenomeFraction10xParser with {self.min_fraction=}", file=sys.stderr)

    def _evaluate_parsed_value(self, value: float) -> bool:
        return value > self.min_fraction

    def _found_correct_line(self, line: str) -> bool:
        if line.startswith(self._PREFIX) and line.endswith(self._SUFFIX):
            return True
        return False

    def _parse_value_from_line(self, line: str) -> float | int:
        float_str = line.removeprefix(self._PREFIX).removesuffix(self._SUFFIX).strip()
        print(f"GenomeFraction10xParser detected {line=} and parsed {float_str}", file=sys.stderr)
        return self.converter.convert(float_str)


@dataclass
class QualimapProcessor:
    filenames: list[str]
    references: list[str]
    parsers: list[AbstractQualimapParser]

    def _evaluate_file_by_parsers(self, filename: str) -> bool:
        for parser in self.parsers:
            if not parser.parse_file(filename):
                return False
        return True

    def evaluate_files(self) -> list[str]:
        genome_results_files = [os.path.join(filename, "genome_results.txt") for filename in self.filenames]
        existing_files = [
            (filename, ref) for filename, ref in zip(genome_results_files, self.references) if os.path.exists(filename)
        ]
        return [ref for filename, ref in existing_files if self._evaluate_file_by_parsers(filename)]


def build_parsers(criteria: dict[str, float | int]) -> list[AbstractQualimapParser]:
    parsers: list[AbstractQualimapParser] = []

    for criterion_name, threshold in criteria.items():
        if criterion_name == "min_mean_coverage":
            parsers.append(MeanCoverageParser(threshold))
        elif criterion_name == "min_genome_fraction_with_10x_coverage":
            parsers.append(GenomeFraction10xParser(threshold))
        else:
            raise ValueError(f"Unknown criterion: {criterion_name}")

    return parsers


def evaluate_mapping_quality(
    qualimap_dirs: list[str], reference_names: list[str], criteria: dict[str, float | int], output_file: str
):
    parsers = build_parsers(criteria)
    qualimap_processor = QualimapProcessor(qualimap_dirs, reference_names, parsers)
    passed_references = qualimap_processor.evaluate_files()

    with open(output_file, "w") as out_file:
        for passed_ref in passed_references:
            out_file.write(passed_ref + "\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    evaluate_mapping_quality(
        snakemake.input.qualimaps,
        snakemake.params.reference_names,
        snakemake.params.criteria,
        snakemake.output.passed_refs,
    )
