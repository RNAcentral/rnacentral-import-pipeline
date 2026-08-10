from __future__ import annotations

import argparse
import csv
import re
import typing as ty
from pathlib import Path

from rnacentral_pipeline import schemas
from rnacentral_pipeline.output_format import is_parquet
from rnacentral_pipeline.parquet_writers import typed_parquet_writer
from rnacentral_pipeline.stopfree.data import StopfreeResult

Record = tuple[str, str]
SEQUENCE_LINE = re.compile(r"^[A-Za-z]+$")


def normalize_sequence(sequence: str) -> str:
    return sequence.upper().replace("U", "T")


def read_fasta_records(path: Path) -> ty.Iterable[Record]:
    identifier: str | None = None
    chunks: list[str] = []

    with path.open("r") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield (identifier, normalize_sequence("".join(chunks)))
                identifier = line[1:].split()[0]
                chunks = []
                continue
            if line.startswith(";"):
                continue
            if not SEQUENCE_LINE.fullmatch(line):
                continue
            chunks.append(line)

    if identifier is not None:
        yield (identifier, normalize_sequence("".join(chunks)))


def build_results(
    records: ty.Iterable[Record],
    stop_free_runs: ty.Iterable[tuple[str, int]],
    gc_contents: ty.Iterable[tuple[str, float]],
    run_probabilities: ty.Iterable[tuple[str, float]],
    max_probability: float,
) -> ty.Iterable[StopfreeResult]:
    gc_by_id = dict(gc_contents)
    probability_by_id = dict(run_probabilities)

    for seq_id, run_length in stop_free_runs:
        yield StopfreeResult.build(
            seq_id,
            run_length,
            gc_by_id.get(seq_id),
            probability_by_id.get(seq_id),
            max_probability,
        )


def calculate_results(
    records: ty.Iterable[Record],
    max_probability: float,
) -> list[StopfreeResult]:
    import stopfree

    normalized_records = list(records)
    stop_free_runs = stopfree.calculate_stop_free_runs_with_ids(normalized_records)
    gc_contents = stopfree.calculate_gc_content(normalized_records)
    run_probabilities = stopfree.calculate_run_probability(stop_free_runs, gc_contents)

    return list(
        build_results(
            normalized_records,
            stop_free_runs,
            gc_contents,
            run_probabilities,
            max_probability,
        )
    )


def write_results(results: ty.Iterable[StopfreeResult], output: Path) -> None:
    if is_parquet():
        with typed_parquet_writer(
            output / "results.parquet", schemas.STOPFREE_RESULTS
        ) as writer:
            for result in results:
                writer.writerow(result.writeable())
        return

    with (output / "results.csv").open("w", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter=",",
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
        )
        for result in results:
            writer.writerow(result.writeable())


def run(input_fasta: Path, output: Path, max_probability: float) -> None:
    records = list(read_fasta_records(input_fasta))
    results = calculate_results(records, max_probability)
    output.mkdir(parents=True, exist_ok=True)
    write_results(results, output)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run stopfree on a FASTA file")
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--max-probability",
        type=float,
        default=0.0001,
        help="Maximum null-model probability still considered protein coding",
    )
    return parser


def main(argv: ty.Sequence[str] | None = None) -> int:
    args = build_argument_parser().parse_args(argv)
    run(args.input_fasta, args.output, args.max_probability)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
