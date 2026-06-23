from __future__ import annotations

import csv
import typing as ty

import attr
from attr.validators import instance_of as is_a


@attr.s()
class StopfreeResult:
    urs = attr.ib(validator=is_a(str))
    stop_free_run_length = attr.ib(validator=is_a(str))
    gc_content = attr.ib(validator=is_a(str))
    run_probability = attr.ib(validator=is_a(str))
    is_protein_coding = attr.ib(validator=is_a(bool))

    @classmethod
    def build(
        cls,
        urs: str,
        stop_free_run_length: ty.Optional[int],
        gc_content: ty.Optional[float],
        run_probability: ty.Optional[float],
        max_probability: float,
    ) -> "StopfreeResult":
        is_protein_coding = (
            run_probability is not None and run_probability <= max_probability
        )
        return cls(
            urs=urs,
            stop_free_run_length=_format_value(stop_free_run_length),
            gc_content=_format_value(gc_content, "NaN"),
            run_probability=_format_value(run_probability, "NaN"),
            is_protein_coding=is_protein_coding,
        )

    def writeable(self) -> ty.List[str]:
        return [
            self.urs,
            self.stop_free_run_length,
            self.gc_content,
            self.run_probability,
            str(self.is_protein_coding),
        ]


@attr.s()
class StopfreeWriter:
    results = attr.ib(
        metadata={
            "csv_options": {
                "delimiter": ",",
                "quotechar": '"',
                "quoting": csv.QUOTE_MINIMAL,
                "lineterminator": "\n",
            }
        }
    )

    def write(self, results: ty.Iterable[StopfreeResult]):
        for result in results:
            self.results.writerow(result.writeable())


def _format_value(
    value: ty.Optional[ty.Union[int, float]],
    na_value: str = "",
) -> str:
    if value is None:
        return na_value
    return str(value)
