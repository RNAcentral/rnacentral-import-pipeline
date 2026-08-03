# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import collections as col
import csv
import itertools as it
import json
import operator as op
import os
import typing as ty
from contextlib import ExitStack
from pathlib import Path

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional
from sqlitedict import SqliteDict

from rnacentral_pipeline import schemas
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.output_format import is_parquet
from rnacentral_pipeline.parquet_writers import typed_parquet_writer

NAME_ALIASES = {
    "common name",
    "equivalent name",
    "genbank common name",
    "genbank synonym",
    "scientific name",
    "synonym",
}


@attr.s(hash=True)
class TaxonomyEntry(object):
    tax_id = attr.ib(validator=is_a(int))
    name = attr.ib(validator=is_a(str))
    lineage = attr.ib(validator=is_a(str))
    aliases = attr.ib(validator=is_a(list), hash=False)
    replaced_by = attr.ib(validator=optional(is_a(int)))
    rank = attr.ib(validator=is_a(str), default="")
    reference_proteome = attr.ib(validator=is_a(bool), default=False)
    is_deleted = attr.ib(validator=is_a(bool), default=False)

    @classmethod
    def build(cls, entry, names, rank="", reference_proteome=False):
        aliases = set()
        for name_entry in names:
            (tax_id, name, _, name_class) = name_entry
            assert tax_id == entry[0]
            if name_class in NAME_ALIASES:
                aliases.add(name)

        return cls(
            tax_id=int(entry[0]),
            name=entry[1],
            lineage=entry[2] + entry[1],
            aliases=sorted(aliases),
            replaced_by=None,
            rank=rank,
            reference_proteome=reference_proteome,
        )

    @classmethod
    def build_from_ena(cls, tax_id):
        """Build an entry for a taxid by resolving it against the ENA taxonomy
        service (with the UniProt fallback baked into the phylogeny helper).
        This is how we recover taxids that RNAcentral references but that are
        absent from the NCBI taxdump file, because ENA assigns/exposes taxids
        (especially environmental/uncultured ones) ahead of the taxdump.

        Raises phy.UnknownTaxonId / phy.FailedTaxonId if the taxid cannot be
        resolved, so callers can fail loudly rather than silently drop it."""
        data = phy.phylogeny(int(tax_id))
        name = data["scientificName"]
        lineage = data.get("lineage") or ""
        return cls(
            tax_id=int(tax_id),
            name=name,
            lineage="{}{}".format(lineage, name),
            aliases=[],
            replaced_by=None,
            rank=data.get("rank") or "",
            reference_proteome=False,
            is_deleted=False,
        )

    @classmethod
    def build_deleted(cls, tax_id):
        """Build an entry for a taxid that NCBI has deleted (listed in
        delnodes.dmp). Deleted taxids are stripped from every other dump, so we
        have no name/lineage for them; we store a sentinel name and flag the row
        as deleted. The loader preserves any name/lineage we already had for a
        taxon that was live in a previous release (see 000__taxonomy.sql)."""
        return cls(
            tax_id=int(tax_id),
            name="deleted taxid {}".format(tax_id),
            lineage="",
            aliases=[],
            replaced_by=None,
            rank="",
            reference_proteome=False,
            is_deleted=True,
        )

    def writeable(self):
        yield [
            self.tax_id,
            self.name,
            self.lineage,
            json.dumps(self.aliases),
            self.replaced_by,
            self.rank,
            self.reference_proteome,
            self.is_deleted,
        ]


def ncbi_reader(handle):
    def cleaned_lines(to_clean):
        for line in to_clean:
            cleaned = line.replace("\t|\n", "\n").replace("\t|\t", "\t")
            yield cleaned

    return csv.reader(cleaned_lines(handle), delimiter="\t")


def grouped_extra(handle, group_idx=0):
    reader = ncbi_reader(handle)
    data = col.defaultdict(list)
    for key, values in it.groupby(reader, op.itemgetter(group_idx)):
        data[key].extend(list(values))
    return data


def parse_nodes(handle):
    reader = ncbi_reader(handle)
    return {row[0]: row[2] for row in reader}


def parse_delnodes(handle):
    reader = ncbi_reader(handle)
    return {row[0] for row in reader if row}


def parse_ref_proteomes(handle):
    reader = csv.reader(handle, delimiter="\t")
    for _ in range(15):
        next(reader)  # skip preamble & header
    return {int(row[1]) for row in reader}


def parse(
    handle,
    names_handle,
    merged_handle,
    nodes_handle,
    delnodes_handle,
    ref_proteomes_handle=None,
):
    lineage = ncbi_reader(handle)
    names = grouped_extra(names_handle)
    merged = grouped_extra(merged_handle, group_idx=1)
    nodes = parse_nodes(nodes_handle)
    deleted = parse_delnodes(delnodes_handle)
    ref_proteomes = (
        parse_ref_proteomes(ref_proteomes_handle) if ref_proteomes_handle else set()
    )

    for raw in lineage:
        possible_names = names.get(raw[0], [])
        rank = nodes.get(raw[0], "")
        is_ref_proteome = int(raw[0]) in ref_proteomes
        entry = TaxonomyEntry.build(
            raw, possible_names, rank=rank, reference_proteome=is_ref_proteome
        )
        yield entry

        for (old_tax_id, replaced) in merged.get(raw[0], []):
            assert int(replaced) == entry.tax_id
            yield attr.evolve(entry, tax_id=int(old_tax_id), replaced_by=entry.tax_id)

    # Taxids NCBI has deleted outright. These are stripped from the other dumps,
    # so they never appear in the lineage loop above, but RNAcentral may still
    # hold active accessions stamped with them. Skip any that are somehow still
    # live nodes as a guard against an inconsistent dump.
    for del_tax_id in deleted:
        if del_tax_id in nodes:
            continue
        yield TaxonomyEntry.build_deleted(del_tax_id)


def parse_directory(
    directory: Path, ref_proteomes_path=None
) -> ty.Iterable[TaxonomyEntry]:
    names = [
        "fullnamelineage.dmp",
        "names.dmp",
        "merged.dmp",
        "nodes.dmp",
        "delnodes.dmp",
    ]
    filenames = [os.path.join(directory, name) for name in names]
    with ExitStack() as stack:
        files = [stack.enter_context(open(f)) for f in filenames]
        if ref_proteomes_path:
            ref_handle = stack.enter_context(open(ref_proteomes_path))
            yield from parse(*files, ref_proteomes_handle=ref_handle)
        else:
            yield from parse(*files)


def write(directory: Path, output, ref_proteomes_path=None):
    if isinstance(output, (str, Path)):
        path = Path(output)
        if is_parquet():
            with typed_parquet_writer(path, schemas.TAXONOMY) as writer:
                for entry in parse_directory(
                    directory, ref_proteomes_path=ref_proteomes_path
                ):
                    for row in entry.writeable():
                        writer.writerow(row)
            return
        with path.open("w") as handle:
            writer = csv.writer(handle)
            for entry in parse_directory(
                directory, ref_proteomes_path=ref_proteomes_path
            ):
                writer.writerows(entry.writeable())
        return

    writer = csv.writer(output)
    for entry in parse_directory(directory, ref_proteomes_path=ref_proteomes_path):
        writer.writerows(entry.writeable())


class UnresolvableTaxids(Exception):
    """Raised when one or more referenced taxids cannot be resolved against
    ENA/UniProt. Carries every failure so the operator sees the full list in
    one go rather than one taxid at a time."""

    def __init__(self, failures):
        self.failures = failures
        detail = ", ".join("{} ({})".format(tid, reason) for tid, reason in failures)
        super().__init__(
            "Could not resolve {} taxid(s) against ENA/UniProt: {}".format(
                len(failures), detail
            )
        )


def resolve_missing(taxids) -> ty.List[TaxonomyEntry]:
    """Resolve an iterable of taxids against ENA. Resolves every taxid before
    returning anything: if any fail, raises UnresolvableTaxids listing them all
    and yields no partial result, so a caller never loads a half-filled set."""
    entries = []
    failures = []
    for taxid in taxids:
        try:
            entries.append(TaxonomyEntry.build_from_ena(taxid))
        except (phy.UnknownTaxonId, phy.FailedTaxonId) as err:
            failures.append((str(taxid), type(err).__name__))
    if failures:
        raise UnresolvableTaxids(failures)
    return entries


def write_missing(taxid_handle, output):
    """Read taxids (one per line) that are missing from rnc_taxonomy, resolve
    them against ENA, and write rows in the rnc_taxonomy load format. Fails
    loudly (UnresolvableTaxids) if any taxid cannot be resolved."""
    taxids = [line.strip() for line in taxid_handle if line.strip()]
    entries = resolve_missing(taxids)
    writer = csv.writer(output)
    for entry in entries:
        writer.writerows(entry.writeable())


def index(directory: Path, output: str):
    mapping = SqliteDict(filename=output)
    for entry in parse_directory(directory):
        mapping[str(entry.tax_id)] = entry
    mapping.commit()
