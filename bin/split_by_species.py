#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import csv
import re
from collections import Counter
from pathlib import Path

import click

KINGDOMS = [
    ("bacteria", "Bacteria"),
    ("archaea", "Archaea"),
    ("viruses", "Viruses"),
    ("plants", "Viridiplantae"),
    ("fungi", "Fungi"),
    ("vertebrates", "Vertebrata"),
    ("metazoa", "Metazoa"),
    ("environmental", "unclassified entries"),
]


def get_kingdom(lineage):
    for name, marker in KINGDOMS:
        if marker in lineage:
            return name
    return "protists"


def safe_name(name):
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def load_species(species_csv, out_dir):
    """Map each taxid to its kingdom directory and file name, creating the directories."""
    with open(species_csv) as handle:
        rows = [
            (taxid, get_kingdom(lineage), safe_name(name))
            for taxid, name, lineage in csv.reader(handle)
        ]

    # Disambiguate only the (kingdom, name) collisions, leaving the common case clean.
    collisions = {
        key for key, n in Counter((k, n) for _, k, n in rows).items() if n > 1
    }

    species = {}
    for taxid, kingdom, name in rows:
        (out_dir / kingdom).mkdir(parents=True, exist_ok=True)
        fname = f"{name}_{taxid}" if (kingdom, name) in collisions else name
        species[taxid] = (kingdom, fname)
    return species


def split(species_csv, fasta, out_dir):
    """
    Write each species' records to <kingdom>/<name>.fasta. The FASTA must be grouped
    by taxid (the species-specific query orders it), so one file is open at a time.
    """
    out_dir = Path(out_dir)
    species = load_species(species_csv, out_dir)
    seen = set()
    current = None
    handle = None

    with open(fasta) as fin:
        for line in fin:
            if line.startswith(">"):
                taxid = line.split()[0].rsplit("_", 1)[-1]
                if taxid != current:
                    if handle:
                        handle.close()
                        handle = None
                    if taxid in seen:
                        raise SystemExit(
                            f"taxid {taxid} appears in more than one block: the FASTA is not grouped by taxid"
                        )
                    seen.add(taxid)
                    current = taxid
                    if taxid in species:
                        kingdom, name = species[taxid]
                        handle = open(out_dir / kingdom / f"{name}.fasta", "w")
            if handle:
                handle.write(line)

    if handle:
        handle.close()


@click.command()
@click.argument("species_csv", type=click.Path(exists=True))
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("out_dir", default=".", type=click.Path())
def cli(species_csv, fasta, out_dir):
    """
    Split a taxid-grouped FASTA into one file per species, organised into kingdom
    directories.
    """
    split(species_csv, fasta, out_dir)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
