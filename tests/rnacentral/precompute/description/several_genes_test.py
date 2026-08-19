# -*- coding: utf-8 -*-

"""
select_with_several_genes takes the lowest gene/locus_tag, but TAIR, SGD and
miRBase accessions can all be missing one. min() then compares None with None
and kills the whole process_range task.
"""

import typing as ty

import attr
import pytest

from rnacentral_pipeline.rnacentral.precompute.description.species_specific import (
    select_with_several_genes,
)


@attr.s()
class FakeAccession:
    """Only the attributes select_with_several_genes touches."""

    description = attr.ib()
    locus_tag = attr.ib(default=None)


def select(accessions):
    return select_with_several_genes(
        accessions,
        "genes",
        r"%s$",
        attribute="locus_tag",
        description_items="locus_tag",
        max_items=6,
    )


def test_handles_every_accession_missing_the_attribute():
    accessions = [
        FakeAccession(description="Arabidopsis thaliana tRNA"),
        FakeAccession(description="Arabidopsis thaliana tRNA-Ala"),
    ]
    assert select(accessions) == "Arabidopsis thaliana tRNA"


def test_ignores_missing_attributes_when_picking_the_candidate():
    accessions = [
        FakeAccession(description="A. thaliana ncRNA"),
        FakeAccession(description="A. thaliana AT2G01021", locus_tag="AT2G01021"),
        FakeAccession(description="A. thaliana AT1G01020", locus_tag="AT1G01020"),
    ]
    # The candidate's own locus_tag is stripped from its description, so the
    # AT1G prefix surviving here proves the null-tagged accession lost.
    assert select(accessions) == "A. thaliana (AT1G1020, AT2G1021)"


def test_appends_the_only_gene_name_it_has():
    accessions = [
        FakeAccession(description="A. thaliana ncRNA"),
        FakeAccession(description="A. thaliana ncRNA", locus_tag="AT1G01020"),
    ]
    assert select(accessions) == "A. thaliana ncRNA (AT1G01020)"
