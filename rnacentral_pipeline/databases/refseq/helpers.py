# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import re

import attr

import rnacentral_pipeline.databases.helpers.embl as embl
from rnacentral_pipeline.databases import data as dat

URL = "https://www.ncbi.nlm.nih.gov/nuccore/{primary_id}.{version}"

NCRNA = {
    "ncRNA",
    "precursor_RNA",
    "rRNA",
}


def optional_id(feature):
    for ref in feature.qualifiers["db_xref"]:
        if ref.startswith("GeneID:"):
            return ref
    return None


def parent_accession(record):
    return record.annotations["accessions"][0]


def url(record):
    return URL.format(
        primary_id=primary_id(record),
        version=embl.seq_version(record),
    )


def description(record, feature):
    product = embl.product(feature)

    gene_feature = record.features[1]
    if product == "other RNA":
        product = gene_feature.qualifiers.get("note", [])
        if len(product) == 1:
            product = product[0]
            if ";" in product:
                product = product[product.index(";") + 1 :].strip()
            # Only lower case things that do not look like genes
            if product[0].isupper() and product[1].islower():
                product = product[0].lower() + product[1:]

    gene = embl.gene(feature) or embl.locus_tag(feature)
    if gene:
        gene = " (%s)" % gene
    else:
        gene = ""

    return "{organism} {product}{gene}".format(
        organism=embl.organism(record.features[0])[0],
        product=product,
        gene=gene,
    )


def primary_id(record):
    return record.annotations["accessions"][0]


def accession(record, feature):
    return "{primary_id}.{version}:{start}..{stop}:{feature_type}".format(
        primary_id=primary_id(record),
        version=embl.seq_version(record),
        start=feature.location.start + 1,
        stop=feature.location.end,
        feature_type=feature.type,
    )


def xref_data(feature):
    feature_annotations = feature.qualifiers.get("db_xref", [])
    return embl.grouped_annotations(feature_annotations, ":")


def rna_type(record, feature):
    rna_type = embl.rna_type(feature)
    if rna_type != "precursor_RNA":
        return rna_type

    products = feature.qualifiers.get("product", [])
    if products:
        prod = products[0].lower()
        if prod.startswith("microrna") or prod.startswith("pre-microrna"):
            return "SO:0001244"

    genes = feature.qualifiers.get("gene", [])
    for gene in genes:
        if gene.lower().startswith("mir"):
            return "SO:0001244"

    syn = feature.qualifiers.get("gene_synonym", [])
    if syn and any("mir" in s.lower() for s in syn):
        return "SO:0001244"
    return rna_type


def as_entry(record, source, feature):
    """
    Create an Entry based upon the record, source feature and ncRNA feature.
    """

    try:
        optional = optional_id(feature)
    except:
        optional = optional_id(source)

    return dat.Entry(
        primary_id=primary_id(record),
        accession=accession(record, feature),
        ncbi_tax_id=embl.taxid(record),
        database="REFSEQ",
        sequence=embl.sequence(record, feature),
        regions=[],
        rna_type=rna_type(record, feature),
        url=url(record),
        seq_version=embl.seq_version(record),
        optional_id=optional,
        note_data={},
        xref_data=xref_data(feature),
        species=embl.species(record),
        common_name=embl.common_name(record),
        lineage=embl.lineage(record),
        gene=embl.gene(feature),
        locus_tag=embl.locus_tag(feature),
        product=embl.product(feature),
        parent_accession=parent_accession(record),
        project=embl.project(record),
        organelle=embl.organelle(source),
        inference=embl.inference(feature),
        standard_name=embl.standard_name(feature),
        description=description(record, feature),
        mol_type=embl.mol_type(source),
        is_composite="N",
        gene_synonyms=embl.gene_synonyms(feature),
        references=embl.references(record),
    )


def ncrna_features(features):
    return [f for f in features if f.type in NCRNA]


def generate_related(entries):
    """
    This goes through all given entries, which are assumed to all be from the
    same gene, and thus splicing variants, and populates the related_sequences
    feature with the required related sequence information.
    """

    for first in entries:
        related = first.related_sequences
        for second in entries:
            # We have to test the object and accession in the case that there
            # are duplicate LOCUS entries and one of them is already processed
            # and has a related mature product.
            if first == second or first.accession == second.accession:
                continue

            relationship = "isoform"
            if first.rna_type == "SO:0001244":
                if second.rna_type == "SO:0000276":
                    relationship = "mature_product"
                else:
                    msg = "Unknown type of relationship between %s\t%s"
                    raise ValueError(msg % (first, second))

            elif first.rna_type == "SO:0000276":
                if second.rna_type == "SO:0001244":
                    relationship = "precursor"
                elif second.rna_type == "SO:0000276":
                    continue
                else:
                    print(first.rna_type)
                    print(second.rna_type)
                    msg = "Unknown type of relationship between %s\t%s"
                    raise ValueError(msg % (first, second))

            related.append(
                dat.RelatedSequence(
                    sequence_id=second.accession,
                    relationship=relationship,
                    coordinates=[],
                    evidence=dat.RelatedEvidence.empty(),
                )
            )
        yield attr.evolve(first, related_sequences=related)
