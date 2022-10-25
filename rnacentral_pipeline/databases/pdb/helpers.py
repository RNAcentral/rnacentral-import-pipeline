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

import csv
import logging
import re
import typing as ty

from rnacentral_pipeline.databases.data import AnyReference, Reference
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pubs
from rnacentral_pipeline.databases.pdb.data import ChainInfo, ReferenceMapping

RIBOSOMES = set(
    [
        "5S",
        "5.8S",
        "16S",
        "18S",
        "23S",
        "28S",
        "30S",
        "40S",
        "60S",
        "80S",
        "LSU",
    ]
)

URL = "https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}"

LOGGER = logging.getLogger(__name__)

ALLOWED = re.compile("^[ABCDFGHIKMNRSTVWXYU]+$", re.IGNORECASE)


class InvalidSequence(Exception):
    pass


class MissingProduct(Exception):
    """
    Raised when no product is available for the given chain.
    """

    pass


class MissingTypeInfo(Exception):
    """
    Raised when the chain molecule_names field is empty and we can't infer type
    """

    pass


def is_mrna(chain: ChainInfo) -> bool:
    mrna_names = [
        "mRNA",
        "Messenger RNA",
    ]
    try:
        name = product(chain)
    except MissingProduct:
        LOGGER.warn(f"No product information for {chain}, so skipping")
        return True

    if re.search("tmRNA", name, re.IGNORECASE):
        return False
    return any(re.search(n, name, re.IGNORECASE) for n in mrna_names)


def is_ncrna(info: ChainInfo) -> bool:
    if not info.molecule_type:
        return False
    return "RNA" in info.molecule_type and not is_mrna(info)


def sequence(info: ChainInfo) -> str:
    """
    Fetches the sequence of the row as DNA.
    """
    sequence = info.sequence.replace("U", "T")
    # In many tRNA's there is a single last amino acid, so we ignore that and
    # check before excluding sequences.
    if not re.match(ALLOWED, sequence):
        raise InvalidSequence("%s appears to be mislabelled protein" % info)
    return sequence


def taxid(info: ChainInfo) -> int:
    """
    Fetch the taxid from the row. This will deal with , or empty taxids.
    """

    # if no taxonomy id, use that of the synthetic construct
    if not info.taxids:
        return 32630  # synthetic construct

    # If there is >1 taxid is is synthetic
    if len(info.taxids) > 1:
        return 32630

    return info.taxids[0]


def as_reference(row):
    # This is pretty dirty but it should work assuming that each name
    # always has a ',' after both the first and last name.

    pmid = None
    if row["pubmed_id"]:
        return pubs.reference(int(row["pubmed_id"]))

    if row["doi"]:
        doi = row["doi"]
        if not doi.lower().startswith("doi:"):
            doi = "doi:" + doi
        return pubs.reference(doi)

    authors = []
    for author in row["author_list"]:
        if author["full_name"]:
            authors.append(author["full_name"].replace(",", ""))
            continue

        current = []
        if author["last_name"]:
            current.append(author["last_name"])
        if author["first_name"]:
            current.append(author["first_name"])
        authors.append(" ".join(current))

    journal = row["journal_info"]["pdb_abbreviation"]

    return Reference(
        authors=", ".join(authors),
        location=journal,
        title=row["title"],
        pmid=pmid,
        doi=row["doi"],
    )


def references_for(info: ChainInfo, mapping: ReferenceMapping) -> ty.List[AnyReference]:
    return mapping.get(reference_mapping_id(info), [])


def compound_rna_type(compound: str) -> str:
    compound = compound.upper()
    for simple_type in ["tRNA", "tmRNA", "snRNA"]:
        if simple_type.upper() in compound:
            return simple_type

    # rRNA
    for ribo_name in RIBOSOMES:
        if ribo_name in compound:
            return "rRNA"

    # SRP
    if "SRP" in compound:
        return "SRP_RNA"

    # Ribozyme
    if "RIBOZYME" in compound and "HAMMERHEAD" not in compound:
        return "ribozyme"

    # Hammerhead ribozyme
    if "RIBOZYME" in compound and "HAMMERHEAD" in compound:
        return "hammerhead_ribozyme"

    # snoRNA
    if "SNORNA" in compound:
        return "ncRNA"

    return "misc_RNA"


def rna_type(info: ChainInfo) -> str:
    if not info.molecule_names:
        raise MissingTypeInfo(
            f'Cannot find RNA type for {info}, falling back to "misc_RNA"'
        )
        return "misc_RNA"
    return compound_rna_type(info.molecule_names[0])


def url(info: ChainInfo) -> str:
    """
    Generate a URL for a given result. It will point to the page for the whole
    structure.
    """
    return URL.format(pdb_id=info.pdb_id.lower())


def note_data(info: ChainInfo) -> ty.Dict[str, str]:
    data = {
        "releaseDate": info.release_day(),
        "structureTitle": info.title,
    }

    if info.resolution is not None:
        data["resolution"] = str(info.resolution)

    if info.experimental_method:
        data["experimentalTechnique"] = info.experimental_method.upper()

    return data


def description(info: ChainInfo, max_length=80) -> str:
    compound = product(info)[:max_length] + (product(info)[max_length:] and "...")
    return "{compound} from {source} (PDB {pdb}, chain {chain})".format(
        compound=compound,
        source=info.organism_scientific_name,
        pdb=info.pdb_id.upper(),
        chain=info.chain_id,
    )


def product(info: ChainInfo) -> str:
    if not info.molecule_names:
        raise MissingProduct(f"No products in: {info}")
    return info.molecule_names[0]


def reference_mapping_id(info: ChainInfo) -> str:
    return info.pdb_id


def lineage(info: ChainInfo) -> str:
    return phy.lineage(taxid(info))


def species(info: ChainInfo) -> str:
    return phy.species(taxid(info))


def load_overrides(handle) -> ty.Set[ty.Tuple[str, str]]:
    """
    Parse TSV file of pdb_id chain and produce a set of (pdb_id, chain). PDB id
    will be lowercased. This is used to ensure all sequences with an Rfam match
    are loaded into the pipeline.
    """
    overrides = set()
    for row in csv.reader(handle, delimiter="\t"):
        overrides.add((row[0].lower(), row[1]))
    return overrides
