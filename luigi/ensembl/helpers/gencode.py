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

from .. import data


def xref_data(entry):
    """
    Get GENCODE specific xref data. This will strip out the standard
    xref_data and add ENSEMBL as an xref.
    """

    xrefs = {k: v for k, v in entry.xref_data.items()}
    xrefs['Ensembl'] = [entry.accession]
    return xrefs


def primary_id(ensembl_entry):
    """
    Get a GENCODE primary id. This will use the ID prefixed with OTT and
    strip out the extra bits to create a GENCODE specific id.
    """
    return ensembl_entry.primary_id


def references(ensembl_entry):
    """
    Get the standard GENCODE reference.
    """

    return [data.Reference(
        authors=(
            "Harrow J, Frankish A, Gonzalez JM, Tapanari E, Diekhans M, "
            "Kokocinski F, Aken BL, Barrell D, Zadissa A, Searle S, Barnes"
            " I, Bignell A, Boychenko V, Hunt T, Kay M, Mukherjee G, Rajan"
            " J,  Despacio-Reyes G, Saunders G, Steward C, Harte R, Lin M,"
            " Howald  C, Tanzer A, Derrien T, Chrast J, Walters N, "
            "Balasubramanian S,  Pei B, Tress M, Rodriguez JM, Ezkurdia "
            "I, van Baren J, Brent M,  Haussler D, Kellis M, Valencia A, "
            "Reymond A, Gerstein M, Guigo  R, Hubbard TJ. "
        ),
        location="Genome Research",
        title=(
            "GENCODE: the reference human genome annotation for "
            "The ENCODE Project"
        ),
        pmid=22955987,
        doi="10.1101/gr.135350.111",
        accession=ensembl_entry.accession,
    )]
