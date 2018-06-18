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


def description_of(rna_type, entries):
    """
    Compute a generic name that works for sequences that have no specific
    taxon id.

    Parameters
    ----------
    rna_type : str
        The current rna_type for the sequence.
    sequence : Rna
        An Rna model that represents the sequence.

    Returns
    -------
    name : str
        A string naming the sequence.
    """

    rna_type = rna_type.replace('_', ' ')
    species_count = len(set(e.taxid for e in entries))
    if species_count == 1:
        species = entries[0].species
        return '{species} {rna_type}'.format(
            species=species,
            rna_type=rna_type
        )

    return '{rna_type} from {species_count} species'.format(
        rna_type=rna_type,
        species_count=species_count
    )
