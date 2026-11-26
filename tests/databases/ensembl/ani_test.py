# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.ensembl import ani


def report(accession, status=None, own=None, best=None):
    found = {"accession": accession}
    if status:
        check = {"taxonomy_check_status": status}
        if own is not None:
            check["submitted_ani_match"] = {"ani": own}
        if best is not None:
            check["best_ani_match"] = {"ani": best}
        found["average_nucleotide_identity"] = check
    return found


def test_only_failed_genomes_are_listed():
    reports = [
        report("GCA_1.1", "OK"),
        report("GCA_2.1", "Inconclusive"),
        report("GCA_3.1", "Failed"),
        report("GCA_4.1"),
    ]
    assert list(ani.failed(reports)) == ["GCA_3.1"]


def test_an_inconclusive_genome_closer_to_another_species_is_listed():
    reports = [
        report("GCA_1.1", "Inconclusive", own=82.76, best=99.99),
        report("GCA_2.1", "Inconclusive", own=96.1, best=99.2),
        report("GCA_3.1", "Inconclusive", own=82.0, best=91.0),
        report("GCA_4.1", "Inconclusive", best=99.0),
    ]
    assert list(ani.failed(reports)) == ["GCA_1.1"]


def test_every_assembly_in_species_json_is_checked():
    data = {
        "species": {
            "A": {"assemblies": {"GCA_2.1": {}, "GCA_1.1": {}}},
            "B": {"assemblies": {"GCA_1.1": {}}},
            "C": {},
        }
    }
    assert ani.assemblies(data) == ["GCA_1.1", "GCA_2.1"]
