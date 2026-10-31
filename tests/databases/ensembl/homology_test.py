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

import gzip

from rnacentral_pipeline.databases.ensembl import homology

HEADER = (
    "ref_species\tref_assembly\tquery_species\tquery_assembly\tref_gene_stable_id\t"
    "ref_gene_name\tquery_gene_stable_id\tquery_gene_name\thomology_type\t"
    "query_perc_id\tquery_perc_cov\n"
)


def mrna_line(transcript, gene, tags):
    return (
        f"1\thavana\tmRNA\t1\t9\t.\t+\t.\tID=transcript:{transcript};"
        f"Parent=gene:{gene};biotype=protein_coding;tag={tags}\n"
    )


def pair(ref_gene, query_gene):
    return f"x\tx\tx\tx\t{ref_gene}\tA\t{query_gene}\tA\thomolog_rbbh\t90.0\t100.0\n"


def test_links_canonical_mrnas_of_homologous_genes_once(tmp_path):
    human = tmp_path / "human.gff3"
    human.write_text(
        mrna_line("ENST1", "ENSG1", "Ensembl_canonical,MANE_Select")
        + mrna_line("ENST2", "ENSG1", "gencode_basic")
    )
    mouse = tmp_path / "mouse.gff3"
    mouse.write_text(mrna_line("ENSMUST1", "ENSMUSG1", "Ensembl_canonical"))

    human_homology = tmp_path / "human.tsv.gz"
    with gzip.open(human_homology, "wt") as out:
        out.write(HEADER + pair("ENSMUSG1", "ENSG1") + pair("ENSBTAG1", "ENSG1"))
    mouse_homology = tmp_path / "mouse.tsv.gz"
    with gzip.open(mouse_homology, "wt") as out:
        out.write(HEADER + pair("ENSG1", "ENSMUSG1"))

    found = list(homology.rows([human_homology, mouse_homology], [human, mouse]))

    assert sorted(t for _, t in found) == ["ENSMUST1", "ENST1"]
    assert found[0][0] == found[1][0]
