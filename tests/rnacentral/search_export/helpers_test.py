# -*- coding: utf-8 -*-

import pytest

from rnacentral_pipeline.rnacentral.search_export import data

INSDC = "Submitted (12-JAN-2020) to the INSDC. Smith lab"
NO_ISSUE = {"has_issue": False}


def refs(items):
    return [(r["attrib"]["dbname"], r["attrib"]["dbkey"]) for r in items]


@pytest.mark.parametrize(
    "args,expected",
    [
        ((9606, ["N"], ["lncRNA"], [], NO_ISSUE), "3"),
        ((1, ["N"], ["lncRNA"], ["HGNC"], NO_ISSUE), "4"),
        ((10090, ["N"], ["lncRNA"], [], NO_ISSUE), "2"),
        ((1, ["N"], ["lncRNA"], [], NO_ISSUE), "1"),
        ((9606, ["Y"], ["lncRNA"], [], NO_ISSUE), "0"),
        ((9606, ["N"], ["misc_RNA"], [], NO_ISSUE), "2.5"),
        ((9606, ["N"], ["lncRNA"], [], {"has_issue": True}), "2.5"),
        ((9606, ["Y"], ["misc_RNA"], [], {"has_issue": True}), "-1.0"),
    ],
)
def test_boost_orders_results(args, expected):
    assert data.boost(*args) == expected


@pytest.mark.parametrize(
    "genes,products,expected",
    [
        (["TP53", None], [], {"TP53"}),
        (["foo"], ["unrelated product"], {"foo"}),
        (
            [],
            ["hsa-mir-21-5p"],
            {"hsa-mir-21-5p", "mir-21-5p", "hsa-mir-21", "mir-21"},
        ),
        (
            [],
            ["microRNA hsa-let-7a-3p"],
            {"hsa-let-7a-3p", "let-7a-3p", "hsa-let-7a", "let-7a", "let-7"},
        ),
        ([], ["hsa-mir-21a"], {"hsa-mir-21a", "mir-21a", "mir-21"}),
        ([], [None, ""], set()),
    ],
)
def test_get_genes_adds_the_generic_mirna_names(genes, products, expected):
    assert data.get_genes(genes, products) == expected


def test_as_authors_splits_and_dedupes():
    assert data.as_authors(["A, B", None, "B, C", ""]) == {"A", "B", "C"}


@pytest.mark.parametrize(
    "note,expected",
    [
        (None, {}),
        ("", {}),
        ('{"url": "x"}', {"url": "x"}),
        (
            "GO:0001 GO:0002 SO:0003 text",
            {"GO": {"GO:0001", "GO:0002"}, "SO": {"SO:0003"}},
        ),
        ("plain text", {}),
        ("{broken", {}),
    ],
)
def test_parse_note(note, expected):
    assert data.parse_note(note) == expected


def test_insdc_titles_are_separated_from_journals():
    assert data.as_journals(["Nature", INSDC, None]) == {"Nature"}
    assert list(data.as_insdc(["Nature", INSDC, None])) == ["Smith lab"]


def test_names_and_short_urs():
    assert data.as_name("URS1", 9606) == "Unique RNA Sequence URS1_9606"
    assert data.short_urs("URS0000000001", 9606) == "URS1_9606"


def test_as_ref_keeps_the_lowercase_taxonomy_name():
    assert data.as_ref("go", "x") == {"attrib": {"dbname": "GO", "dbkey": "x"}}
    assert data.as_ref("ncbi_taxonomy_id", 9606)["attrib"] == {
        "dbname": "ncbi_taxonomy_id",
        "dbkey": "9606",
    }


MIRBASE = {
    "name": "miRBase",
    "optional_id": "MI1",
    "non_coding_id": None,
    "external_id": "ext",
    "parent_accession": "A1",
    "accession": "acc",
}


def test_standard_references_for_an_expert_database():
    assert refs(data.standard_references([MIRBASE])) == [
        ("MIRBASE", "MI1"),
        ("MIRBASE", "ext"),
        ("ENA", "A1"),
    ]


def test_standard_references_add_the_accession_for_hgnc():
    hgnc = dict(
        MIRBASE, name="HGNC", optional_id=None, external_id="HGNC:5", accession="ENSG1"
    )

    assert refs(data.standard_references([hgnc])) == [
        ("HGNC", "HGNC:5"),
        ("ENA", "A1"),
        ("HGNC", "ENSG1"),
    ]


def test_references_are_deduplicated_and_end_with_the_taxonomy():
    result = data.references(9606, [MIRBASE], ["1", None, "1"], ["10.1/x"], [])

    assert refs(result) == [
        ("MIRBASE", "MI1"),
        ("MIRBASE", "ext"),
        ("ENA", "A1"),
        ("PUBMED", "1"),
        ("DOI", "10.1/x"),
        ("ncbi_taxonomy_id", "9606"),
    ]


def test_note_references_read_ontology_ids():
    notes = ['{"ontology": {"GO": ["GO:1"], "SO": ["SO:2"]}}']

    assert refs(data.note_references(notes)) == [("GO", "GO:1"), ("SO", "SO:2")]


def test_rfam_problems():
    status = {"has_issue": True, "id": 1, "messages": [], "bad_fit": True, "ok": False}

    assert data.rfam_problems(status) == ["bad_fit"]
    assert data.rfam_problems({"has_issue": False}) == ["none"]
    assert data.problem_found({"has_issue": 0}) == "False"


def test_popular_species_and_rna_type_names():
    assert data.as_popular(9606) == 9606
    assert data.as_popular(1) is None
    assert data.normalize_rna_type(["misc_RNA"]) == "misc RNA"
    assert data.as_active(["Y"]) == "Obsolete"
    assert data.as_active(["N"]) == "Active"


def test_go_annotation_fields():
    annotations = [
        {
            "qualifier": "part_of",
            "go_term_id": "GO:1",
            "go_name": "n",
            "assigned_by": "B",
        },
        {
            "qualifier": "enables",
            "go_term_id": "GO:2",
            "go_name": "m",
            "assigned_by": "A",
        },
        {
            "qualifier": "other",
            "go_term_id": "GO:3",
            "go_name": "o",
            "assigned_by": None,
        },
    ]

    assert data.part_of(annotations) == ["GO:1", "n"]
    assert data.has_go_annotations(annotations) == "True"
    assert data.has_go_annotations(annotations[2:]) == "False"
    assert data.go_source(annotations) == ["A", "B"]
    assert data.go_source([{"assigned_by": None}]) == ["Not Available"]


def test_interacting_partners():
    proteins = [
        {
            "interacting_protein_id": "UniProt:P1",
            "label": "L",
            "synonyms": ["S1"],
            "methods": ["m2"],
        },
        None,
        {
            "interacting_protein_id": None,
            "label": None,
            "synonyms": None,
            "methods": None,
        },
    ]
    rnas = [{"interacting_rna_id": "RNA:R1", "urs": "URS1", "methods": ["m1"]}]

    assert data.interacting_proteins(proteins) == ["L", "P1", "S1"]
    assert data.interacting_rnas(rnas + [None]) == {"R1", "URS1"}
    assert data.interacting_evidence([proteins[0], proteins[2]], rnas) == ["m1", "m2"]
    assert data.has_interacting_proteins([]) == "False"
    assert data.has_interacting_rnas(rnas) == "True"


def test_notes_give_diseases_urls_and_synonyms():
    notes = ['{"diseases": ["d1", {"name": "d2"}], "url": "http://x"}']

    assert data.diseases(notes) == ["d1", "d2"]
    assert data.urls(notes) == ["http://x"]
    assert data.gene_synonyms({"a;b", "c, d", None, "e"}) == {"a", "b", "c", "d", "e"}


def test_pdb_ids_keep_the_entry_and_chain():
    xrefs = [{"name": "PDBE", "accession": "1ABC_A_B"}, {"name": "X", "accession": "q"}]

    assert data.pdb_ids(xrefs) == ["1ABC_B"]


@pytest.mark.parametrize(
    "tax_string,expected",
    [
        ("Eukaryota; Metazoa; Vertebrata", "vertebrates"),
        ("Eukaryota; Metazoa; Arthropoda", "metazoa"),
        ("Bacteria; Proteobacteria", "bacteria"),
        ("Eukaryota; Alveolata", "protists"),
    ],
)
def test_kingdom(tax_string, expected):
    assert data.kingdom(tax_string) == expected


def test_kingdoms_are_sorted_and_unique():
    assert data.kingdoms(["Eukaryota; Fungi", "Bacteria", "Bacteria; x"]) == [
        "bacteria",
        "fungi",
    ]


def test_editing_event_fields():
    events = [
        {
            "chromosome": "1",
            "genomic_location": "1:5",
            "repeat_type": "Alu",
            "reference": "A",
            "edit": "G",
        }
    ]

    assert data.edit_chromosome(events) == "1"
    assert data.edit_locations(events) == ["1:5"]
    assert data.edit_repeat_type(events) == "Alu"
    assert data.edit_ref_to_edit(events) == ["A->G"]
    assert data.has_editing_event(events) == "True"
    assert data.edit_chromosome([]) == []
    assert data.edit_ref_to_edit([]) == []


def test_getters_and_value_helpers():
    assert data.get_or_empty("k")({"k": "v"}) == ["v"]
    assert data.get_or_empty("k")({"k": True}) == ["True"]
    assert data.get_or_empty("k")(None) == []
    assert data.first(["a"]) == "a"
    assert data.first([5]) == "5"
    assert data.first("abc") == "abc"
    assert data.only(5) == ["5"]
    assert data.unique(["a", "", None, "a"]) == {"a"}
    assert data.unique_lower(["A", None]) == {"a"}
    assert data.has_value([]) == ["False"]
    assert data.has_publications([1]) == "True"
