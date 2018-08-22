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

import json

from lxml import etree as ET

from rnacentral_pipeline.rnacentral.search_export import data


def test_create_tag_with_complete_dict():
    root = ET.Element('a')
    entry = {'attrib': {'name': 'bob'}, 'text': 'nice'}
    element = data.create_tag(root, 'example', entry)
    assert element.tag == 'example'
    assert element.attrib == {'name': 'bob'}
    assert element.text == 'nice'
    assert element in list(root)


def create_tag_can_handle_unicode():
    root = ET.Element('a')
    entry = u'tRNA Asp ⊄UC'
    element = data.create_tag(root, 'example', entry)
    assert element.tag == 'example'
    assert element.attrib == {}
    assert element.text == u'tRNA Asp ⊄UC'
    assert element in list(root)



def test_create_tag_can_use_simple_text():
    root = ET.Element('a')
    element = data.create_tag(root, 'example', 'hi')
    assert element.tag == 'example'
    assert element.attrib == {}
    assert element.text == 'hi'
    assert element in list(root)


def test_create_tag_can_use_simple_text_and_attribs():
    root = ET.Element('a')
    attrib = {'type': 'important'}
    element = data.create_tag(root, 'example', 'hi', attrib=attrib)
    assert element.tag == 'example'
    assert element.attrib == {'type': 'important'}
    assert element.text == 'hi'
    assert element in list(root)


def test_create_tag_overrides_given_attributes_with_tag_dict():
    root = ET.Element('a')
    entry = {'attrib': {'name': 'bob'}, 'text': 'nice'}
    attrib = {'name': 'steve', 'help': 'never'}
    element = data.create_tag(root, 'example', entry, attrib=attrib)
    assert element.attrib == {'name': 'bob', 'help': 'never'}
    assert element.text == 'nice'
    assert element in list(root)


def test_parse_note_can_parse_space_note():
    assert data.parse_note('GO:3333 SO:5000 extra') == {
        'GO': {'GO:3333'},
        'SO': {'SO:5000'},
    }


def test_parse_note_can_handle_json_note():
    ans = {'ontology': {'GO': ['GO:100']}}
    assert data.parse_note(json.dumps(ans)) == ans


def test_can_build_references_with_complex_note():
    refs = data.note_references([
        '6911',
        '6911',
        'Quality score:93.74 Alignment:reference ECO:0000053 SO:0001000 GO:0003735 GO:0005840',
    ])
    assert sorted(refs) == [
        {'attrib': {'dbname': 'ECO', 'dbkey': 'ECO:0000053'}},
        {'attrib': {'dbname': 'GO', 'dbkey': 'GO:0003735'}},
        {'attrib': {'dbname': 'GO', 'dbkey': 'GO:0005840'}},
        {'attrib': {'dbname': 'SO', 'dbkey': 'SO:0001000'}},
    ]

def test_references_build_correct_data():
    xrefs = []
    notes = [
        '',  # Old style empty note
        'GO:3333 SO:5000 extra',  # Old style note
        None,  # New style empty note
        json.dumps({'ontology': {'GO': ['GO:100']}}),
        json.dumps({'GO': ['GO:500', 'GO:100']}),
    ]
    val = data.references(
        10,
        xrefs,
        [10000, None, 50000],
        [None, '10.1093/nar/9.13.2999', '', '10.1016/s0923-1811(03)00030-6'],
        notes
    )
    assert val == [
        {'attrib': {'dbname': 'PUBMED', 'dbkey': '10000'}},
        {'attrib': {'dbname': 'PUBMED', 'dbkey': '50000'}},
        {'attrib': {'dbname': 'DOI', 'dbkey': '10.1093/nar/9.13.2999'}},
        {'attrib': {'dbname': 'DOI', 'dbkey': '10.1016/s0923-1811(03)00030-6'}},
        {'attrib': {'dbname': 'GO', 'dbkey': 'GO:3333'}},
        {'attrib': {'dbname': 'SO', 'dbkey': 'SO:5000'}},
        {'attrib': {'dbname': 'GO', 'dbkey': 'GO:100'}},
        {'attrib': {'dbname': 'GO', 'dbkey': 'GO:500'}},
        {'attrib': {'dbname': 'ncbi_taxonomy_id', 'dbkey': '10'}},
    ]

def notes_can_parse_weird_notes():
    notes = [
        'Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence NAS] [PMID 9159473]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence TAS] [PMID 9443958]; GO_process: GO:0006414 - translational elongation [Evidence NAS] [PMID 9159473]',
        '{"text": ["Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence IC] [PMID 9023104]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence ISM] [PMID 9023104]; GO_process: GO:0006414 - translational elongation [Evidence IC] [PMID 9023104]"]}',
        '{"isotype": "Val", "anticodon_positions": [{"relative_stop": 37, "relative_start": 35}], "anticodon": "AAC", "url": "http://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/genes/tRNA-Val-AAC-2-1.html", "score": 66.3}',
        'Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence NAS] [PMID 9159473]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence TAS] [PMID 9443958]; GO_process: GO:0006414 - translational elongation [Evidence NAS] [PMID 9159473]',
        '{"text": ["Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence IC] [PMID 9023104]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence ISM] [PMID 9023104]; GO_process: GO:0006414 - translational elongation [Evidence IC] [PMID 9023104]"]}',
        '{"text": ["Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence IC] [PMID 9023104]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence ISM] [PMID 9023104]; GO_process: GO:0006414 - translational elongation [Evidence IC] [PMID 9023104]"]}',
        '{"text": ["Covariance Model: Eukaryote; CM Score: 67.95", "Legacy ID: chrXV.trna15-ValAAC"], "ontology": ["ECO:0000202", "GO:0030533", "SO:0000253"]}',
        'Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence NAS] [PMID 9159473]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence TAS] [PMID 9443958]; GO_process: GO:0006414 - translational elongation [Evidence NAS] [PMID 9159473]',
        '{"text": ["Covariance Model: Eukaryote; CM Score: 67.95", "Legacy ID: chrXV.trna15-ValAAC"], "ontology": ["ECO:0000202", "GO:0030533", "SO:0000253"]}',
        'Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence NAS] [PMID 9159473]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence TAS] [PMID 9443958]; GO_process: GO:0006414 - translational elongation [Evidence NAS] [PMID 9159473]',
        '{"text": ["Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence IC] [PMID 9023104]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence ISM] [PMID 9023104]; GO_process: GO:0006414 - translational elongation [Evidence IC] [PMID 9023104]"]}',
        'Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence NAS] [PMID 9159473]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence TAS] [PMID 9443958]; GO_process: GO:0006414 - translational elongation [Evidence NAS] [PMID 9159473]',
        '{"text": ["Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence IC] [PMID 9023104]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence ISM] [PMID 9023104]; GO_process: GO:0006414 - translational elongation [Evidence IC] [PMID 9023104]"]}',
        '{"text": ["Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence IC] [PMID 9023104]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence ISM] [PMID 9023104]; GO_process: GO:0006414 - translational elongation [Evidence IC] [PMID 9023104]"]}',
        '{"text": ["Valine tRNA (tRNA-Val), predicted by tRNAscan-SE analysis; GO_component: GO:0005829 - cytosol [Evidence IC] [PMID 9023104]; GO_function: GO:0030533 - triplet codon-amino acid adaptor activity [Evidence ISM] [PMID 9023104]; GO_process: GO:0006414 - translational elongation [Evidence IC] [PMID 9023104]"]}',
    ]
    refs = data.note_references(notes)
    assert refs
