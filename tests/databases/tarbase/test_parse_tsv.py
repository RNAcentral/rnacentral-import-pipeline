import pytest

from rnacentral_pipeline.databases.data.entry import Entry
from rnacentral_pipeline.databases.data.related import (
    RelatedCoordinate,
    RelatedEvidence,
    RelatedSequence,
)
from rnacentral_pipeline.databases.tarbase.parser import parse


@pytest.mark.tarbase
def test_tarbase_tsv_parse_existing():
    """
    Try to parse tsv for an existing accession and get the data we expect
    Take interation btwn hsa-miR-103a-3p and ENSG00000000460
    """

    # Accession line: "TARBASE:hsa-miR-103a-3p","","1","1","24","ncRNA","","TARBASE","hsa-miR-103a-3p","","Homo sapiens (human) hsa-miR-103a-3p","","","","","TARBASE:hsa-miR-103a-3p","","","","miRNA","{""url"": ""http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=hsa-miR-103a-3p""}","","","{}","SO:0000276","http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=hsa-miR-103a-3p"
    # Related Seq: "TARBASE:hsa-miR-103a-3p","ENSEMBL:ENSG00000000460","target_protein","{""Microarrays""}"

    related_sequences = [
        RelatedSequence(
            sequence_id=f"ENSEMBL:ENSG00000000460",
            relationship="target_protein",
            evidence=RelatedEvidence(methods=["PAR-CLIP"]),
        )
    ]
    ## NB: the method changed between version 8 & 9: Microarrays -> PAR-CLIP
    url = "https://dianalab.e-ce.uth.gr/tarbasev9/interactions?gene=C1orf112&mirna=hsa-miR-103a-3p"
    correct_entry = Entry(
        primary_id="hsa-miR-103a-3p",
        accession="TARBASE:hsa-miR-103a-3p",
        ncbi_tax_id=9606,
        database="TARBASE",
        sequence="AGCAGCATTGTACAGGGCTATGA",
        regions=[],  ## Don't create regions for protein target
        rna_type="SO:0000276",
        url=url,
        seq_version=1,
        optional_id="MIMAT0000101",
        description="Homo sapiens (human) hsa-miR-103a-3p",
        note_data={"url": url},
        xref_data={},
        related_sequences=related_sequences,
    )

    parsed_entry = parse("data/tarbase/tarbase.tsv")[0]

    assert parsed_entry == correct_entry
