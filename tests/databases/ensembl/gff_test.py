# -*- coding: utf-8 -*-

from rnacentral_pipeline.databases.ensembl import gff

GENE_ONLY = (
    "##gff-version 3\n"
    "##sequence-region Chromosome 1 1000\n"
    "# genome-build EBI ASM996v1\n"
    "Chromosome\tena\tncRNA_gene\t86216\t86293\t.\t-\t.\t"
    "ID=gene:EBG00000000222;Name=tRNA-Val;biotype=tRNA;gene_id=EBG00000000222\n"
)


def load(tmp_path, text):
    path = tmp_path / "genes.gff3"
    path.write_text(text)
    return gff.load_coordinates(path, assembly_id="ASM996v1")


# Thermococcus kodakarensis lists its ncRNAs as genes alone, with no transcript or exon.
def test_an_ncrna_gene_without_transcripts_is_its_own_transcript(tmp_path):
    info = load(tmp_path, GENE_ONLY)["EBG00000000222"]
    assert info.so_rna_type == "SO:0000253"
    (region,) = info.regions
    assert (region.chromosome, region.strand.value) == ("Chromosome", -1)
    assert [(e.start, e.stop) for e in region.exons] == [(86216, 86293)]


def test_an_unknown_biotype_falls_back_to_a_generic_ncrna(tmp_path):
    text = GENE_ONLY.replace("biotype=tRNA", "biotype=mystery")
    assert load(tmp_path, text)["EBG00000000222"].so_rna_type == "SO:0000655"
