SELECT
    genome.assembly_id as assembly_id,
    genome.assembly_full_name as assembly_full_name,
    genome.gca_accession as gca_accession,
    genome.assembly_ucsc as assembly_ucsc,
    genome.common_name as common_name,
    genome.taxid as taxid,
    genome.ensembl_url as ensembl_url,
    genome.division as division,
    genome.blat_mapping as blat_mapping,
    genome.example_chromosome as example_chromosome,
    genome.example_start as example_start,
    genome.example_end example_end
FROM ensembl_assembly genome
;
