TRUNCATE table ensembl_assembly;

INSERT INTO ensembl_assembly (
    assembly_id,
    assembly_full_name,
    gca_accession,
    assembly_ucsc,
    common_name,
    taxid,
    ensembl_url,
    division,
    blat_mapping,
    example_chromosome,
    example_end,
    example_start,
    subdomain
) (
SELECT
    assembly_id,
    assembly_full_name,
    gca_accession,
    assembly_ucsc,
    common_name,
    taxid,
    ensembl_url,
    division,
    blat_mapping,
    example_chromosome,
    example_end,
    example_start,
    SUBDOMAIN
FROM load_assemblies
) ON CONFLICT (assembly_id) DO UPDATE
SET
    assembly_full_name = EXCLUDED.assembly_full_name,
    gca_accession = EXCLUDED.gca_accession,
    assembly_ucsc = EXCLUDED.assembly_ucsc,
    common_name = EXCLUDED.common_name,
    taxid = EXCLUDED.taxid,
    ensembl_url = EXCLUDED.ensembl_url,
    division = EXCLUDED.division,
    blat_mapping = EXCLUDED.blat_mapping,
    example_chromosome = EXCLUDED.example_chromosome,
    example_end = EXCLUDED.example_end,
    example_start = EXCLUDED.example_start,
    subdomain = EXCLUDED.subdomain
;
