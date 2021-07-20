BEGIN TRANSACTION;

INSERT INTO ensembl_pseudogene_regions (
  gene,
  region_name,
  chromosome,
  strand,
  region_start,
  region_stop,
  assembly_id,
  exon_count
) (
SELECT
  load.gene,
  load.region_name,
  load.chromosome,
  load.strand,
  load.exon_start,
  load.exon_stop,
  load.assembly_id,
  load.exon_count,
from load_ensembl_pseudogenes load
);

INSERT INTO ensembl_pseudogene_exons (
  pseudogene_region_id,
  exon_start,
  exon_stop
) (
SELECT
  pseudo.id,
  load.exon_start,
  load.exon_stop
FROM load_ensembl_pseudogenes
);

COMMIT;
