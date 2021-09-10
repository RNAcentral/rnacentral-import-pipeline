\timing

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
SELECT distinct
  load.gene,
  load.region_name,
  load.chromosome,
  load.strand,
  load.exon_start,
  load.exon_stop,
  load.assembly_id,
  load.exon_count
from load_ensembl_pseudogenes load
join ensembl_assembly assem
on
  assem.assembly_id = load.assembly_id
) ON CONFLICT (md5(region_name)) DO NOTHING;

INSERT INTO ensembl_pseudogene_exons (
  region_id,
  exon_start,
  exon_stop
) (
SELECT
  pseudo.id,
  load.exon_start,
  load.exon_stop
FROM load_ensembl_pseudogenes load
join ensembl_pseudogene_regions pseudo
on
   pseudo.region_name = load.region_name
) ON CONFLICT (region_id, exon_start, exon_stop) DO NOTHING ;

COMMIT;
