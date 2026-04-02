\timing

BEGIN TRANSACTION;

create index if not exists ix_load_rnc_sequence_regions__accession on load_rnc_sequence_regions(accession);

-- Update the table to include urs_taxid and pretty database name
update load_rnc_sequence_regions regions
set
  urs_taxid = xref.upi || '_' || xref.taxid
from xref
where
  xref.ac = regions.accession
  and xref.deleted = 'N'
;

-- Ensure we have found all URS/taxids
alter table load_rnc_sequence_regions alter column urs_taxid set not null;

-- Update name to include URS taxid
update load_rnc_sequence_regions regions
set
  region_name = urs_taxid || region_name
;

-- Map UCSC genome accesions to ensembl ones
UPDATE load_rnc_sequence_regions regions
SET
  assembly_id = genomes.assembly_id
FROM ensembl_assembly genomes
WHERE
  regions.assembly_id = genomes.assembly_ucsc
;

--- Delete all regions that come from an assembly we do not know about
DELETE from load_rnc_sequence_regions load
WHERE
  NOT EXISTS (select 1 from ensembl_assembly assem WHERE assem.assembly_id = load.assembly_id)
;

-- Workaround for mis-selection of mapped entries
update rnc_sequence_regions reg
set
        was_mapped = false
from rnc_accession_sequence_region map
where
        map.region_id = reg.id
;

-- drop indices before doing further updates/deletes
DROP INDEX IF EXISTS idx_rnc_sequence_regions_id;
DROP INDEX IF EXISTS idx_rnc_accession_sequence_region_region_id;

-- Make a copy of the data I will delete from rnc_accession_sequence_region into a backup table
-- Hopefully we can then just drop it...
drop table if exists rnc_ac_sr_backup;

-- Delete all mapped locations that are redundant with a given, but not yet
-- loaded location. These will have the same region_name/assembly as a known
-- location. It is possible that the overall region has the same endpoints but
-- different exon/intron boundaries because of mapping. So we delete the mapped
-- coordinates that will be overwritten by the given locations to be load.

-- New problems introduced by the genes infrastructure:
-- We can't just delete a region now, because of fk constraints in the rnc_gene_members
-- table. We also can't just delete genes, so we have to back up the regions
-- that were in genes, then delete them from the membership table
-- Only after this can we delete the region so it can be reloaded
-- As a precaution, we are going to mark any gene with zero members as inactive
-- at this step. They will hopefully be recovered later in the pipeline when the
-- gene prediction runs.

DROP TABLE IF EXISTS tmp_deleted_gene_members;
CREATE TABLE tmp_deleted_gene_members AS
SELECT gm.*
FROM rnc_gene_members gm
JOIN rnc_sequence_regions regions ON gm.locus_id = regions.id
JOIN load_rnc_sequence_regions load
  ON load.region_name = regions.region_name
  AND load.assembly_id = regions.assembly_id
WHERE regions.was_mapped = true;

DELETE FROM rnc_gene_members gm
USING rnc_sequence_regions regions,
      load_rnc_sequence_regions load
WHERE gm.locus_id = regions.id
  AND load.region_name = regions.region_name
  AND load.assembly_id = regions.assembly_id
  AND regions.was_mapped = true;


DELETE FROM rnc_sequence_regions regions
USING load_rnc_sequence_regions load
WHERE
  load.region_name = regions.region_name
  AND load.assembly_id = regions.assembly_id
  AND regions.was_mapped = true
;

-- We need to re-count the transcript numbers for these altered genes
UPDATE rnc_genes g
SET num_transcripts = (
  SELECT COUNT(*)
  FROM rnc_gene_members gm
  WHERE gm.rnc_gene_id = g.id
)
WHERE g.id IN (
  SELECT DISTINCT rnc_gene_id
  FROM tmp_deleted_gene_members
);

-- Any gene where there are 0 regions is now marked as inactive
UPDATE rnc_genes g
SET is_active = false
WHERE g.id IN (
  SELECT DISTINCT rnc_gene_id
  FROM tmp_deleted_gene_members
)
AND num_transcripts = 0;


-- Upsert regions table with needed info. Note the max's (for all but
-- exon_start/stop) in the select are meaningless as they will all be the same
-- for a given region_id/assembly. We also index both region_id and assembly to
-- support dealing with more than one assembly per species, ie hg19 and hg38.
insert into rnc_sequence_regions (
  urs_taxid,
  region_name,
  chromosome,
  strand,
  region_start,
  region_stop,
  assembly_id,
  exon_count,
  was_mapped,
	identity
) (
select
  max(load.urs_taxid),
  load.region_name,
  max(load.chromosome),
  max(load.strand),
  min(load.exon_start),
  max(load.exon_stop),
  load.assembly_id,
  max(load.exon_count),
  false,
  null
from load_rnc_sequence_regions load
join ensembl_assembly ensembl on ensembl.assembly_id = load.assembly_id
group by load.region_name, load.assembly_id
)
ON CONFLICT (md5(region_name), assembly_id) do UPDATE
set
  was_mapped = excluded.was_mapped,
  "identity" = excluded.identity
;

-- Populate all exons
insert into rnc_sequence_exons (
  region_id,
  exon_start,
  exon_stop
) (
select
  regions.id,
  load.exon_start,
  load.exon_stop
from load_rnc_sequence_regions load
join rnc_sequence_regions regions on regions.region_name = load.region_name
join ensembl_assembly ensembl on ensembl.assembly_id = load.assembly_id
) ON CONFLICT (region_id, exon_start, exon_stop) DO NOTHING
;

-- Populate the link from accession to sequence regions
INSERT into rnc_accession_sequence_region (
  accession,
  region_id
) (
select distinct
  load.accession,
  regions.id
from load_rnc_sequence_regions load
join rnc_sequence_regions regions on regions.region_name = load.region_name
) ON CONFLICT (accession, region_id) DO NOTHING;

do $$
declare no_exons int;
begin
	select into no_exons count(distinct t.id) from (
		select regions.id
		from rnc_sequence_regions regions
		left join rnc_sequence_exons exons on exons.region_id = regions.id
		group by regions.id
		having regions.exon_count != count(exons.*)
		) t;
	assert no_exons = 0, 'Some regions ' || no_exons || ' are missing exons';
end $$;

-- Mark all regions as having a coordinate in rnc_rna_precomputed for later
update rnc_rna_precomputed pre
set
  has_coordinates = true
from load_rnc_sequence_regions load
where
  load.urs_taxid = pre.id
;

drop table load_rnc_sequence_regions;

REFRESH MATERIALIZED VIEW active_urs_taxids;
ANALYZE active_urs_taxids;

DROP INDEX IF EXISTS ix_rnc_sequence_regions_active_provided__urs_taxid;
DROP INDEX IF EXISTS ix_rnc_sequence_regions_active_provided__assembly_id;
DROP INDEX IF EXISTS ix_rnc_sequence_regions_active_provided__urs_taxid_assembly_id;

REFRESH MATERIALIZED VIEW rnc_sequence_regions_active_provided;

CREATE INDEX ix_rnc_sequence_regions_active_provided__urs_taxid
    ON rnc_sequence_regions_active_provided(urs_taxid);
CREATE INDEX ix_rnc_sequence_regions_active_provided__assembly_id
    ON rnc_sequence_regions_active_provided(assembly_id);
CREATE INDEX ix_rnc_sequence_regions_active_provided__urs_taxid_assembly_id
    ON rnc_sequence_regions_active_provided(urs_taxid, assembly_id);

ANALYZE rnc_sequence_regions_active_provided;


DROP INDEX IF EXISTS ix_rnc_sequence_regions_active_mapped__urs_taxid;
DROP INDEX IF EXISTS ix_rnc_sequence_regions_active_mapped__assembly_id;

REFRESH MATERIALIZED VIEW rnc_sequence_regions_active_mapped;

CREATE INDEX ix_rnc_sequence_regions_active_mapped__urs_taxid
    ON rnc_sequence_regions_active_mapped(urs_taxid);
CREATE INDEX ix_rnc_sequence_regions_active_mapped__assembly_id
    ON rnc_sequence_regions_active_mapped(assembly_id);

ANALYZE rnc_sequence_regions_active_mapped;

DROP INDEX IF EXISTS ix_rnc_sequence_regions_active__urs_taxid;
DROP INDEX IF EXISTS ix_rnc_sequence_regions_active__assembly_id;

REFRESH MATERIALIZED VIEW rnc_sequence_regions_active;

CREATE INDEX ix_rnc_sequence_regions_active__urs_taxid
    ON rnc_sequence_regions_active(urs_taxid);
CREATE INDEX ix_rnc_sequence_regions_active__assembly_id
    ON rnc_sequence_regions_active(assembly_id);

ANALYZE rnc_sequence_regions_active;

CREATE INDEX idx_rnc_sequence_regions_id ON rnc_sequence_regions(id);
CREATE INDEX idx_rnc_accession_sequence_region_region_id ON rnc_accession_sequence_region(region_id);

COMMIT;
