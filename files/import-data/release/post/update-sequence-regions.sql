-- This approach does have one edge case that will not be handled. That is if
-- the overall start/stop endpoints are consitent but the internal exon/intron
-- structure is different. If this happens then we will load weird stuff
-- basically.

create index ix_load_rnc_sequence_regions__accession on load_rnc_sequence_regions(accession);

-- Update the table to include urs_taxid and pretty database name
update load_rnc_sequence_regions regions
set
  urs_taxid = xref.upi || '_' || xref.taxid,
  providing_database = db.display_name
from xref, rnc_database db
where
  xref.ac = regions.accession
  and xref.deleted = 'N'
  and db.id = xref.dbid
;

-- This should have no effect but it is here for sanity.
delete from load_rnc_sequence_regions where urs_taxid is null;

-- Update name to include URS taxid
update load_rnc_sequence_regions regions
set
  region_name = urs_taxid || region_name
;

-- Change the providing_databases in rnc_sequence_features to remove all
-- mentions of the databases we have in the load table. We do this because some
-- regions may come from more than one database and we don't update all
-- databases at once so we have allow regions to remain if some other database
-- supports them.
UPDATE rnc_sequence_regions regions
set
  providing_databases = array_remove(regions.providing_databases, load.providing_database)
from load_rnc_sequence_regions load
where
  load.region_name = regions.region_name
  and load.assembly_id = regions.assembly_id
;

-- Delete all regions and exons where there are no providing databases
DELETE FROM rnc_sequence_regions WHERE was_mapped = false and providing_databases = '{}'::text[];

-- Delete all mapped locations that have the same region_name/assembly as a
-- known location. It is possible that the overall region has the same endpoints
-- but different exon/intron boundries because of mapping. So we delete the
-- mapped coordinates that will be overwritten by the given locations to be
-- loaded.
DELETE FROM rnc_sequence_regions regions
USING load_rnc_sequence_regions load
WHERE
  load.region_name = regions.region_name
  AND regions.was_mapped = true
;

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
  was_mapped,
	identity,
	providing_databases
) (
select
  max(load.urs_taxid),
  load.region_name,
  max(load.chromosome),
  max(load.strand),
  min(load.exon_start),
  max(load.exon_stop),
  load.assembly_id,
  false,
  null,
  array_agg(distinct load.providing_database)
from load_rnc_sequence_regions load
join ensembl_assembly ensembl on ensembl.assembly_id = load.assembly_id
group by load.region_name, load.assembly_id
)
ON CONFLICT (region_name, assembly_id) do UPDATE
set
  was_mapped = excluded.was_mapped,
  "identity" = excluded.identity,
  providing_databases = rnc_sequence_regions.providing_databases || excluded.providing_databases
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
)
;
