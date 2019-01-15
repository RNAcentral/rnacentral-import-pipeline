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

-- Ensure we have found all URS/taxids
alter table load_rnc_sequence_regions alter column urs_taxid set not null;

-- Update name to include URS taxid
update load_rnc_sequence_regions regions
set
  region_name = urs_taxid || region_name
;

-- Change the providing_databases in rnc_sequence_features to remove all
-- mentions of the databases we have in the load table. We do this because some
-- regions may come from more than one database and we don't update all
-- databases at once so we have allow regions to remain if some other database
-- supports them. This has a possible edge case, if there is a case where a
-- sequence is moved, that is the overall sequence is the same but it is found
-- in a new location, this will not remove the old location. I suspect this to be
-- very rare but if it happens we will end up including an outdated location.
UPDATE rnc_sequence_regions regions
set
  providing_databases = array_remove(regions.providing_databases, load.providing_database)
from load_rnc_sequence_regions load
where
  load.region_name = regions.region_name
  and load.assembly_id = regions.assembly_id
;

-- Delete all regions and exons where there are no providing databases. Note
-- that we are relying upon cascading deletes to ensure we do not fill the exon
-- table with orphan rows.
DELETE FROM rnc_sequence_regions WHERE was_mapped = false and providing_databases = '{}'::text[];

-- Delete all mapped locations that are redundant with a given, but not yet
-- loaded location. These will have the same region_name/assembly as a known
-- location. It is possible that the overall region has the same endpoints but
-- different exon/intron boundaries because of mapping. So we delete the mapped
-- coordinates that will be overwritten by the given locations to be load.
DELETE FROM rnc_sequence_regions regions
USING load_rnc_sequence_regions load
WHERE
  load.region_name = regions.region_name
  AND load.assembly_id = regions.assembly_id
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
  exon_count,
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
  max(load.exon_count),
  false,
  null,
  array_agg(distinct load.providing_database)
from load_rnc_sequence_regions load
join ensembl_assembly ensembl on ensembl.assembly_id = load.assembly_id
group by load.region_name, load.assembly_id
)
ON CONFLICT (md5(region_name), assembly_id) do UPDATE
set
  was_mapped = excluded.was_mapped,
  "identity" = excluded.identity,
  providing_databases = rnc_sequence_regions.providing_databases || excluded.providing_databases
;

-- Ensure all providing databases are distinct
UPDATE rnc_sequence_regions
SET
  providing_databases = ARRAY(SELECT DISTINCT unnest(providing_databases))
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

drop table load_rnc_sequence_regions;
