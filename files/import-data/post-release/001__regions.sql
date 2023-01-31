\timing

BEGIN TRANSACTION;

create index if not exists ix_load_rnc_sequence_regions__accession on load_rnc_sequence_regions(accession);

-- Update the table to include urs_taxid
update load_rnc_sequence_regions regions
set
  urs_taxid = xref.upi || '_' || xref.taxid,
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
  "identity" = excluded.identity,
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

COMMIT;
