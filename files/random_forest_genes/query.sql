COPY( select
distinct region_name,
urs_taxid,
gene,
regions.chromosome,
regions.assembly_id,
strand,
region_start,
region_stop,
exon_count,
exon_start,
exon_stop,
pc.so_rna_type,
region_stop - (region_start-1)  as length,
xref.taxid,
xref.dbid

-- *
from rnc_sequence_regions regions
left join rnc_accession_sequence_region acc_map on acc_map.region_id = regions.id
left join xref on xref.ac = acc_map.accession
left join rnc_accessions acc on acc.accession = acc_map.accession
join rnc_sequence_exons ex on
	ex.region_id = regions.id
join rnc_rna_precomputed pc on
	pc.upi = xref.upi
where
	xref.deleted = 'N'
and
/* These are all Ensembl flavours */
  xref.dbid in (25, 31, 34, 35, 36, 47)
and
	not regions.was_mapped

) TO STDOUT with CSV HEADER
