COPY (
  SELECT
    urs
  FROM xref
  where
    xref.deleted = 'N'
    -- Whole Ensembl, SGD, ENA and RegulonDB mRNAs (dbid 60-63) are not scanned, only their UTRs
    and not (
      xref.dbid IN (60, 61, 62, 63)
      and exists (
        select 1 from rnc_accessions acc
        where acc.accession = xref.ac and acc.rna_type = 'SO:0000234'
      )
    )
) TO STDOUT WITH (FORMAT CSV);
