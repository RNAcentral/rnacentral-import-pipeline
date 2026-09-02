COPY (
  select distinct xref.urs_taxid
  from xref
  join rnc_database db on db.id = xref.dbid
  where
    db.descr in (:dbs)
    and xref.deleted = 'N'
) TO STDOUT
