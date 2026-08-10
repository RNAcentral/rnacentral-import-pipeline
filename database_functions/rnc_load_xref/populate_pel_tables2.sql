CREATE OR REPLACE FUNCTION rnc_load_xref.populate_pel_tables2()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
    -- load_retro_tmp is single-dbid, so this is the dbid for the whole table.
    v_dbid bigint := (SELECT in_dbid FROM load_retro_tmp LIMIT 1);
BEGIN

  insert into xref_pel_not_deleted
  (
    dbid,
    created,
    upi,
    version_i,
    timestamp,
    userstamp,
    ac,
    version,
    last,
    deleted,
    taxid
  )
  select l.in_dbid, l.in_load_release created, l.comparable_prot_upi upi,
         case
           when (x.upi != l.comparable_prot_upi) then x.version_i + 1
           else x.version_i
         end,
         clock_timestamp() as "timestamp",
         user userstamp, l.in_ac ac, l.in_version as "version", l.in_load_release as "last", 'N' deleted,
         in_taxid taxid
  from load_retro_tmp l, xref x
  where x.ac = l.in_ac
  and   x.dbid = l.in_dbid
  and   x.dbid = v_dbid   -- partition pruning hint; redundant (load_retro_tmp is single-dbid)
  and   x.upi = l.comparable_prot_upi
  and   l.comparable_prot_upi is not null
  and   not -- this condition differentiates this procedure from populate_pel_tables1
        (x.last < l.in_load_release and x.upi = l.comparable_prot_upi and (x.version = l.in_version or (x.version is null and l.in_version is null)))
  and   x.deleted = 'N';

END;
$function$

