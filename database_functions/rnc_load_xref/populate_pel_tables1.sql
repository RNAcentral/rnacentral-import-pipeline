CREATE OR REPLACE FUNCTION rnc_load_xref.populate_pel_tables1(v_previous_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
    -- load_retro_tmp is single-dbid, so this is the dbid for the whole table.
    v_dbid bigint := (SELECT in_dbid FROM load_retro_tmp LIMIT 1);
BEGIN

with l_data AS
(
  select x.dbid, x.created, x.urs, x.version_i, x.timestamp,
         x.userstamp, x.ac, x.version, case
           when (x.last < l.in_load_release and x.urs = l.comparable_prot_upi and (x.version = l.in_version or (x.version is null and l.in_version is null))) then l.in_load_release
           else coalesce(v_previous_release,x.last)
         end as last, case
           when (x.last < l.in_load_release and x.urs = l.comparable_prot_upi and (x.version = l.in_version or (x.version is null and l.in_version is null))) then 'N'
           else 'Y'
         end as deleted,
         case
           when (x.last < l.in_load_release and x.urs = l.comparable_prot_upi and (x.version = l.in_version or (x.version is null and l.in_version is null))) then coalesce(l.in_taxid,x.taxid)
           else x.taxid
         end taxid
  from load_retro_tmp l,
       xref x
  where x.ac = l.in_ac
  and   x.dbid = l.in_dbid
  and   x.dbid = v_dbid   -- partition pruning hint; redundant (load_retro_tmp is single-dbid)
  and   x.urs = l.comparable_prot_upi
  and   l.comparable_prot_upi is not null
  and   x.deleted = 'N'
),
ins1 AS
(
  insert into xref_pel_deleted
  (
    dbid,
    created,
    urs,
    version_i,
    timestamp,
    userstamp,
    ac,
    VERSION,
    last,
    deleted,
    taxid
  )
  select dbid, created, urs, version_i, timestamp,
         userstamp, ac, VERSION, last, deleted,
         taxid
  from l_data
  where deleted = 'Y'
)
insert into xref_pel_not_deleted
(
  dbid,
  created,
  urs,
  version_i,
  timestamp,
  userstamp,
  ac,
  version,
  last,
  deleted,
  taxid
)
select dbid, created, urs, version_i, timestamp,
       userstamp, ac, version, last, deleted,
       taxid
from l_data
where deleted = 'N';

END;
$function$
