CREATE OR REPLACE FUNCTION rnc_load_rna.set_comparable_prot_upi()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    UPDATE load_retro_tmp l
       SET comparable_prot_upi = n.prot_upi
      FROM load_md5_new_sequences n
     WHERE n.in_md5 = l.in_md5
       AND l.comparable_prot_upi IS NULL;

    EXECUTE 'create index load_retro_tmp$ac_dbid_upi on load_retro_tmp(in_ac, in_dbid, comparable_prot_upi)';
    execute 'analyze load_retro_tmp';

  END;
$function$

