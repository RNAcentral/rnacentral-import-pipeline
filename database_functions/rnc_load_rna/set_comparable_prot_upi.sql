CREATE OR REPLACE FUNCTION rnc_load_rna.set_comparable_prot_upi()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    -- Now COMPARABLE_PROT_UPI IS NOT NULL also for new sequences

    UPDATE LOAD_RETRO_TMP L
       SET COMPARABLE_PROT_UPI = (
        SELECT N.PROT_UPI
          FROM load_md5_new_sequences N
         WHERE N.IN_MD5 = L.IN_MD5
       ) 
    WHERE COMPARABLE_PROT_UPI IS NULL;


    --COMMIT;

    EXECUTE 'create index load_retro_tmp$ac_dbid_upi on load_retro_tmp(in_ac, in_dbid, comparable_prot_upi)';
    execute 'analyze load_retro_tmp';

  END;

  /*

  * Deposit the new sequences into the main
  * table containing all sequences from all sources.
  */

$function$

