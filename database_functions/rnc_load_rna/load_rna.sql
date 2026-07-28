CREATE OR REPLACE FUNCTION rnc_load_rna.load_rna(p_in_dbid bigint, p_in_load_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
BEGIN

    perform rnc_load_rna.load_retro_tmp_table(p_in_dbid, p_in_load_release);
    perform rnc_load_rna.load_md5_stats_table ();
    perform rnc_load_rna.load_md5_collisions_table ();
    perform rnc_load_rna.load_md5_new_sequences_table ();
    perform rnc_load_rna.set_comparable_prot_upi ();
    perform rnc_load_rna.store_new_sequences ();

  END;

$function$

