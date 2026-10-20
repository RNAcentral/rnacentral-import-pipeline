CREATE OR REPLACE FUNCTION rnc_load_rna.load_rna(p_in_dbid bigint, p_in_load_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
BEGIN

    -- One notice per step. The gap between two lines in the log is how long the
    -- previous step took; without them the whole load is one silent block of hours.
    RAISE NOTICE 'load_rna 1/6: load_retro_tmp_table';
    perform rnc_load_rna.load_retro_tmp_table(p_in_dbid, p_in_load_release);
    RAISE NOTICE 'load_rna 2/6: load_md5_stats_table';
    perform rnc_load_rna.load_md5_stats_table ();
    RAISE NOTICE 'load_rna 3/6: load_md5_collisions_table';
    perform rnc_load_rna.load_md5_collisions_table ();
    RAISE NOTICE 'load_rna 4/6: load_md5_new_sequences_table';
    perform rnc_load_rna.load_md5_new_sequences_table ();
    RAISE NOTICE 'load_rna 5/6: set_comparable_prot_upi';
    perform rnc_load_rna.set_comparable_prot_upi ();
    RAISE NOTICE 'load_rna 6/6: store_new_sequences';
    perform rnc_load_rna.store_new_sequences ();

  END;

$function$
