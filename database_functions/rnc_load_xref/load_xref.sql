CREATE OR REPLACE FUNCTION rnc_load_xref.load_xref(p_previous_release bigint, p_in_dbid bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    perform rnc_load_xref.prepare_pel_tables();
    perform rnc_load_xref.populate_pel_tables1(p_previous_release);
    perform rnc_load_xref.populate_pel_tables2();
    perform rnc_load_xref.load_upi_max_versions_table(p_in_dbid);
    perform rnc_load_xref.load_max_versions_table();
    perform rnc_load_xref.populate_pel_tables3(p_in_dbid);
    perform rnc_load_xref.populate_pel_tables4(p_in_dbid, p_previous_release);
    perform rnc_load_xref.do_pel_exchange(p_in_dbid);
    -- do_checks intentionally NOT called here; run() calls it once after the
    -- per-database loop (it scans the whole xref table regardless of dbid).
END;
$function$

