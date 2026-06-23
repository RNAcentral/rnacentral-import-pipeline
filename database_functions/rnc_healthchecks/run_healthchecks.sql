CREATE OR REPLACE FUNCTION rnc_healthchecks.run_healthchecks()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    perform rnc_healthchecks.check_unique_upis();
    perform rnc_healthchecks.check_negative_taxid();
    perform rnc_healthchecks.check_unique_ac_staging();
    perform rnc_healthchecks.check_unique_ac_xref();
    perform rnc_healthchecks.check_rnas_without_xrefs();
    perform rnc_healthchecks.check_xrefs_without_lit_refs();
    perform rnc_healthchecks.check_xrefs_without_ac_data();
    perform rnc_healthchecks.check_expert_db_accessions();
    perform rnc_healthchecks.check_xref_id_not_null();

  END;

$function$

