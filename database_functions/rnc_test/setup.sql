CREATE OR REPLACE FUNCTION rnc_test.setup()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    perform rnc_test.initialize_releases;
    perform rnc_test.initialize_rna_table;
    perform rnc_test.initialize_xref_table;

  END;

$function$

