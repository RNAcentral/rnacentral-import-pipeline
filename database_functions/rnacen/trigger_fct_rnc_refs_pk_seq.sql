CREATE OR REPLACE FUNCTION rnacen.trigger_fct_rnc_refs_pk_seq()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
  SELECT nextval('rnc_refs_pk_seq') INTO   NEW.id;
RETURN NEW;
END;
$function$

