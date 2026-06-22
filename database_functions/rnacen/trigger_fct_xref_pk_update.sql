CREATE OR REPLACE FUNCTION rnacen.trigger_fct_xref_pk_update()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
  NEW.id := nextval('xref_pk_seq');
RETURN NEW;
END;
$function$

