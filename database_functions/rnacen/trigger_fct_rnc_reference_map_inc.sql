CREATE OR REPLACE FUNCTION rnacen.trigger_fct_rnc_reference_map_inc()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
  SELECT nextval('rnc_reference_map_seq') INTO   NEW.id;
RETURN NEW;
END;
$function$

