CREATE OR REPLACE FUNCTION rnacen.trigger_fct_rna_audit()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
BEGIN
   NEW.timestamp := LOCALTIMESTAMP;
   NEW.userstamp := USER;
RETURN NEW;
END;
$function$

