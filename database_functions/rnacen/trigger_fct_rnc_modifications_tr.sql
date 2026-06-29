CREATE OR REPLACE FUNCTION rnacen.trigger_fct_rnc_modifications_tr()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
        SELECT rnc_modifications_sq.nextval
        INTO NEW.id;
    RETURN NEW;
END;
$function$

