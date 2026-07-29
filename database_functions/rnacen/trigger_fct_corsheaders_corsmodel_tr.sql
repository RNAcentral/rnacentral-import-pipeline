CREATE OR REPLACE FUNCTION rnacen.trigger_fct_corsheaders_corsmodel_tr()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
        SELECT corsheaders_corsmodel_sq.nextval
        INTO NEW.id;
    RETURN NEW;
END;
$function$

