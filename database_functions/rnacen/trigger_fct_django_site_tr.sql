CREATE OR REPLACE FUNCTION rnacen.trigger_fct_django_site_tr()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
        SELECT django_site_sq.nextval
        INTO NEW.id;
    RETURN NEW;
END;
$function$

