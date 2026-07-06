CREATE OR REPLACE FUNCTION rnacen.trigger_fct_django_content_type_tr()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
        SELECT django_content_type_sq.nextval
        INTO NEW.id;
    RETURN NEW;
END;
$function$

