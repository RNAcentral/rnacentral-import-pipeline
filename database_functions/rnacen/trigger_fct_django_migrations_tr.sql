CREATE OR REPLACE FUNCTION rnacen.trigger_fct_django_migrations_tr()
 RETURNS trigger
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
        SELECT django_migrations_sq.nextval
        INTO NEW.id;
    RETURN NEW;
END;
$function$

