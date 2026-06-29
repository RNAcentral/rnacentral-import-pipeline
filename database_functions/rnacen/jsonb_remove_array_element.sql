CREATE OR REPLACE FUNCTION rnacen.jsonb_remove_array_element(arr jsonb, element jsonb)
 RETURNS jsonb
 LANGUAGE sql
 IMMUTABLE
AS $function$
    select arr- (
        select ordinality- 1
        from jsonb_array_elements(arr) with ordinality
        where value = element)::int
$function$

