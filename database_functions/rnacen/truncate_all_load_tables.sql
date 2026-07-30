CREATE OR REPLACE FUNCTION rnacen.truncate_all_load_tables()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE


 c_get_load_tables CURSOR FOR
  select distinct tablename
  from pg_tables
  where schemaname = 'rnacen'
  and   tablename like 'load_%';

BEGIN
    for load_table in c_get_load_tables loop
      RAISE NOTICE 'truncate table % ', load_table.tablename;
      EXECUTE 'truncate table '||load_table.tablename;
    end loop;
  end;
$function$

