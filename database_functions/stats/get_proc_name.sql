CREATE OR REPLACE FUNCTION stats.get_proc_name(a_number integer)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
declare
_context text;
fcesig text;
fcoid oid;
my_count integer;
begin
  get diagnostics _context = pg_context;
  fcesig := substring(_context from 'function (.*?) line');
  fcoid := to_regprocedure(fcesig::cstring);
  
  raise notice 'before query %', clock_timestamp();

  select count (*) 
    into my_count
    from rnacen.xref;

  raise notice 'after query %', clock_timestamp();

--  raise notice '% oid: % started %', fcesig, fcoid, clock_timestamp();    
  
  perform pg_sleep (2 * $1);
  
--  raise notice '% oid: % ended %', fcesig, fcoid, clock_timestamp();      
end;
$function$

