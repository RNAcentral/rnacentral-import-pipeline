CREATE OR REPLACE FUNCTION stats.call_a_func(a_number integer)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
declare
_context text;
fcesig text;
fcoid oid;
begin
  get diagnostics _context = pg_context;
  fcesig := substring(_context from 'function (.*?) line');
  fcoid := to_regprocedure(fcesig::cstring);
  raise notice '% oid: % started %', fcesig, fcoid, clock_timestamp();    

  perform stats.get_proc_name (a_number);
  
  perform pg_sleep (3 * $1);

  perform stats.get_proc_name (a_number);

  raise notice '% oid: % ended %', fcesig, fcoid, clock_timestamp();      
end;
$function$

