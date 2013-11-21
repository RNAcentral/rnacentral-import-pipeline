set define off

create or replace PROCEDURE        "GATHER_SCHEMA_STATISTICS"
IS

BEGIN

  DBMS_STATS.gather_schema_stats('RNACEN');

END gather_schema_statistics;
/
set define on
