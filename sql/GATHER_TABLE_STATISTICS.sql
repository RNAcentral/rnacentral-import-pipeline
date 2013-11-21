set define off

create or replace PROCEDURE        "GATHER_TABLE_STATISTICS" (table_id IN VARCHAR2)
IS

  table_name varchar2(30);

BEGIN
  table_name := upper(table_id);
  DBMS_STATS.gather_table_stats('RNACEN', table_name);

END gather_table_statistics;
/
set define on
