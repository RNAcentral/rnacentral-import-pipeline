set define off

create or replace procedure add_column is

   CURSOR c_names IS
   SELECT table_name FROM user_tables
   WHERE table_name like 'LOAD_%';

begin
   FOR v_names IN c_names LOOP

       if v_names.table_name != 'LOAD_ENSEMBL_ARMADILLO' THEN

       execute immediate


       'alter table '||v_names.table_name ||' add TAXID VARCHAR2(30)';

       end if;

   END LOOP;

end add_column;
/
set define on
