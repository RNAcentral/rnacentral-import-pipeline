-- Copyright [2009-2014] EMBL-European Bioinformatics Institute
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--      http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

set define off

create or replace PROCEDURE REBUILD_TABLE_INDEXES (table_id IN VARCHAR2)
IS

    CURSOR index_list
    IS
        SELECT index_name
          FROM sys.all_indexes
         WHERE table_name = table_id
           AND owner in ('UNIPARC')
           AND INDEX_TYPE = 'NORMAL'
           AND PARTITIONED ='NO'
           AND TEMPORARY = 'N';

    CURSOR partitioned_index_list
    IS
        SELECT ai.index_name, aip.partition_name

          FROM sys.all_indexes ai, sys.all_ind_partitions aip
         WHERE ai.table_name = table_id
           AND ai.owner = 'UNIPARC'
           AND aip.index_owner ='UNIPARC'
           AND ai.partitioned = 'YES'
           AND TEMPORARY = 'N'
           AND SUBPARTITION_COUNT = 0
           AND ai.index_name = aip.index_name
      ORDER BY aip.partition_name asc;

    table_type VARCHAR2(3) := NULL;

    BEGIN


        SELECT PARTITIONED
          INTO table_type
          FROM sys.all_tables
         WHERE table_name = table_id
           AND OWNER = 'UNIPARC'
           AND TEMPORARY = 'N';

     IF table_type = 'NO'
     THEN
        FOR i IN index_list
        LOOP
            EXECUTE IMMEDIATE

            'ALTER INDEX ' || i.index_name || ' REBUILD TABLESPACE UNIIND';
        END LOOP;

     ELSE
        FOR p IN partitioned_index_list
        LOOP
            EXECUTE IMMEDIATE
            'ALTER INDEX ' || p.index_name || ' REBUILD PARTITION ' || p.partition_name;
        END LOOP;
     END IF;

END rebuild_table_indexes;
/
set define on
