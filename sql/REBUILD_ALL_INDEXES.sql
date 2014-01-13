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

create or replace PROCEDURE REBUILD_ALL_INDEXES
IS
    CURSOR index_list
    IS
        SELECT index_name
          FROM sys.all_indexes
         WHERE owner in ('UNIPARC')
           AND INDEX_TYPE = 'NORMAL'
           AND PARTITIONED ='NO'
           AND TEMPORARY = 'N'
      ORDER BY table_name desc;

      CURSOR partitioned_index_list

      IS
        SELECT ai.index_name, ai.table_name, aip.partition_name
          FROM sys.all_indexes ai, sys.all_ind_partitions aip
         WHERE ai.owner = 'UNIPARC'
           AND aip.index_owner ='UNIPARC'
           AND ai.index_type = 'NORMAL'
           AND ai.partitioned = 'YES'
           AND TEMPORARY = 'N'
           -- AND SUBPARTITION_COUNT = 0
           AND ai.index_name = aip.index_name
      ORDER BY aip.partition_name asc;



      CURSOR subpartitioned_index_list
      IS
        SELECT ai.index_name, ai.table_name, ais.subpartition_name
          FROM sys.all_indexes ai, sys.all_ind_subpartitions ais
         WHERE ai.owner = 'UNIPARC'
           AND ais.index_owner ='UNIPARC'
           AND ai.index_type= 'NORMAL'
           AND ai.partitioned = 'YES'
           AND ai.TEMPORARY = 'N'
           AND ai.index_name = ais.index_name
      ORDER BY ais.subpartition_name asc;

     BEGIN

        FOR i IN index_list
        LOOP
            EXECUTE IMMEDIATE
            'ALTER INDEX ' || i.index_name || ' REBUILD TABLESPACE UNIIND';
        END LOOP;

        FOR p IN partitioned_index_list
        LOOP
            EXECUTE IMMEDIATE
            'ALTER INDEX ' || p.index_name || ' REBUILD PARTITION ' || p.partition_name;
        END LOOP;

        FOR s IN subpartitioned_index_list

        LOOP
            EXECUTE IMMEDIATE
            'ALTER INDEX ' || s.index_name || ' REBUILD SUBPARTITION ' || s.subpartition_name;
        END LOOP;

END rebuild_all_indexes;
/
set define on
