set define off

create or replace PROCEDURE REBUILD_TABLE_INDEXES (table_id IN VARCHAR2)
IS

/***********************************************************************
/* Package Name:  REBUILD_TABLE_INDEXES
/* Author:        Steven Rosanoff
/* Date Created:  13/12/2012
/* Description:   The package allows the user to rebuild indexes
/*                on a specific table. You must specify a table name
/*                when the procedure is called.
/* Automated:     No. This procedure is run manually.
/**********************************************************************/


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
