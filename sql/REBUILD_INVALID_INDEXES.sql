set define off

create or replace PROCEDURE        "REBUILD_INVALID_INDEXES"
IS
/***********************************************************************
/* Package Name:  REBUILD_INVALID_INDEXES
/* Author:        Steven Rosanoff
/* Date Created:  13/12/2012
/* Description:   This package rebuilds any indexes marked as unusable.
/* Automated:     NO. This procedure is called when the new release is
/*                prepared. The script that calls this procedure is:
/**********************************************************************/

    CURSOR index_list
    IS

      SELECT INDEX_NAME
      FROM sys.ALL_INDEXES
      WHERE status= 'UNUSABLE';

  BEGIN

      FOR ind IN index_list
      LOOP
        EXECUTE IMMEDIATE
        'ALTER INDEX ' || ind.index_name || ' REBUILD TABLESPACE UNIIND';
      END LOOP;

  END REBUILD_INVALID_INDEXES;
/
set define on
