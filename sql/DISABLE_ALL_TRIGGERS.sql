set define off

create or replace PROCEDURE        "DISABLE_ALL_TRIGGERS"
IS
/***********************************************************************
/* Package Name:  DISABLE_ALL_TRIGGERS
/* Author:        Steven Rosanoff
/* Date Created:  26/10/2012
/* Description:   This package enables all triggers belonging to UNIPARC.
/* Automated:     NO
/**********************************************************************/

    CURSOR trigger_list
    IS
      SELECT trigger_name

        FROM sys.all_triggers
       WHERE owner = 'UNIPARC';

  BEGIN

      FOR trig IN trigger_list
      LOOP
        EXECUTE IMMEDIATE
        'ALTER TRIGGER ' || trig.trigger_name || ' DISABLE';
      END LOOP;

  END DISABLE_ALL_TRIGGERS;
/
set define on
