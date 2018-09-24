COPY (
select
  descr
from rnc_database
where
  alive = 'Y'
);
