select
  rfam_acc
from family  f
where
  (f.type like '%riboswitch%' or f.type like 'Gene%' or f.type like '%IRES%' or f.type like '%thermoregulator%')
  and (f.type not like '%lncRNA%')
