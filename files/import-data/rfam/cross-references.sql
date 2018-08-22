select distinct
    link.*,
    family.type
from database_link link
join family on family.rfam_acc = link.rfam_acc
