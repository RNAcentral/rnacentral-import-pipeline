select
    clan.clan_acc 'id',
    clan.description 'name',
    clan.comment 'description',
    count(distinct membership.rfam_acc) 'family_count'
from clan
join clan_membership membership
on
    membership.clan_acc = clan.clan_acc
group by clan.clan_acc
