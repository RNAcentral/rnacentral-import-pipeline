COPY(
select upi from xref

where dbid in (11, 16)
and deleted = 'N'
and EXTRACT (DAY FROM (CURRENT_TIMESTAMP - timestamp)) < 7

) TO STDOUT (FORMAT CSV)
