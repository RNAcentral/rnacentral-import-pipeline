COPY(
select upi from xref

where deleted = 'N'
and EXTRACT (DAY FROM (CURRENT_TIMESTAMP - timestamp)) < 7

) TO STDOUT (FORMAT CSV)
