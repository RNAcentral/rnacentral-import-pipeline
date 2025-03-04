select
    rfam_acc,
    rfam_id,
    f.type,
    fr.type as 'sequence_type',
    f.description,
    rs.ncbi_id,
    fr.rfamseq_acc,
    fr.seq_start,
    fr.seq_end,
    group_concat(distinct concat(dl.db_id,":",dl.db_link)) as dbxrefs,
    group_concat(distinct concat("PMID:",fl.pmid)) as PMIDS,
    rs.version,
    tx.species,
    tx.tax_string
from family f
join full_region fr using (rfam_acc)
join rfamseq rs using (rfamseq_acc)
join taxonomy tx using (ncbi_id)
join database_link dl using (rfam_acc)
join family_literature_reference fl using (rfam_acc)
where
    f.rfam_acc = @family
    and fr.type != 'seed'
    and fr.is_significant = true
group by rfam_acc, rfam_id, f.type, fr.type, f.description,rs.ncbi_id, fr.rfamseq_acc, fr.seq_start, fr.seq_end, '', rs.version, tx.species, tx.tax_string
;
