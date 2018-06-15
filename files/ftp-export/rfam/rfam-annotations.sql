select
    hits.upi,
    hits.rfam_model_id,
    score,
    e_value,
    sequence_start,
    sequence_stop,
    model_start,
    model_stop,
    models.long_name
from rfam_model_hits hits
join rna_active active on active.upi = hits.upi
join rfam_models models on models.rfam_model_id = hits.rfam_model_id
order by hits.upi, hits.sequence_start, hits.rfam_model_id

