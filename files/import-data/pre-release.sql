create table if not exists load_rnc_sequence_features (
    accession varchar(100) NOT NULL,
    taxid int not null,
    start int not null,
    stop int not null,
    feature_name varchar(50),
    metadata jsonb
);

create table if not exists load_rnc_secondary_structure (
    rnc_accession_id varchar(100),
    secondary_structure text,
    md5 varchar(32)
);
