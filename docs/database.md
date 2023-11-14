# Database design

This is some documentation on the structure of the [RNAcentral] database.
This is a cleaned up version of the output of `\d tablename` with some notes on the meaning of the columns and some general notes on the table.

# Tables

## Notes on indexing

There are several ways to assign coordinates in a sequence.
If this is new to you please read this nice [blog post](https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/).
In our data we have two types, the `1-based, half open` and `0-based, half open` ones.
The first is what many genomic coordinates use, second is what most programming languages use.
If a column is called an `index` it is the `0-based, half open`, while a `position` is a `1-based, half open index`.
Note that bugs have happened and sometimes data is mixed or incorrect, we are working to fix this.

## `go_term_annotations`

## `go_term_publication_map`

## `rfam_go_terms`

A table that maps [Rfam] families to the [Gene Ontology (GO)](https://geneontology.org/) terms that are assigned to each [Rfam] family.
The terms may be filtered slightly to provide "better" terms for [RNAcentral] usage.

|       Column       |  Type    | Description                                   |
|--------------------|----------|-----------------------------------------------|
| `rfam_go_term_id`  | integer  | A id for this row.                            |
| `rfam_model_id`    | text     | A reference to the `rfam_model.id` field.     |
| `ontology_term_id` | text     | A reference to the `ontology_terms.id` field. |

## `ontology_terms`

## `rna`

This table stores all sequences [RNAcentral] has ever seen.
No data is ever deleted from this table, even in the case of bugs.
No data should ever be deleted as unpredictable things may occur, also as we promise to keep all sequences public this could cause issues if the ids are public.
This serves as where the [`URS`](https://rnacentral.org/help#rnacentral-identifiers) ids are created and stored.
These ids are what we consistently use to identify a sequence to the public.
However, internally we are not consistent with the naming of columns.
So there are columns called `urs`, `upi` and sometimes `rna_id` which all may be the same thing.
This is for historical reasons, basically the database was built of an existing one, [UniProt], which used `upi`, and not all columns were renamed.
Then future work was not consistent in the name selected.
Currently, columns should be named `urs` if it contains a `urs` and `urs_taxid` if it contains a `urs_taxid`.

Another fun note is that this table contains two sequence columns `seq_short` and `seq_long`.
Only one will ever be populated.
There are two columns because the database used to be oracle and supposedly oracle prefers if large and small strings are stored in different columns.
The sequences are stored as DNA.
In general sequences should contain no more than 10% uncertainty characters, N, Y, etc.
For details of what those codes mean see: <https://www.bioinformatics.org/sms/iupac.html>

| Column      |            Type             | Description                                                                                                              |
|-------------|-----------------------------|--------------------------------------------------------------------------------------------------------------------------|
| `id`        | integer                      | A numeric ID for each sequence. Not generally used elsewhere.                                                            |
| `upi`       | text                        | This is the `URS` for the given sequence.                                                                                |
| `timestamp` | timestamp without time zone | When the sequence was created                                                                                            |
| `userstamp` | text                        | The user that created it, never important.                                                                               |
| `crc64`     | text                        | The [CRC](https://en.wikipedia.org/wiki/Cyclic_redundancy_check) checksum for the sequence. Not generally used anywhere. |
| `len`       | integer                     | The length of the sequence.                                                                                              |
| `seq_short` | text                        | The sequence as DNA, if it is under 4,000 nts.                                                                           |
| `seq_long`  | text                        | The sequence as DNA, if it is over 4,000 nts.                                                                            |
| `md5`       | text                        | The MD5 hash of the sequence.                                                                                            |

## `rnc_rna_precomputed`

## `rnc_accessions`

This table represents an entry the metadata about a sequence in a particular database.
The columns are basically the possible fields from a EMBL formatted file, because very early [RNAcentral] was based on importing from ENA.
This is no longer true, but the fundamental structure remains.

|      Column       | Type     | Description                                                                                                                        |
|-----------------  |-------   |--------------                                                                                                                      |
|  `id`             | integer  | A numeric id for each accession, not used anywhere.                                                                                |
|  `accession`      | text     | A unique accession for each sequence. Generally this is formatted as `DATABASE:DATABASE-id`.                                       |
|  `parent_ac`      | text     |
|  `seq_version`    | integer  | Version of the sequence. For sequences like `ENST00000001.1`, this is `1`. This should match the sequence version field in `xref`. |
|  `feature_start`  | integer  |
|  `feature_end`    | integer  |
|  `feature_name`   | text     |
|  `ordinal`        | integer  |
|  `division`       | text     |
|  `keywords`       | text     |
|  `description`    | text     |
|  `species`        | text     |
|  `organelle`      | text     |
|  `classification` | text     |
|  `project`        | text     |
|  `is_composite`   | text     |
|  `non_coding_id`  | text     |
|  `database`       | text     | The name of the database which this accession is for. This is actually a reference to `rnc_database.descr` field.                  |
|  `external_id`    | text     |
|  `optional_id`    | text     |
|  `common_name`    | text     |
|  `allele`         | text     |
|  `anticodon`      | text     | The anticodon for this sequence, if any. Not used and is replaced by the anticodon entries in `rnc_sequence_features`. |
|  `chromosome`     | text     | The chromosome this accession is found on, if any. |
|  `experiment`     | text     |
|  `function`       | text     |
|  `gene`           | text     | Name of the gene the database assigns to this sequence                                                                             |
|  `gene_synonym`   | text     | A comma separated list of synonyms for the gene.                                                                                   |
|  `inference`      | text     | The phylogeny that the providing database inferred. This only comes from two databases, [SILVA] and [MGnify] right now.            |
|  `locus_tag`      | text     |
|  `map`            | text     |
|  `mol_type`       | text     |
|  `ncrna_class`    | text     |
|  `note`           | text     | A catch all field for data from the database. Newer accessions are JSON formatted, older ones are ad-hoc.                          |
|  `old_locus_tag`  | text     |
|  `operon`         | text     |
|  `product`        | text     |
|  `pseudogene`     | text     |
|  `standard_name`  | text     |
|  `db_xref`        | text     |
|  `rna_type`       | text     | The [Sequence Ontology (SO)](http://www.sequenceontology.org/) term for this sequence.                                             |

## `rnc_secondary_structure`

This is an outdated table, scheduled for removal/simplification.
It was used to represent secondary structures provided by some expert databases.

## `rnc_secondary_structure_layout`

This table stores the secondary structures [RNAcentral] knows about.
These secondary structures are always computed by [R2DT].
This table does not include the SVG layouts the program produces, as that it stored on the file system.
The table *assumes* that a sequence may only have 1 secondary structure, this is likely to have to change eventually.

One thing to now about this table is it includes a pair of columns `inferred_should_show` and `assigned_should_show`, which are used to decide if the secondary structure should be shown to users.
We run [R2DT] on all sequences, but it doesn't always do a good job.
To detect cases where R2DT doesn't do well, we have some logic in the pipeline to compute if a secondary structure is worth storing.

|       Column           |       Type       | Description                                                                    |
|------------------------|------------------|------------------------------------------------------------------------------  |
| `urs`                  | text             | The `URS` that the secondary structure is assigned to.                         |
| `secondary_structure`  | text             | A [dot-bracket] string of the secondary structure.                             |
| `overlap_count`        | integer          | Number of overlaps in the diagram, as computed by [R2DT].                      |
| `basepair_count`       | integer          | Number of base pairs in the diagram, as computed by [R2DT].                    |
| `model_start`          | integer          | The start index in the model matched by the model.                             |
| `model_stop`           | integer          | The end index in the model matched by the model.                               |
| `sequence_start`       | integer          | The start index in the sequence matched by the model.                          |
| `sequence_stop`        | integer          | The end index in the sequence matched by the model.                            |
| `sequence_coverage`    | double precision | Fraction of the sequence that is covered by the secondary structure.           |
| `model_id`             | integer          | Reference to the matching model.                                               |
| `id`                   | integer           | A numeric id for the row.                                                      |
| `inferred_should_show` | boolean          | A computed flag to show if this secondary structure should be shown to user.   |
| `model_coverage`       | double precision | Fraction of the model that is matched by this secondary structure.             |
| `stk`                  | text             | The contents of the STK file computed by [R2DT].                               |
| `assigned_should_show` | boolean          | The manually assigned should_show value.                                       |

## `rnc_taxonomy`

This stores the taxonomic information [RNAcentral] knows about.
This is effectively a mirror of data from the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) and is used to create foreign keys and ensure taxids we process are known.
Each row is a taxon and all taxons should have a row.

|  Column       |  Type   | Description                                                   |
|---------------|---------|---------------------------------------------------------------|
| `id`          | integer | The NCBI taxon id of the taxon.                               |
| `name`        | text    | The name of the taxon.                                        |
| `lineage`     | text    | A `; ` separated string of the lineage of the taxon.          |
| `aliases`     | text[]  | A list of known aliases for the taxon.                        |
| `replaced_by` | integer | The id of the taxon which replaced this one, if any.          |
| `common_name` | text    | The common name, if any for this taxon.                       |
| `is_deleted`  | boolean | A flag if this taxon was deleted.                             |

## `protein_info`

This table is meant to provide [RNAcentral] with just enough information to display some summaries of proteins from [UniProt].
This doesn't have much information, and it should stay that way, we only need to be able to display a few labels and provide links.

|       Column        |  Type  | Description |
|---------------------|--------|-----------|
| `protein_accession` | text   |
| `description`       | text   |
| `label`             | text   |
| `synonyms`          | text[] | A list of synonyms for the protein. |

[MGnify]: https://www.ebi.ac.uk/metagenomics
[R2DT]: http://r2dt.bio/
[RNAcentral]: https://rnacentral.org/
[Rfam]: https://rfam.org/
[SILVA]: https://www.arb-silva.de/
[UniProt]: https://www.uniprot.org/
[dot-bracket]: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/io/rna_structures.html
