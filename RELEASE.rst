Release Notes
=============

This is a summary of changes to the internals of the database and how we
process data. The format is based off of `Keep a Changelog
<http://keepachangelog.com/en/1.0.0/>`_. Restructured text is being used because
it provides more structure to the text.

Release 10
----------

Added
`````

- JSON schema based import.

  This will now import lncipedia data using our `JSON schema
  <https://github.com/RNAcentral/rnacentral-data-schema>`_. In theory everything
  we currently extract can be represented in it, however, it seems that not
  everything is. Notably, we will have fewer publications and no product
  information. On the plus side we will have coordinates for hg19 and hg38 for
  this database.

- A ``qa_status`` table.

  This table is where we will store all QA status. This is a generalization and
  cleanup of the ``rfam_problems`` field from ``rnc_rna_precomputed``. This is a
  cleaner, easier to query method of representing the same data. The old column
  is not yet deleted, but is not being populated.

Release 9
---------

Added
`````

- Imported RGD data.

  This import may be a one off task. One issue with this is that they do not
  provide their sequences on their FTP site. These sequences can be extracted
  via the gff files and chromosomes they provide, however I didn't want to add
  that complexity. I asked them to provide sequences on an on going basis, but
  they may have considered it a one off task. We will have to follow up with
  them on this for long term updates.

- Create Rfam annotation export. This is now officially part of our FTP export.

- Added ``rnc_genomic_mapping`` table. This table stores the infered locations
  for particular sequences.

- Added a ``has_coordinates`` column to ``rnc_rna_precomputed``. This column is
  meant to reflect the if a UPI/taxid pair has any known genomic locations. It
  summarizes if there any entries in ``rnc_coordinates`` and
  ``rnc_genomic_mapping``. It defaults to ``false`` and isn't currently updated
  for entries where taxid is null.

Changed
```````

- GENCODE accessions are now namespaced with ``GENCODE:`` prefixed to the
  corresponding Ensembl id to produce the GENCODE accession.

- All ENA parsing is done in python. This has several effects:

  1. We do not write out both the original ENA entry and the expert database
     entry for sequences that come from an expert database. We modify the
     original ENA entry and only produce data for the expert database.

  2. We determine what sequences come from expert databases via the `XREF
     service <https://www.ebi.ac.uk/ena/browse/xref-service-rest>`_ instead of
     parsing ``DR`` lines. We think XREF service more reliable than ``DR``
     lines as it seems to be updated automatically and much faster, while
     historically ``DR`` lines do not.

  3. This may change how some columns are populated. The ``locus_tag`` field
     will now be populated with data from the ``locus`` qualifier from the EMBL
     file. For TAIR the value in the ``locus`` qualifier was used for
     ``optional_id`` in the database.

  4. Fields in ``rnc_accessions`` which used to be the empty string are now
     ``NULL``.

  5. The data in the ``note`` and ``db_xref`` are now JSON data structures. For
     ``note`` the data structure will have ``text`` and ``ontology`` fields. The
     ``ontology`` field is a mapping from a selected ontology (GO, SO, ECO) to a
     list of the terms that have been extracted. The ``text`` field contains a
     list of all strings of all text that could not be recognized as an from an
     ontology. As an example:

     .. code:: json

       {
         "ontology": [
             "ECO:0000202",
             "GO:0030533",
             "SO:0000253"
         ],
         "text": [
             "Covariance Model: Bacteria; CM Score: 87.61",
             "Legacy ID: chr.trna3-GlyGCC"
         ]
       }

     The ``db_xref`` field will contain data from the ``db_xref`` qualifier of
     the record. As well as data from the comment field. Additionally there with
     will be a ``db_xref`` field which contains the information from the ``DR``
     lines, excluding the MD5 xref. The keys will be the database in upper case,
     and the values will be a list of primary id and secondary id will be
     ``null`` or the secondary id for the database, if present.

     .. code:: json

       {
         "ena_refs": {
           "WORMBASE": ["WBGene00196009", "F19C6.9"],
           "BIOSAMPLE": ["SAMEA3138177", null]
         }
       }

  6. More publications are extracted. This will pull publication data from the
     ``experiment`` qualifier. Sometimes this qualifier contains a string like
     ``PMID: 10, 12``, in which case those numbers are extracted and the used to
     lookup the publication information from Europe PMC.

  7. The ``anticodon`` field is more filled out. This will try to pull the
     sequence from the ``anticodon`` qualifier, or the ``gene`` qualifier, or
     from the ``note`` qualifier, using a variety of patterns.

- Database mappings are generated by this pipeline instead of the webapp. There
  is only 1 intentional change from the previous format, and that is adding
  chain ids to PDBe mappings as per `rnacentral-webcode/288
  <https://github.com/RNAcentral/rnacentral-webcode/issues/288>`_.

- Other exports are now done by this pipeline and are expected to produce the
  same results. They are:

  - MD5 mapping export
  - FASTA export
