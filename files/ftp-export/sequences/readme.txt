
===================================================================
RNAcentral Sequence Data
===================================================================

This directory contains sequences with RNAcentral ids in FASTA format.

* rnacentral_active.fasta.gz
Current set of sequences that are present in at least one expert database.

* rnacentral_species_specific_ids.fasta.gz
Current set of sequences that are present in at least one expert database using
the species specific URS ID's.

* rnacentral_inactive.fasta.gz
All RNAcentral sequences that used to be present in one or more expert database
but don't have any current cross-references.

* by-database/
Active sequences split into individual files, one per expert database
(e.g. mirbase.fasta).

* by-species/
Species-specific URS sequences split into individual files, one per species 
(e.g. homo_sapiens.fasta.gz).

* example.txt
A small example file showing the format of rnacentral_active.fasta.gz
and rnacentral_inactive.fasta.gz.
            