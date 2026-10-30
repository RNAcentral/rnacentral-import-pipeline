# Ensembl switch-over: manual steps

These steps happen after the code on `update/Ensembl-import` has passed a test run.
Ensembl is then one `ENSEMBL` database, built from
`https://ftp.ebi.ac.uk/pub/ensemblorganisms/species.json`, and it replaces the
Plants, Metazoa, Fungi and Protists divisions.

## In the release that first loads the unified ENSEMBL

1. **Retire the old division xrefs.** Nothing in the pipeline retires them. If
   you skip this, every plant, metazoan, fungal and protist sequence shows up
   twice: once under its old division and once under `ENSEMBL`. Retire
   `ENSEMBL_PLANTS`, `ENSEMBL_METAZOA`, `ENSEMBL_FUNGI` and `ENSEMBL_PROTISTS`
   in the same release. Count the rows first and commit only if the count
   looks right.
2. **Run `files/precompute/deactivate-stale.sql` after precompute.** Release-mode
   precompute doesn't reselect a sequence whose only xrefs were the retired
   division ones, so its precompute row would stay active.

## Database clean-up

3. **Drop `ensembl_coordinate_systems`.** Nothing loads it any more and the web
   team confirmed nothing reads it. Also drop `load_coordinate_info` if it is
   still in the database.
4. **Review assemblies that genome mapping skips.** Genome mapping now takes
   genomes from species.json and logs a count of `ensembl_assembly` rows with
   `selected_genome = true` that aren't listed there. Decide whether to set
   those rows to `false`.

## Search

5. **Ask EBI Search to add `kingdom` as a facet** on the RNAcentral index. The
   export writes the field (bacteria, archaea, viruses, plants, fungi,
   vertebrates, metazoa, environmental or protists, the same split as the FTP
   folders), but facets and their labels are set in EBI Search's configuration.
   Until they add it, the field can be searched but won't appear as a facet.

## Website (`rnacentral-webcode`, branch `feature/ensembl-unified-database`)

6. **Deploy the webcode branch in the same release.** It archives the four
   division entries, redirects their pages to the Ensembl page, describes Ensembl
   as covering every kingdom, and adds the kingdom facet.
7. **Set the Ensembl `version`** in `portal/config/expert_databases.py`. It still
   says `116`, and species.json has no release number. One option is the
   species.json `last_updated` date.

## Ensembl coding (`ENSEMBL_CODING`, branch `update/Ensembl-mrna`)

### Before the test run

9. **Check that dbid 60 is free in the target database.** The public copy tops
   out at 59, but production may be ahead of it.
10. **Create the database row and its xref partitions.** `xref` is partitioned
    by `dbid`, then by `deleted`. The legacy `xref_insert_trigger` isn't used.
    This copies the layout of dbid 59:

    ```sql
    BEGIN;
    INSERT INTO rnacen.rnc_database (id, timestamp, userstamp, descr, alive, display_name, full_descr)
    VALUES (60, now(), 'RNACEN', 'ENSEMBL_CODING', 'Y', 'Ensembl coding', 'Ensembl coding');

    CREATE TABLE rnacen.xref_p60 PARTITION OF rnacen.xref FOR VALUES IN (60) PARTITION BY LIST (deleted);
    CREATE TABLE rnacen.xref_p60_deleted PARTITION OF rnacen.xref_p60 FOR VALUES IN ('Y');
    CREATE TABLE rnacen.xref_p60_not_deleted PARTITION OF rnacen.xref_p60 FOR VALUES IN ('N');

    DO $$
    DECLARE t text; flag text;
    BEGIN
      FOR t, flag IN VALUES ('xref_p60_deleted', 'Y'), ('xref_p60_not_deleted', 'N') LOOP
        EXECUTE format('ALTER TABLE rnacen.%I ADD CONSTRAINT %I CHECK (dbid = 60 AND deleted = %L)', t, t || '_check', flag);
        EXECUTE format('ALTER TABLE rnacen.%I ADD CONSTRAINT %I FOREIGN KEY (created) REFERENCES rnacen.rnc_release(id)', t, t || '_fk1');
        EXECUTE format('ALTER TABLE rnacen.%I ADD CONSTRAINT %I FOREIGN KEY (dbid) REFERENCES rnacen.rnc_database(id)', t, t || '_fk2');
        EXECUTE format('ALTER TABLE rnacen.%I ADD CONSTRAINT %I FOREIGN KEY (last) REFERENCES rnacen.rnc_release(id)', t, t || '_fk3');
        EXECUTE format('ALTER TABLE rnacen.%I ADD CONSTRAINT %I FOREIGN KEY (urs) REFERENCES rnacen.rna(urs)', t, t || '_fk4');
        EXECUTE format('CREATE INDEX %I ON rnacen.%I (ac)', t || '$ac', t);
        EXECUTE format('CREATE INDEX %I ON rnacen.%I (created)', t || '$created', t);
        EXECUTE format('CREATE INDEX %I ON rnacen.%I (dbid)', t || '$dbid', t);
        EXECUTE format('CREATE UNIQUE INDEX %I ON rnacen.%I (id)', t || '$id', t);
        EXECUTE format('CREATE INDEX %I ON rnacen.%I (last)', t || '$last', t);
        EXECUTE format('CREATE INDEX %I ON rnacen.%I (taxid)', t || '$taxid', t);
        EXECUTE format('CREATE INDEX %I ON rnacen.%I (urs)', t || '$upi', t);
      END LOOP;
    END $$;
    COMMIT;
    ```

11. **Add the three relationship types the UTR links use.** Run these outside a
    transaction:

    ```sql
    ALTER TYPE rnacen.related_sequence_relationship ADD VALUE IF NOT EXISTS 'five_prime_utr';
    ALTER TYPE rnacen.related_sequence_relationship ADD VALUE IF NOT EXISTS 'three_prime_utr';
    ALTER TYPE rnacen.related_sequence_relationship ADD VALUE IF NOT EXISTS 'mrna';
    ```

12. **Run the import with only `databases.ensembl_coding.run = true`.** It already
    limits itself to the species in `databases.ensembl_coding.species`.

### After the test run

13. **Check the counts:** one mRNA per Ensembl `protein_coding` transcript on the
    primary assembly, each with up to two UTRs, and links both ways in
    `rnc_related_sequences`. UTRs under 10 nt are dropped.
14. **IGV relies on run order.** For human and mouse, the coordinate export
    replaces genome mapping's `<species>.<assembly>.ensembl.gff3.gz` with mRNA gene
    models. That only works because the export runs after genome mapping.
15. **Move the mRNA test fixture into the `test-data` submodule.** The branch
    predates the submodule.
16. **Still to decide:** whether mRNAs go into text search, sequence search and
    the FTP FASTA and ID-mapping exports, and how the website shows UTR links and
    tags (web team).

## SGD coding (`SGD_CODING`, branch `update/SGD-coding`)

17. **Create the database row and xref partitions for dbid 61.** Run item 10's SQL
    with `61`, `'SGD_CODING'` and `'SGD coding'` in place of `60`,
    `'ENSEMBL_CODING'` and `'Ensembl coding'`. Item 11's link types are shared.
18. **Run the import with `databases.sgd_coding.run = true`.** It reads SGD's GFF,
    which carries the Pelechano 2013 transcript ends and the genome sequence.
    Dubious ORFs are left out.

## Decisions to confirm

8. **Description weight for Ensembl bacteria and archaea is 10.** Ensembl never
   imported them before, so there was no existing weight. It is set low so they
   don't take over descriptions from sources that already win. See
   `rnacentral/genes/random_forest/data.py`.
