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
   looks right. Then set `alive = 'N'` on those four `rnc_database` rows (ids
   31, 34, 35, 36). A retirement by hand never refreshes their
   `num_sequences`, so they would stay in the FTP export's database list and
   each would get an empty FASTA file.
2. **Run `files/precompute/deactivate-stale.sql` after precompute.** Release-mode
   precompute doesn't reselect a sequence whose only xrefs were the retired
   division ones, so its precompute row would stay active.

## Delta loading

The first release that loads ENSEMBL has no stored manifest, so every geneset is parsed and the
load is FULL. That release seeds `pipeline_tracking_import` for ENSEMBL, ENSEMBL_GENCODE and TAIR, and
`pipeline_tracking_import_files` for the genesets. Later releases fetch only genesets whose
`md5sum.txt` entries moved. `--force_full_import` reparses everything.

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
8. **Check the web code no longer reads `ensembl_assembly.subdomain`.** The
   import stopped writing `division` and `subdomain`, since Ensembl dropped its
   divisions and subdomains. `division` was already nullable. The first run
   makes `subdomain` nullable as well (`000__assemblies.sql`). Existing rows
   keep their old values, and new rows leave both columns blank. Neither is
   used for xrefs, so the division retirement in step 1 is a separate job. Once
   nothing reads them, drop both columns.

9. **Load compara separately when its files are published.** The main run leaves
   `ensembl.compara.run` false, so it skips compara. Once the files exist, run
   with `ensembl.run = false` and `ensembl.compara.run = true`, with the release
   on. It fetches and parses only compara, then loads it in the release. It
   has to come after the first Ensembl load, since it maps transcripts to
   `ENSEMBL` xrefs (`dbid = 25`).

## Ensembl coding (`ENSEMBL_CODING`, branch `update/Ensembl-mrna`)

### Before the test run

9. **Check that dbid 60 is free in the target database.** The public copy tops
   out at 59, but production may be ahead of it.
10. **Create the database row and its xref partitions.** `xref` is partitioned
    by `dbid`, then by `deleted`. The legacy `xref_insert_trigger` isn't used.
    This copies the layout of dbid 59:

    ```sql
    BEGIN;
    INSERT INTO rnacen.rnc_database (id, timestamp, userstamp, descr, alive, display_name, full_descr, description, url, reference)
    VALUES (60, now(), 'RNACEN', 'ENSEMBL_CODING', 'Y', 'Ensembl coding', 'Ensembl coding',
      ': protein-coding mRNAs and their 5'' and 3'' UTRs from Ensembl',
      'https://www.ensembl.org',
      '[{"title":"Ensembl 2025","authors":"Dyer SC, Austine-Orimoloye O, Azov AG, et al.","journal":"Nucleic Acids Research, 2025","pubmed_id":"39656687"}]'::jsonb);

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

## Assemblies: review before the first production run

The import now keeps each species on the exact assembly RNAcentral selects today
(`files/import-data/ensembl/current-assemblies.txt`, the `selected_genome = true`
rows of the release-28 public database). The `ensembl_assembly` table needs a
proper review and clean-up. Checked against species.json on 2026-10-30:

- **Regenerate `current-assemblies.txt` and `assembly-ids.tsv` from production**
  before the run. The public copy may lag behind it. `assembly-ids.tsv` is
  `gca_accession<TAB>assembly_id` for every `ensembl_assembly` row, selected rows
  first: `copy (select gca_accession, assembly_id from ensembl_assembly where
  gca_accession <> '' order by 1, selected_genome desc, 2) to stdout`.
- **Renamed assemblies: handled.** 103 genomes keep their GCA but Ensembl now
  names them differently (Plasmodium `GCA000002765v3` → `GCA_000002765`, Aedes
  `AaegL5` → `AaegL5.0`), and 14 short names (`v1.0`, `EF1`) belong to different
  genomes. The import takes each genome's id from `assembly-ids.tsv` by exact GCA,
  keeps it for a RefSeq copy of a held GenBank assembly, and gives a clashing name
  a GCA-based id (`GCA014634605v1`).
- **Versions Ensembl no longer offers.** 10 genomes only have a different
  version in species.json now (e.g. Amphiprion percula `.1` → `.2`), so they
  move assembly and their old coordinates stay behind.
- **Leftover rows.** Old assemblies keep their coordinates after a species moves
  on, e.g. Plasmodium's `ASM276v2` (1,499 coordinates, not selected). Nothing
  removes them, because the clean-up only drops assemblies with no coordinates.
- **The legacy list.** 1,252 of the 2,119 entries in `legacy-assemblies.txt` have
  no row in RNAcentral, and 98 have several versions or none selected.
- **Wheat.** RNAcentral's selected wheat assembly (PGSBv2.0) has no coordinates.
  Decide whether IWGSC RefSeq is the one to use.
- **Newer references not taken up.** Zebrafish (GRCz12tu) and Arabidopsis
  (TAIR12) have newer genomes with better UTR coverage. Moving a species is a
  whole-species decision, ncRNAs included.

## Decisions to confirm

8. **Description weight for Ensembl bacteria and archaea is 10.** Ensembl never
   imported them before, so there was no existing weight. It is set low so they
   don't take over descriptions from sources that already win. See
   `rnacentral/genes/random_forest/data.py`.
