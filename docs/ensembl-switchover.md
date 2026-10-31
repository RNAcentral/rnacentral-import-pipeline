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

## Assemblies: review before the first production run

The import now keeps each species on the exact assembly RNAcentral selects today
(`files/import-data/ensembl/current-assemblies.txt`, the `selected_genome = true`
rows of the release-28 public database). The `ensembl_assembly` table needs a
proper review and clean-up. Checked against species.json on 2026-10-30:

- **Regenerate `current-assemblies.txt` from production** before the run. The
  public copy may lag behind it.
- **Renamed assemblies.** 34 genomes keep the same GCA but Ensembl's
  files now give a different assembly name. Examples: Plasmodium falciparum
  (`GCA000002765v3` → `GCA_000002765`), P. vivax, vervet (`ChlSab1.1` →
  `Chlorocebus_sabeus_1.1`), several Aspergillus and Candida, and names with a
  new `gb` suffix. Folsomia candida could not be checked: its GTF is missing
  from the FTP. Unless the new names are mapped to RNAcentral's existing ones,
  the import creates a second row and duplicate coordinates for each.
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
