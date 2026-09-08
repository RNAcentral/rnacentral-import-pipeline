# Training the r2dt "should show" model

`files/r2dt/should-show/model.joblib` is a `RandomForestClassifier` that
decides, for each R2DT secondary structure diagram, whether it's good enough
to display on the website. It is used by `rnac r2dt should-show compute`
(called from `workflows/r2dt.nf`) to label every newly-mapped URS.

All of the commands below live under `rnac r2dt should-show`
(`rnacentral_pipeline/cli/r2dt.py`, logic in
`rnacentral_pipeline/rnacentral/r2dt/should_show.py`). `rnac` is the CLI
entrypoint alias used elsewhere in this repo:

```bash
rnac() { uv run python -c 'from rnacentral_pipeline.cli import cli; cli()' "$@"; }
```

## The features

The model is trained on, per URS/model pairing:

- `source_index` — which R2DT pipeline produced the diagram (crw, ribovision,
  gtrnadb, rnase_p, rfam)
- `sequence_length`, `model_length`, `model_basepair_count`
- `diagram_bps`, `diagram_overlap_count`
- `diagram_sequence_length`, `diagram_model_length` — derived from the
  diagram's sequence/model start-stop coordinates

These come from joining `rna`, `r2dt_results` and `r2dt_models` in the
database (`should_show.fetch_modeled_data`); `should_show.infer_columns`
derives the two length columns and the source index from the raw start/stop
columns.

## The quickest path: retrain from the existing labelled corpus

`data/r2dt/should-show/labelled-corpus.csv` is a frozen, already-featured,
hand-labelled snapshot (see its header for the exact columns, and
`tests/rnacentral/r2dt/should_show_test.py`'s module docstring for how it was
produced). If you just want to retrain on the existing labels — e.g. after a
scikit-learn upgrade, or a tweak to the training code — you don't need a
database at all:

```bash
rnac r2dt should-show build-model \
    data/r2dt/should-show/labelled-corpus.csv \
    files/r2dt/should-show/model.joblib
```

`build-model` looks at the input CSV's header: if it has a `label` column
(as the labelled corpus does), it's loaded directly. Otherwise the input is
treated as a plain `urs,flag` list and `--db-url` (or `$PGDATABASE`) is
required to fetch the features live from the database.

Training does a stratified train/test split first, runs 5-fold cross
validation on the training split purely to estimate generalisation (F1 and
accuracy), then refits a fresh model on the *entire* training split and
scores that final model against the held-out test split. Both sets of
metrics are printed:

```
5-fold CV: f1=0.9551 (+/- 0.0089) accuracy=0.9485 (+/- 0.0101)
Test set: f1=0.9462 accuracy=0.9391
```

## Adding new labelled examples first

If you want to grow or refresh the labelled set before retraining:

1. **`fetch-data`** — build a plain features CSV for a list of URS (one per
   line), fetched from the database:

   ```bash
   rnac r2dt should-show fetch-data --db-url $PGDATABASE urs-list.csv features.csv
   ```

   This is also useful just to explore feature values before committing to
   a final training approach.

2. **`inspect-data`** — intended to fetch a richer sheet (link to the site,
   model source/name/SO term, plus an empty "Labeled Should show" column)
   for a person to fill in by hand, so it can be pasted into a spreadsheet
   with embedded SVGs for side-by-side comparison and labelled there.

   ```bash
   rnac r2dt should-show inspect-data --db-url $PGDATABASE urs-list.csv to-label.csv
   ```

3. **`convert-sheet`** — once the spreadsheet has been labelled and
   downloaded as CSV, convert it back into the `urs,flag` training format.
   Requires a `urs` column and a `Labeled Should show` column
   (`true`/`false`, case-insensitive; blank or unrecognised values are
   skipped and logged):

   ```bash
   rnac r2dt should-show convert-sheet labelled-sheet.csv new-labels.csv
   ```

   Append the result to `data/r2dt/should-show/training-data.csv`.

4. Re-run `fetch-data` over the full `training-data.csv` to build a
   fully-featured CSV, hand-add the `label` column (0/1, matching the flags
   in `training-data.csv`) to turn it into a corpus shaped like
   `labelled-corpus.csv`, then train as in the quickest-path section above.

## Running the trained model

`compute` applies a model to a batch of URS, fetching their features live
from the database — this is the production path used by
`workflows/r2dt.nf`:

```bash
rnac r2dt should-show compute --db-url $PGDATABASE \
    files/r2dt/should-show/model.joblib urs-list.csv should-show.csv
```

## After retraining

- `tests/rnacentral/r2dt/should_show_test.py` hardcodes the shipped model's
  accuracy and error counts (`EXPECTED_ACCURACY`, `WRONGLY_SHOWN`,
  `WRONGLY_HIDDEN`) as a change detector, and a frozen `prediction` column in
  `labelled-corpus.csv` for the prediction-stability test. Both need
  regenerating against the new model — see that test file's module
  docstring and the `corpus()` helper for how the checks are computed.
- Don't commit a retrained `model.joblib` without re-running
  `tests/rnacentral/r2dt/should_show_test.py` and updating those baselines
  to match.
