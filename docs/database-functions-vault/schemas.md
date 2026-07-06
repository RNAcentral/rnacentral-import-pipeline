# Schemas

The four schemas exported into [`database_functions/`](../../database_functions):

- [[rnc_update]] — orchestrates a release load end-to-end.
- [[rnc_load_rna]] — loads new sequences into the `rna` table (UPI assignment).
- [[rnc_load_xref]] — full-release xref load via partition build-and-exchange.
- [[release]] — release-table query and mutation helpers.

See [[external-dependencies]] for schemas referenced but not present here.
Back to [[_index]].
