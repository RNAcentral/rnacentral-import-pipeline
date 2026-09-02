// The divisions of the (old) Ensembl import, and whether any of them are set to
// run. Derived on demand rather than stored in a param: entry scripts like
// import-data.nf are run directly, so main.nf is not evaluated and anything
// computed there is simply absent.

def ensembl_divisions() {
  ['plants', 'fungi', 'protists', 'metazoa', 'vertebrates']
}

def any_ensembl_division_runs() {
  ensembl_divisions().any { d -> params.databases.ensembl[d]?.run }
}
