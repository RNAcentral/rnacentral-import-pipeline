process psicquic {
  output:
  path('*.{csv,parquet}')

  when: params.databases.psicquic?.run

  script:
  """
  echo -e "#ID(s) interactor A\tID(s) interactor B\tAlt. ID(s) interactor A\tAlt. ID(s) interactor B\tAlias(es) interactor A\tAlias(es) interactor B\tInteraction detection method(s)\tPublication 1st author(s)\tPublication Identifier(s)\tTaxid interactor A\tTaxid interactor B\tInteraction type(s)\tSource database(s)\tInteraction identifier(s)\tConfidence value(s)\tExpansion method(s)\tBiological role(s) interactor A\tBiological role(s) interactor B\tExperimental role(s) interactor A\tExperimental role(s) interactor B\tType(s) interactor A\tType(s) interactor B\tXref(s) interactor A\tXref(s) interactor B\tInteraction Xref(s)\tAnnotation(s) interactor A\tAnnotation(s) interactor B\tInteraction annotation(s)\tHost organism(s)\tInteraction parameter(s)\tCreation date\tUpdate date\tChecksum(s) interactor A\tChecksum(s) interactor B\tInteraction Checksum(s)\tNegative\tFeature(s) interactor A\tFeature(s) interactor B\tStoichiometry(s) interactor A\tStoichiometry(s) interactor B\tIdentification method participant A\tIdentification method participant B" > raw.tsv

  for svc in ${params.databases.psicquic.services.join(' ')}; do
    curl -sf "https://www.ebi.ac.uk/Tools/webservices/psicquic/\$svc/webservices/current/search/query/rnacentral?format=tab27" >> raw.tsv
  done

  rnac psicquic parse raw.tsv .
  """
}
