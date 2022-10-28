include { crw } from './databases/crw'
include { ena } from './databases/ena'
include { ensembl } from './databases/ensembl'
include { expressionatlas } from './databases/expressionatlas'
include { five_s_rrnadb } from './databases/5srrnadb'
include { flybase } from './databases/flybase'
include { genecards_suite } from './databases/genecards_suite'
include { gtrnadb } from './databases/gtrnadb'
include { hgnc } from './databases/hgnc'
include { intact } from './databases/intact'
include { lncbase } from './databases/lncbase'
include { lncbook } from './databases/lncbook'
include { lncipedia } from './databases/lncipedia'
include { mirbase } from './databases/mirbase'
include { mirgenedb } from './databases/mirgenedb'
include { pdbe } from './databases/pdbe'
include { pirbase } from './databases/pirbase'
include { plncdb } from './databases/plncdb'
include { pombase } from './databases/pombase'
include { psicquic } from './databases/psicquic'
include { quickgo } from './databases/quickgo'
include { refseq } from './databases/refseq'
include { rfam } from './databases/rfam'
include { rgd } from './databases/rgd'
include { ribovision } from './databases/ribovision'
include { sgd } from './databases/sgd'
include { silva } from './databases/silva'
include { snodb } from './databases/snodb'
include { snorna_database } from './databases/snorna_database'
include { tarbase } from './databases/tarbase'
include { zfin } from './databases/zfin'
include { zwd } from './databases/zwd'

process build_context {
  memory '6GB'
  when { params.needs_taxonomy }

  output:
  path('context.db')

  """
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
  tar xvf new_taxdump.tar.gz
  mkdir taxdump
  mv *.dmp taxdump
  rnac context build taxdump context.db
  """
}

workflow parse_databases {
  emit: data
  main:

    build_context | set { context }

    Channel.empty() \
    | mix(
      crw(),
      five_s_rrnadb(),
      ena(),
      ensembl(),
      expressionatlas(),
      flybase(),
      genecards_suite(),
      gtrnadb(context),
      hgnc(),
      intact(),
      lncbase(),
      lncbook(),
      lncipedia(),
      mirbase(),
      mirgenedb(),
      pdbe(),
      pirbase(),
      plncdb(),
      pombase(),
      psicquic(),
      quickgo(),
      refseq(),
      rfam(),
      rgd(),
      ribovision(),
      sgd(),
      silva(context),
      snodb(),
      snorna_database(),
      tarbase(),
      zfin(),
      zwd(context),
    ) \
    | flatten \
    | set { data }
}
