include { circatlas } from './databases/circatlas'
include { circpedia } from './databases/circpedia'
include { crw } from './databases/crw'
include { ena } from './databases/ena'
include { ensembl } from './databases/ensembl'
include { evlncrnas } from './databases/evlncrnas'
include { expressionatlas } from './databases/expressionatlas'
include { five_s_rrnadb } from './databases/5srrnadb'
include { flybase } from './databases/flybase'
include { genecards_suite } from './databases/genecards_suite'
include { gtrnadb } from './databases/gtrnadb'
include { hgnc } from './databases/hgnc'
include { intact } from './databases/intact'
include { japonicusdb } from './databases/japonicusdb'
include { lncbase } from './databases/lncbase'
include { lncbook } from './databases/lncbook'
include { lncipedia } from './databases/lncipedia'
include { mgi } from './databases/mgi'
include { mgnify } from './databases/mgnify'
include { mirbase } from './databases/mirbase'
include { mirgenedb } from './databases/mirgenedb'
include { mirtrondb } from './databases/mirtrondb'
include { modomics } from './databases/modomics'
include { noncode } from './databases/noncode'
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
include { tmrna } from './databases/tmrna'
include { zfin } from './databases/zfin'
include { zwd } from './databases/zwd'
include {ribocentre } from './databases/ribocentre'

process build_context {
  memory '6GB'
  errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  maxRetries 5

  output:
  path('context.db')

  script:
  """
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
  tar xvf new_taxdump.tar.gz
  mkdir taxdump
  mv *.dmp taxdump
  rnac context build taxdump context.db
  """
}

workflow parse_databases {
  main:

    if (params.get('needs_taxonomy', false)) {
      build_context | set { context }
    }
    else {
      channel.empty() | set { context }
    }

    channel.empty() \
    | mix(
      circatlas(),
      circpedia(),
      crw(),
      five_s_rrnadb(),
      ena(),
      ensembl(),
      evlncrnas(),
      expressionatlas(),
      flybase(),
      genecards_suite(),
      gtrnadb(context),
      hgnc(),
      intact(),
      japonicusdb(),
      lncbase(),
      lncbook(),
      lncipedia(),
      mgi(),
      mirbase(),
      mgnify(),
      mirgenedb(),
      mirtrondb(),
      modomics(),
      noncode(),
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
      ribocentre(),
      sgd(),
      silva(context),
      snodb(),
      snorna_database(),
      tarbase(),
      tmrna(),
      zfin(),
      zwd(context),
    ) \
    | flatten \
    | set { data }
  emit: data
}
