include { five_s_rrnadb } from './databases/5srrnadb'
include { ena } from './databases/ena'
include { ensembl } from './databases/ensembl'
include { ensembl_fungi } from './databases/ensembl_fungi'
include { ensembl_metazoa } from './databases/ensembl_metazoa'
include { ensembl_plants } from './databases/ensembl_plants'
include { ensembl_protists } from './databases/ensembl_protists'
include { flybase } from './databases/flybase'
include { genecards_suite } from './databases/genecards_suite'
include { gtrnadb } from './databases/gtrnadb'
include { intact } from './databases/intact'
include { lncbase } from './databases/lncbase'
include { lncbook } from './databases/lncbook'
include { lncipedia } from './databases/lncipedia'
include { mirbase } from './databases/mirbase'
include { mirgenedb } from './databases/mirgenedb'
include { pdbe } from './databases/pdbe'
include { pombase } from './databases/pombase'
include { quickgo } from './databases/quickgo'
include { refseq } from './databases/refseq'
include { rfam } from './databases/rfam'
include { rgd } from './databases/rgd'
include { sgd } from './databases/sgd'
include { silva } from './databases/silva'
include { snodb } from './databases/snodb'
include { snorna_database } from './databases/snorna_database'
include { tarbase } from './databases/tarbase'
include { zfin } from './databases/zfin'
include { zwd } from './databases/zwd'

workflow parse_databases {
  emit: data
  main:
    Channel.empty()
    | mix(
      five_s_rrnadb(),
      ena(),
      ensembl(),
      ensembl_fungi(),
      ensembl_metazoa(),
      ensembl_plants(),
      ensembl_protists(),
      flybase(),
      genecards_suite(),
      gtrnadb(),
      intact(),
      lncbase(),
      lncbook(),
      lncipedia(),
      mirbase(),
      mirgenedb(),
      pdbe(),
      pombase(),
      quickgo(),
      refseq(),
      rfam(),
      rgd(),
      sgd(),
      silva(),
      snodb(),
      snorna_database(),
      tarbase(),
      zfin(),
      zwd(),
    ) \
    | flatten \
    | set { data }
}
