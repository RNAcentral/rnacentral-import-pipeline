process regulondb {
  memory '4GB'

  output:
  path('*.{csv,parquet}')

  when: params.databases.regulondb?.run

  script:
  """
  curl -sfG --data-urlencode 'query={ tus: getDataOfFile(fileName: "TUSet") { content } promoters: getDataOfFile(fileName: "PromoterSet") { content } terminators: getDataOfFile(fileName: "TerminatorSet") { content } }' ${params.databases.regulondb.graphql} > regulondb.json
  wget -O genome.embl '${params.databases.regulondb.genome}'
  rnac regulondb parse regulondb.json genome.embl .
  """
}
