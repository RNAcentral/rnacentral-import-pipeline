process modomics {
  when: { params.databases.modomics.run }

  output:
  path('*.csv')

  """
  if [[ "${params.databases.modomics.remote}" = /* ]]; then
    cp "${params.databases.modomics.remote}" modomics.json
  else
    wget -O modomics.json ${params.databases.modomics.remote}
  fi
  rnac modomics parse modomics.json .
  """
}
