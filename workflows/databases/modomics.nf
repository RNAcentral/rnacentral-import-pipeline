process modomics {
  when: { params.databases.modomics.run }
  containerOptions "${params.common_container} --bind /nfs/production/agb/rnacentral/provided-data:/nfs/production/agb/rnacentral/provided-data"

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
