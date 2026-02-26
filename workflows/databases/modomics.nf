process modomics {
  when: { params.databases.modomics.run }

  output:
  path('*.csv')

  """
  wget -O modomics.json ${params.databases.modomics.remote}
  rnac modomics parse modomics.json .
  """
}
