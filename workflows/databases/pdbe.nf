process fetch_and_process {
  when { params.databases.pdb.run }

  output:
  path('*.csv')

  """
  rnac fetch pdb data pdb.json
  rnac fetch pdb extra pdb-extra.json
  rnac external pdb pdb.json pdb-extra.json .
  """
}

workflow mirbase {
  emit: data
  main:
    fetch_and_process | set { data }
}
