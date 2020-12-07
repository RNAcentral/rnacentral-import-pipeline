process pdbe {
  when { params.databases.pdb.run }

  output:
  path('*.csv')

  """
  rnac pdb fetch data pdb.json
  rnac pdb fetch extra pdb-extra.json
  rnac pdb parse pdb.json pdb-extra.json .
  """
}
