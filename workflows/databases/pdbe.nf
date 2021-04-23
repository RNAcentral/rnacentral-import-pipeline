process pdbe {
  when { params.databases.pdb.run }

  output:
  path('*.csv')

  """
  rnac pdb fetch data pdb.json
  rnac pdb fetch extra pdb.json pdb-extra.pickle
  rnac pdb parse pdb.json pdb-extra.pickle .
  """
}
