process pdbe {
  when { params.databases.pdb.run }

  output:
  path('*.csv')

  """
  rnac pdb generate .
  """
}
