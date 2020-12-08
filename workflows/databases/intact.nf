process intact {
  when { params.databases.intact.run }

  output:
  path('*.csv')

  """
  mkdir fetched
  pushd fetched
  wget -O intact.zip $params.databases.intact.remote
  gunzip intact.zip
  popd

  head -1 fetched/intact.txt > intact.txt
  grep -i rnacentral fetched/intact.txt > intact.txt

  rnac external intact intact.txt .
  """
}
