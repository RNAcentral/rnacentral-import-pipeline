process fetch_and_process {
  when { params.databases.intact.run }

  output:
  path('*.csv')

  """
  pushd fetched
  wget -O intact.zip $params.databases.intact.remote
  gunzip intact.zip
  popd

  head -1 fetched/intact.txt > intact.txt
  grep -i rnacnetral fetched/intact.txt > intact.txt

  rnac external intact intact.txt .
  """
}

workflow intact {
  emit: data
  main:
    fetch_and_process | set { data }
}
