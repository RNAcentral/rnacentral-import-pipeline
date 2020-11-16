process fetch_and_process {
  when: { params.databases["5srrnadb"].run }

  output:
  path('*.csv')

  """
  scp ${params.databases["5srrnadb"].remote} 5s.json
  rnac external 5srrnadb 5s.json .
  """
}

workflow five_s_rrnadb {
  emit: fetch_and_process.out
  main: 
    fetch_and_process()
}
