process fetch_and_process {
  tag { "$name" }
  input:
  tuple val(name), path(data), val(column_name)

  output:
  path('*.csv')

  """
  rnac lookup genecards $column_name $data urs-info.pickle
  rnac external $name $data urs-info.pickle .
  """
}
