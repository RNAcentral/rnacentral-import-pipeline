process send_completion_report {
  executor 'local'

  """
  rnac notify report
  """
}

workflow {
  send_completion_report()
}
