process send_completion_report {
  executor 'local'

  script:
  """
  rnac notify report
  """
}

workflow {
  send_completion_report()
}
