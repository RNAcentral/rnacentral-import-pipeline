process slack_message {
  input:
  val(message)

  when: params.get('notify', false)

  script:
  params.notify ? """
  rnac notify step "Import Workflow" "$message"
  """ : "true"
}


process slack_file {
  input:
  path(message)

  when: params.get('notify', false)

  script:
  params.notify ? """
  rnac notify file "$message"
  """ : "true"
}


// A groovy function for use in closures - uses groovy's own URL class to make the request
def slack_closure(msg) {
  if (!params.notify) {
    return
  }

  def configFile = new File("secrets.json")
  if (!configFile.exists()) {
    log.warn "Slack notify is enabled but secrets.json was not found; skipping notification"
    return
  }

  try {
    def config = new groovy.json.JsonSlurper().parseFile(configFile, 'UTF-8')
    def webhook = config.SLACK_WEBHOOK
    if (!webhook) {
      log.warn "Slack notify is enabled but SLACK_WEBHOOK is missing from secrets.json; skipping notification"
      return
    }

    def post = new URL(webhook).openConnection()
    post.setRequestMethod("POST")
    post.setDoOutput(true)
    post.setRequestProperty("Content-Type", "application/json")

    def payload = groovy.json.JsonOutput.toJson([text: msg])
    post.getOutputStream().write(payload.getBytes("UTF-8"))

    def postRC = post.getResponseCode()
    if (postRC != 200) {
      log.warn "Slack webhook returned ${postRC}: ${post.errorStream?.text}"
    }
  } catch (Exception e) {
    log.warn "Could not send Slack notification: ${e}"
  }
}
