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

  def configFile = new File("secrets.json");
  def config = new groovy.json.JsonSlurper().parseFile(configFile, 'UTF-8');

  def post = new URL(config.SLACK_WEBHOOK).openConnection();
  post.setRequestMethod("POST")
  post.setDoOutput(true);
  post.setRequestProperty("Content-Type", "application/json");

  def  payload = "{\"text\" : \"$msg\" }"


  post.getOutputStream().write(payload.getBytes("UTF-8"));
  def postRC = post.getResponseCode();
  if (postRC != 200) {
    println("Something went wrong calling slack webhook!");
    println(post.getInputStream().getText());
  }

}
