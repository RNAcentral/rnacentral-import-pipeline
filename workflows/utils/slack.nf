process slack_message {

  input:
  val(message)

  """
  rnac notify step "Import Workflow" "$message"
  """

}


process slack_file {

  input:
  path(message)

  """
  rnac notify file "$message"
  """

}


import groovy.json.JsonSlurper

// A groovy function for use in closures - uses groovy's own URL class to make the request
def slack_closure(msg) {
  def configFile = new File("secrets.json");
  def config = new JsonSlurper().parseFile(configFile, 'UTF-8');

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
