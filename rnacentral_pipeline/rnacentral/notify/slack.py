"""
Send a notification to slack.

NB: The webhook should be configured in the nextflow profile

"""

import os
import requests

from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError
import psycopg2

REPORT_QUERY = """
SELECT display_name, count(taxid) FROM xref
JOIN rnc_database db ON xref.dbid = db.id
WHERE xref.deleted = 'N'
AND EXTRACT (DAY FROM (CURRENT_TIMESTAMP - xref.timestamp)) < 7
GROUP BY display_name
ORDER BY display_name
"""

def send_notification(title, message, plain=False):
    """
    Send a notification to the configured slack webhook.
    """
    SLACK_WEBHOOK = os.getenv('SLACK_WEBHOOK')
    if SLACK_WEBHOOK is None:
        try:
            from rnacentral_pipeline.secrets import SLACK_WEBHOOK
        except:
            raise SystemExit("SLACK_WEBHOOK environment variable not defined, and couldn't find a secrets file")

    if plain:
        slack_json = {
            "text" : title + ':  ' + message
        }
    else:
        slack_json = {
            "text" : title,
            "blocks" : [
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": message
                    },
                },
            ]
        }
    try:
        response = requests.post(SLACK_WEBHOOK,
                    json=slack_json,
                    headers={'Content-Type':'application/json'}
                    )
        response.raise_for_status()
    except Exception as request_exception:
        raise SystemExit from request_exception

def pipeline_report():
    """
    Generates a nicely formatted report of the number of sequences imported from
    each DB. This uses the slack_sdk, rather than a webhook, and uses the
    blockkit to format the message nicely.

    TODO: What else should go in this? Maybe parsing the log file to get the
    run duration? 
    """
    db_url = os.getenv('PGDATABASE')
    client_token = os.getenv('SLACK_CLIENT_TOKEN')
    channel = os.getenv('SLACK_CHANNEL')

    client = WebClient(token=client_token)

    lock_text_template = "New sequences from *{0}* {1:,}"


    summary_blocks = [{
    "type": "header",
    "text": {
      "type": "plain_text",
      "text": "Workflow Completion report"
      }
    },
    {
    "type": "divider"
    }]
    running_total = 0
    with psycopg2.connect(db_url) as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.DictCursor) as cur:
            cur.execute(REPORT_QUERY)
            res = cur.fetchall()
            for r in res:
                running_total += r[1]
                summary_blocks.append(
                {
                  "type": "section",
                  "text": {
                    "type": "mrkdwn",
                    "text": block_text_template.format(r[0].ljust(30), r[1])
                  }
                }
                )
                summary_blocks.append({
    			"type": "divider"
        		}
                )
            summary_blocks.append(
            {
              "type": "section",
              "text": {
                "type": "mrkdwn",
                "text": f"Total sequences imported: *{running_total:,}*"
              }
            })

    try:
        response = client.chat_postMessage(
            channel=channel,
            text="Workflow completion report",
            blocks=summary_blocks)

        print(response)
    except SlackApiError as e:
        assert e.response["error"]
