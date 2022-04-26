"""
Send a notification to slack.

NB: The webhook should be configured in the nextflow profile

"""

import os
import requests


def send_notification(title, message, plain=False):
    """
    Send a notification to the configured slack webhook.
    """
    SLACK_WEBHOOK = os.getenv('SLACK_WEBHOOK')
    if SLACK_WEBHOOK is None:
        raise SystemExit("SLACK_WEBHOOK environment variable not defined")

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
