# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from pathlib import Path

import click
from furl import furl

import requests

from rnacentral_pipeline.databases.plncdb import parser, get_urls
from rnacentral_pipeline.writers import entry_writer
from rnacentral_pipeline.rnacentral.notify.slack import send_notification

@click.group("plncdb")
def cli():
    """
    A group of commands dealing with PLncDB data.
    """
    pass


@cli.command("parse")
@click.argument("data", type=click.Path(dir_okay=True, readable=True, file_okay=False))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def parse(data, output):
    entries = parser.parse(Path(data))
    with entry_writer(Path(output)) as writer:
        try:
            writer.write(entries)
        except ValueError as e:
            print(e)


@cli.command("fetch-data")
@click.argument("urls", type=click.Path(writable=False, file_okay=True, dir_okay=False))
@click.argument("destination", type=click.Path(writable=True, file_okay=False, dir_okay=True), default='.')
def fetch_data(urls, destination):
    url_dict = {}
    with open(urls, 'r') as url_file:
        for url_line in url_file:
            url_dict[url_line.split(',')[0]] = url_line.split(',')[1:]


    for dir_name in url_dict.keys():
        print(f"Getting data for {dir_name}")

        send_notification("PLncDB Download", f"Getting data for {dir_name}")

        target_path = Path(destination) / dir_name
        target_path.mkdir(exist_ok=True, parents=True)
        for url in url_dict[dir_name]:
            download_file(url, target_path)
        print(f"All data for {dir_name} is downloaded")


@cli.command("get-urls")
@click.argument("remote", type=str)
@click.argument(
    "destination",
    default="urls.txt",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
)
def get_urls_cli(remote, destination):
    get_urls(remote, destination)


def download_file(url, destination=Path('.')):
    local_filename = url.split('/')[-1]
    if (destination / local_filename).exists():
        return local_filename
    # NOTE the stream=True parameter below
    with requests.get(url.strip(), stream=True) as r:
        r.raise_for_status()
        with open(destination / local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                f.write(chunk)
    return local_filename
