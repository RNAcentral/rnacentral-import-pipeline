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

import pyppeteer
from pyppeteer import launch
import asyncio

import requests

from rnacentral_pipeline.databases.plncdb import parser
from rnacentral_pipeline.writers import entry_writer


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
        writer.write(entries)


@cli.command("fetch")
@click.argument("remote", type=str)
@click.argument(
    "destination",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def fetch_data(remote, destination):
    ## fight the javascript app to get the right ftp urls for everything
    urls_dict = asyncio.get_event_loop().run_until_complete(render_urls(remote))

    ## Create directories, then download data into them
    for dir_name in urls_dict.keys():
        target_path = Path(destination) / dir_name
        target_path.mkdir(exist_ok=True)
        for url in urls_dict[dir_name]:
            download_file(url, target_path)


async def render_urls(remote):
    browser = await launch(headless=True)
    page = await browser.newPage()

    await page.goto(remote)
    await asyncio.sleep(5) ## This is needed because the site takes a little while to load, but doesn't tell us

    possible_links = await page.querySelectorAll(".span-name")
    link_names = await page.querySelectorAllEval('.span-name', '(nodes => nodes.map(n => n.innerText))')

    for link, name in zip(possible_links, link_names):
        if "For_RNAcentral" in name:
            await link.click()
            break
    await asyncio.sleep(5) ## This is needed because the site takes a little while to load, but doesn't tell us

    urls = {}
    species_names = set()
    link_names = await page.querySelectorAllEval('.span-name', '(nodes => nodes.map(n => n.innerText))')
    for name in link_names:
        short_name = ".".join(name.split('.')[0:2])
        species_names.add(short_name)
        if short_name in urls.keys():
            urls[short_name].append(f"{remote.replace('/plncdb', '')}/For_RNAcentral/{name.strip()}")
        else:
            urls[short_name] = [f"{remote.replace('/plncdb', '')}/For_RNAcentral/{name.strip()}"]

    # info_file_urls = len(species_names)*[remote.replace("/plncdb", "")]
    info_file_urls = {sn: remote.replace("/plncdb", "") for sn in species_names}

    await page.goto(remote, waitUntil='networkidle2')
    await _click_link_by_name(page, "PLncDB_v2.0")

    link_names = await page.querySelectorAllEval('.span-name', '(nodes => nodes.map(n => n.innerText))')
    # info_file_urls = [a +"PLncDB_v2.0/"+ b.strip().replace(" ", "") for a,b in zip(info_file_urls, link_names)]
    for link_name in link_names:
        info_file_urls[long_2_short_name(link_name)] += ( "PLncDB_v2.0/" + link_name.strip().replace(" ", "") )

    for s_name in link_names:
        await _click_link_by_name(page, s_name)
        f_name = await _find_link_with_name(page, "lncRNA_info")
        if f_name is not None:
            info_file_urls[long_2_short_name(s_name)] += f_name.strip()
        else:
            info_file_urls[[long_2_short_name(s_name)]] = None
        # await _click_link_by_name(page, "lncRNA_info")
        await page.goto(remote, waitUntil='networkidle2')
        await _click_link_by_name(page, "PLncDB_v2.0")

        break

    ## Merge the info urls into the data urls
    all_urls = merge_dicts(urls, info_file_urls)


    return all_urls


async def _click_link_by_name(page, target_name):
    possible_links = await page.querySelectorAll(".span-name")
    link_names = await page.querySelectorAllEval('.span-name', '(nodes => nodes.map(n => n.innerText))')

    for link, name in zip(possible_links, link_names):
        if target_name in name:
            await link.click()
            break
    await asyncio.sleep(2)

async def _find_link_with_name(page, target_name):
    possible_links = await page.querySelectorAll(".span-name")
    link_names = await page.querySelectorAllEval('.span-name', '(nodes => nodes.map(n => n.innerText))')

    for link, name in zip(possible_links, link_names):
        if target_name in name:
            return name

def long_2_short_name(long_name):
    return long_name[0:5]

def merge_dicts(dol1, dos):
  keys = set(dol1).union(dos)
  no = []
  return dict((k, dol1.get(k, no) + [dos.get(k, None)] ) for k in keys)

def download_file(url, destination=Path('.')):
    local_filename = url.split('/')[-1]
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(destination / local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                f.write(chunk)
    return local_filename
