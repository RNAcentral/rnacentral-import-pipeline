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

import urllib
from pathlib import Path
import typing as ty

from furl import furl
import requests
from bs4 import BeautifulSoup


def base_url(url: furl) -> furl:
    base = furl(url)
    base.path.segments = base.path.segments[:-1]
    return base


def extract_urls(base: furl, document: str) -> ty.List[furl]:
    soup = BeautifulSoup(document)
    urls = []
    links = soup.find("table").find_all("a")
    for link in links:
        href = link.get("href")
        if href.endswith("json.gz"):
            urls.append(base / href)
    return urls


def find_urls(url: furl):
    response = requests.get(url.url)
    response.raise_for_status()
    return extract_urls(base_url(url), response.text)
