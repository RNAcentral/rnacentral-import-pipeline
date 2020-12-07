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

import json
import tempfile
from ftplib import FTP
import typing as ty
from contextlib import contextmanager

from rnacentral_pipeline.databases.ensembl.data import Division


def list_releases(ftp: FTP) -> ty.List[str]:
    return [f for f in ftp.nlst() if f.startswith('release-')]


def latest_release(releases: ty.List[str]) -> str:
    return max(releases, key=lambda r: int(r.split('-')[1]))


@contextmanager
def species_info(ftp: FTP, division: Division, release: str):
    info_path = f'{release}/species_metadata_{division.division_name}.json'
    with tempfile.NamedTemporaryFile() as tmp:
        ftp.retrbinary(f"RETR {info_path}", tmp.write)
        tmp.flush()
        tmp.seek(0)
        yield tmp


def generate_paths(base: str, release: str, handle) -> ty.Iterable[ty.Tuple[str, str, str]]:
    _, release_id = release.split('-', 1)
    data = json.load(handle)
    data_types = ['gff3', 'embl']
    for entry in data:
        info = entry['organism']
        name = info['name']
        url_name = info['url_name']
        assembly = entry['assembly']['assembly_default']
        organism_name = f"{url_name}.{assembly}.{release_id}"
        gff_path = f"{base}/{release}/gff3/{name}/{organism_name}.gff3.gz"
        data_files = f"{base}/{release}/embl/{name}/{organism_name}.*.dat.gz"
        yield (name, data_files, gff_path)


def urls_for(division: Division, host: str):
    with FTP(host) as ftp:
        ftp.login()
        ftp.cwd(f'pub/{division.name}/')
        releases = list_releases(ftp)
        latest = latest_release(releases)
        with species_info(ftp, division, latest) as info:
            yield from generate_paths(f'ftp://{host}/pub/{division.name}', latest, info)
