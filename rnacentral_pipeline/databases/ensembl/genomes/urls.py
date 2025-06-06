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
import logging
import re
import tempfile
import typing as ty
from contextlib import contextmanager
from ftplib import FTP

from rnacentral_pipeline.databases.ensembl.data import Division, FtpInfo

LOGGER = logging.getLogger(__name__)


def list_releases(ftp: FTP) -> ty.List[str]:
    return [f for f in ftp.nlst() if f.startswith("release-")]


def latest_release(releases: ty.List[str], ftp: FTP) -> str:
    ## Parse the readme for the current release to avoid getting a half baked release
    readme_lines = []
    ftp.retrlines("RETR current_README", readme_lines.append)
    cur_readme = "\n".join(readme_lines)
    pattern = r"Ensembl Release (\d+) Databases."
    release = re.search(pattern, cur_readme).group(1)
    print(f"Ensembl release {release}")
    return f"release-{release}"


@contextmanager
def species_info(ftp: FTP, division: Division, release: str):
    info_path = f"{release}/species_metadata_{division.division_name}.json"
    print(info_path)
    with tempfile.NamedTemporaryFile() as tmp:
        ftp.retrbinary(f"RETR {info_path}", tmp.write)
        tmp.flush()
        tmp.seek(0)
        yield tmp


def generate_paths(
    ftp: FTP, division: Division, base: str, release: str, handle
) -> ty.Iterable[FtpInfo]:
    _, release_id = release.split("-", 1)
    data = json.load(handle)
    for entry in data:
        info = entry["organism"]
        name = info["name"]

        # So sometimes Ensembl users lower case names for the assemblies but
        # not always, so we try both and generate a name for the existing one.
        assemblies = [
            entry["assembly"]["assembly_default"],
            entry["assembly"]["assembly_default"].lower(),
        ]
        for assembly in assemblies:
            organism_name = f"{info['url_name']}.{assembly}.{release_id}"

            # This detects, and skips things that are part of a collection. I'm not
            # sure what that means right now and those seem to be things that have
            # other genomes that aren't nested in a collection.
            if not any(db["dbname"].startswith(name) for db in entry["databases"]):
                continue

            gff_path = f"{base}/{release}/gff3/{name}/{organism_name}.gff3.gz"
            data_files = f"{base}/{release}/embl/{name}/{organism_name}.*.dat.gz"

            yield FtpInfo(
                division=division,
                species=name,
                data_files=data_files,
                gff_file=gff_path,
            )
            break
        else:
            LOGGER.warn("No files found for %s", info)


def urls_for(division: Division, host: str) -> ty.Iterable[FtpInfo]:
    with FTP(host) as ftp:
        ftp.login()
        print("LOGIN")
        ftp.cwd(f"pub/{division.name}/")
        releases = list_releases(ftp)
        latest = latest_release(releases, ftp)
        with species_info(ftp, division, latest) as info:
            url_base = f"ftp://{host}/pub/{division.name}"
            yield from generate_paths(ftp, division, url_base, latest, info)
