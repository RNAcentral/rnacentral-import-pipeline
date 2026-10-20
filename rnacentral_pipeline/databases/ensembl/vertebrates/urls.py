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
import typing as ty
from contextlib import contextmanager
from ftplib import FTP

from rnacentral_pipeline.databases.ensembl.data import Division, FtpInfo


def latest_release(ftp: FTP) -> str:
    ## VERSION rather than a release-* listing, to avoid picking up a half baked
    ## release. current_README, which we used to parse, is gone from the FTP site.
    lines: ty.List[str] = []
    ftp.retrlines("RETR VERSION", lines.append)
    release = "".join(lines).strip()
    if not release.isdigit():
        raise ValueError(f"Could not determine latest Ensembl release, got {release!r}")
    return f"release-{release}"


@contextmanager
def species_info(ftp: FTP, release: str):
    info_path = f"{release}/species_metadata_EnsemblVertebrates.json"
    # dir="." because singularity --contain gives us a tiny in-memory /tmp
    with tempfile.NamedTemporaryFile(dir=".") as tmp:
        ftp.retrbinary(f"RETR {info_path}", tmp.write)
        tmp.flush()
        tmp.seek(0)
        yield tmp


def generate_paths(base: str, release: str, handle) -> ty.Iterable[FtpInfo]:
    _, release_id = release.split("-", 1)
    data = json.load(handle)
    for entry in data:
        info = entry["organism"]
        name = info["name"]
        assembly = entry["assembly"]["assembly_default"]
        organism_name = f"{name[0].upper() + name[1:]}.{assembly}.{release_id}"
        gff_path = f"{base}/{release}/gff3/{name}/{organism_name}.gff3.gz"
        data_files = f"{base}/{release}/embl/{name}/{organism_name}.*.dat.gz"
        yield FtpInfo(
            division=Division.vertebrates,
            species=name,
            data_files=data_files,
            gff_file=gff_path,
        )


def urls_for(host: str) -> ty.Iterable[FtpInfo]:
    ftp = FTP(host)
    try:
        ftp.login()
        ftp.cwd("pub")
        latest = latest_release(ftp)
        with species_info(ftp, latest) as info:
            yield from generate_paths(f"ftp://{host}/pub", latest, info)
    finally:
        # close() not quit(); a failed transfer makes QUIT raise over the real error
        ftp.close()
