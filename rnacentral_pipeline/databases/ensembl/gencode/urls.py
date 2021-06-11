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

import re
from ftplib import FTP
import typing as ty


def human_key(release):
    match = re.match("^Gencode_human/release_(\d+)(\w*)$", release)
    if not match:
        raise ValueError(f"Could not parse {release}")
    return (int(match.group(1)), match.group(2))


def mouse_key(release):
    parts = release.split("_", 2)
    return int(parts[2][1:])


def latest_human(ftp: FTP) -> str:
    names = ftp.nlst("Gencode_human")
    names = [f for f in names if f.startswith("Gencode_human/release_")]
    best = max(names, key=human_key)
    release_id = "".join(str(p) for p in human_key(best))
    return f"{best}/gencode.v{release_id}.annotation.gff3.gz"


def latest_mouse(ftp: FTP) -> str:
    names = ftp.nlst("Gencode_mouse")
    names = [f for f in names if f.startswith("Gencode_mouse/release_M")]
    best = max(names, key=mouse_key)
    release_id = f"M{mouse_key(best)}"
    return f"{best}/gencode.v{release_id}.annotation.gff3.gz"


def urls_for(host: str) -> ty.Iterable[ty.Tuple[str]]:
    path = "pub/databases/gencode"
    with FTP(host) as ftp:
        ftp.login()
        ftp.cwd(path)
        human_release = latest_human(ftp)
        mouse_release = latest_mouse(ftp)
        yield (f"ftp://{host}/{path}/{human_release}",)
        yield (f"ftp://{host}/{path}/{mouse_release}",)
