# -*- coding: utf-8 -*-

"""
Copyright [2009-present] EMBL-European Bioinformatics Institute
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
from contextlib import contextmanager
from ftplib import FTP
from typing import Dict, List, IO

LOGGER = logging.getLogger(__name__)
FTP_SERVER = "ftp.ebi.ac.uk"


class NoFastaFile(Exception):
    """
    Could not find any file to be used by IGV
    """


@contextmanager
def ftp(host):
    conn = FTP(host)
    conn.login()

    yield conn

    try:
        conn.quit()
    except Exception as err:
        LOGGER.info("Failed to close FTP connection")
        LOGGER.exception(err)


def check_url(file_list: List[str], name: str, path: str, assembly_id: str) -> Dict:
    result = {}
    fasta_file = [file for file in file_list if file == name + ".fa"]
    fasta_index = [file for file in file_list if file == name + ".fa.fai"]
    ensembl = [file for file in file_list if file == name + ".ensembl.gff3.gz"]
    ensembl_index = [file for file in file_list if
                     file == name + ".ensembl.gff3.gz.tbi" or file == name + ".ensembl.gff3.gz.csi"]
    rnacentral = [file for file in file_list if file == name + ".rnacentral.gff3.gz"]
    rnacentral_index = [file for file in file_list if
                        file == name + ".rnacentral.gff3.gz.tbi" or file == name + ".rnacentral.gff3.gz.csi"]

    if (ensembl or rnacentral) and fasta_file and fasta_index:
        result["name"] = assembly_id
        result["fastaFile"] = "https://" + FTP_SERVER + "/" + path + "/" + fasta_file[0]
        result["fastaIndex"] = "https://" + FTP_SERVER + "/" + path + "/" + fasta_index[0]

        if ensembl and ensembl_index:
            result["ensemblTrack"] = "https://" + FTP_SERVER + "/" + path + "/" + ensembl[0]
            result["ensemblIndex"] = "https://" + FTP_SERVER + "/" + path + "/" + ensembl_index[0]

        if rnacentral and rnacentral_index:
            result["rnacentralTrack"] = "https://" + FTP_SERVER + "/" + path + "/" + rnacentral[0]
            result["rnacentralIndex"] = "https://" + FTP_SERVER + "/" + path + "/" + rnacentral_index[0]

    return result


def create_json(species: str, assembly_id: str, output: IO[str]) -> None:
    with ftp(FTP_SERVER) as conn:
        path = "pub/databases/RNAcentral/.genome-browser"
        conn.cwd(path)
        file_list = conn.nlst()
        name = f"{species}.{assembly_id}"
        result = check_url(file_list, name, path, assembly_id)

        if result:
            json.dump(result, output)
        else:
            raise NoFastaFile(f"{species}.{assembly_id}")
