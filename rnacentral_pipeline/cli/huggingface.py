# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.db import cursor
from rnacentral_pipeline.rnacentral.huggingface_export.metadata import (
    get_active_sequence_count,
)
from rnacentral_pipeline.rnacentral.huggingface_export.upload import (
    create_dataset,
    create_or_find_collection,
    upload_dataset,
)


@click.group("huggingface")
def cli():
    """
    This is a set of commands for publishing data to huggingface
    """


@cli.command("create-collection")
@click.argument("release", type=str)
def create_collection(release):
    """
    Creates the collection for this release. All datasets will be grouped under release-collections
    If the collection exists, this will return the slug of the existing collection, so beware.

    Collections are by default private, we will need to switch to public once the release is completed
    """
    slug_str = create_or_find_collection(release)

    with open("collection_slug", "w") as slug_fh:
        slug_fh.write(slug_str)


@cli.command("create-dataset")
@click.argument("release", type=str)
def create_hf_dataset(release: str):
    """
    Creates a dataset repo under RNAcentral, i.e. RNAcentral/{name}.

    Will be private in the first instance and you must run the publish function (rnac huggingface publish)
    to update all repo visibilities in the given collection
    """

    repo_url = create_dataset(release)
    if repo_url is None:
        """
        This means the huggingface hub found a repo with the name already in existence.
        So we need to give up and try again/delete the existing dataset
        """
        return 1

    with open("dataset_url", "w") as dataset_fh:
        dataset_fh.write(repo_url)


@cli.command("create-readme")
@click.argument("release")
@click.argument("update_specs")
@click.option("--db-url", envvar="PGDATABASE")
def create_hf_readme(release, update_specs, db_url):
    """
    Create a template README with some metadata from the database
    """
    with cursor(db_url) as cur:
        seq_count = get_active_sequence_count(cur)

    dataset_size = "10M<n<100M" if seq_count < 100_000_000 else "100M<n<1B"
    template = """---
language: rna
license: cc0-1.0
tags:
    - Biology
    - ncRNA
size_categories:
    - {size_cat}
task_categories:
    - text-generation
    - fill-mask
task_ids:
    - language-modeling
    - masked-language-modeling
pretty_name: RNAcentral
---

# RNAcentral Release {release_no}
This is an official HuggingFace version of the RNAcentral Release {release_no} available at https://rnacentral.org/

Included in this dataset are:
- ncRNA sequences for {seq_count:,} sequences from our member databases

# Using the data
All RNAcentral data is provided with a CC-0 license, so you are free to do as you with with it.

If you do make use of our data, please consider citing our paper: https://doi.org/10.1093/nar/gkaf1329

## Detailed database versions
As part of our releases, we import updates from our expert member databases. Here is the detailed breakdown of
database versions for all those that have been updated this time around.

{database_update_specs}
"""

    database_specs = Path(update_specs).read_text()
    readme_content = template.format(
        size_cat=dataset_size,
        release_no=release,
        seq_count=seq_count,
        database_update_specs=database_specs,
    )

    readme_fh = Path("README.md")

    readme_fh.write_text(readme_content)


@cli.command("upload-data")
@click.argument("release", type=str)
@click.argument("dataset_path", type=str)
def upload_hf_dataset(release: str, dataset_path: str):
    """
    Creates a dataset repo under RNAcentral, i.e. RNAcentral/{name}.

    Will be private in the first instance and you must run the publish function (rnac huggingface publish)
    to update all repo visibilities in the given collection
    """
    upload_dataset(release, dataset_path)


@cli.command("publish")
@click.argument("release")
def publish__hf_release(release):
    """
    Gets the release collection slug and uses it to find all the datasets in it and make them public
    then makes the collection public.
    """
    collection_slug = create_collection(release)
