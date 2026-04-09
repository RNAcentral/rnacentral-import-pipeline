# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

import logging
from pathlib import Path

import huggingface_hub as hf
from huggingface_hub import HfApi
from huggingface_hub.errors import HfHubHTTPError

LOGGER = logging.getLogger(__name__)


def create_or_find_collection(release: str) -> str:
    """
    Create a release on hugingface and return the release slug.
    If the release already exists, will just get the slug and return it.
    """

    title = f"RNAcentral Release {release}"

    description = f"Datasets from RNAcentral's release {release}"

    collection = hf.create_collection(
        title=title,
        description=description,
        private=True,
        exists_ok=True,
    )

    return collection.slug


def create_dataset(release: str) -> str | None:
    """
    Creates a repo of the correct type for this release
    """

    api = HfApi()
    dataset_name = f"release-{release}"
    repo_name = f"RNAcentral/{dataset_name}"
    collection_slug = create_or_find_collection(release)

    try:

        repo_url = hf.create_repo(
            repo_id=repo_name, repo_type="dataset", private=True, exist_ok=False
        )

        api.add_collection_item(
            collection_slug=collection_slug, item_id=repo_name, item_type="dataset"
        )

    except HfHubHTTPError:
        LOGGER.error(f"Repo with name {repo_name} already exists, can't overwrite it")
        repo_url = None

    return repo_url


def upload_dataset(release, dataset_path):
    """
    Uploads a written parquet file to the huggingface hub
    """
    repo_name = f"release-{release}"
    repo_id = f"RNAcentral/{repo_name}"
    dataset_name = Path(dataset_path).name
    api = HfApi()
    api.upload_file(
        path_or_fileobj=dataset_path,
        path_in_repo=dataset_name,
        repo_id=repo_id,
        repo_type="dataset",
    )


def publish(collection_slug):
    api = HfApi()
    collection = api.get_collection(collection_slug)
    for item in collection.items:
        repo_id = item.item_id
        api.update_repo_settings(repo_id=repo_id, private=False)
    api.update_collection_metadata(collection_slug=collection_slug, private=False)
