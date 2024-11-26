# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

import numpy as np
import polars as pl
from datasets import load_dataset
from sklearn.model_selection import train_test_split


def train_test_val(
    input_data, train_path, test_path, val_path, test_frac, val_frac, seed, hub_repo
):
    all_examples = pl.read_parquet(input_data)

    ## Set up one random state for both splits
    random_state = np.random.RandomState(seed)

    ## Double randomised splitting
    train, test_val = train_test_split(
        all_examples, test_size=(test_frac + val_frac), random_state=random_state
    )
    test, val = train_test_split(test_val, test_size=0.5, random_state=random_state)

    ## Now write the output
    train.write_parquet(train_path)
    val.write_parquet(val_path)
    test.write_parquet(test_path)

    if hub_repo is not None:
        dataset = load_dataset(
            "parquet",
            data_files={"train": train_path, "validation": val_path, "test": test_path},
        )
        dataset.push_to_hub(hub_repo)
