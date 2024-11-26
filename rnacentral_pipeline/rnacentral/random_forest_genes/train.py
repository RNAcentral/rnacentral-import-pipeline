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

import pickle as pkl
from pathlib import Path
from tempfile import mkdtemp

import polars as pl
import sklearn
import skops
from huggingface_hub import HfApi
from sklearn.ensemble import RandomForestClassifier
from skops import card, hub_utils
from skops import io as sio


def no_kfold(training_data, model_file, n_estimators, exclude, seed, hub_path):
    ## Load the data
    if not (training_data.endswith("pq") or training_data.endswith("parquet")):
        print("Convert your file to parquet format please")

    train = pl.read_parquet(training_data)

    ## Handle exclusion of the label and anything else we should exclude
    excludes = ["label"]
    exclude = list(exclude)
    excludes.extend(exclude)

    ## Extract as numpy arrays. Note - this won't currently work if the dataframe
    ## contains list columns
    X = train.select(pl.exclude(excludes)).to_numpy()
    y = train.select("label").to_numpy().flatten()

    ## Build the model
    model = RandomForestClassifier(n_estimators=n_estimators, random_state=seed)

    ## Train the model
    model.fit(X, y)

    ## Write the model
    if hub_path is None:
        with open(model_file, "wb") as model_weights:
            pkl.dump(model, model_weights)
    else:
        api = HfApi()
        temp_dir = Path(mkdtemp())
        hub_dir = Path(mkdtemp())
        sio.dump(model, temp_dir / model_file)
        model_card = card.Card(model=model)

        requirements = [
            f"scikit-learn=={sklearn.__version__}",
            f"polars=={pl.__version__}",
            f"skops=={skops.__version__}",
        ]

        hub_utils.init(
            model=temp_dir / model_file,
            requirements=requirements,
            dst=hub_dir,
            task="tabular-classification",
            data=train.to_pandas(),
        )
        metadata = card.metadata_from_config(hub_dir / "config.json")
        model_card.metadata = metadata
        model_card.save(hub_dir / "README.md")
        try:
            api.create_repo(hub_path, repo_type="model")
        except:
            pass
        api.upload_folder(folder_path=hub_dir, repo_id=hub_path, repo_type="model")


def with_kfold(
    training_data, output_folder, folds, exclude, seed, basename, n_estimators
):
    ## Load the data
    if not (training_data.endswith("pq") or training_data.endswith("parquet")):
        print("Convert your file to parquet format please")

    train = pl.read_parquet(training_data)
    print(train.head())

    ## Handle exclusion of the label and anything else we should exclude
    excludes = ["label"]
    exclude = list(exclude)
    excludes.extend(exclude)

    ## Set up one random state for both splits
    random_state = np.random.RandomState(seed)

    ## Extract as numpy arrays. Note - this won't currently work if the dataframe
    ## contains list columns
    X = train.select(pl.exclude(excludes)).to_numpy()
    y = train.select("label").to_numpy().flatten()

    kf = KFold(n_splits=folds, random_state=random_state, shuffle=True)

    feat_importances = []
    for i, (train_idx, test_idx) in enumerate(kf.split(X)):
        # print(f"Training model on fold {i}")
        model = RandomForestClassifier(
            n_estimators=n_estimators, oob_score=True, n_jobs=-1
        )
        model.fit(X[train_idx], y[train_idx])
        if not pathlib.Path(output_folder).exists():
            pathlib.Path(output_folder).mkdir(parents=True)

        model_output_path = pathlib.Path(output_folder) / f"{basename}{i}.pkl"
        names = train.select(pl.exclude(excludes)).columns
        feat_importances.append(
            {n: f for n, f in zip(names, model.feature_importances_)}
        )

        with open(model_output_path, "wb") as output_file:
            pkl.dump(model, output_file)

    importances = pl.DataFrame(feat_importances)
    print(importances)
    importances.write_parquet(pathlib.Path(output_folder) / "importances.pq")
