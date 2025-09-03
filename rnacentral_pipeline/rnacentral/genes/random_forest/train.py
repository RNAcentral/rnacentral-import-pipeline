# -*- coding: utf-8 -*-

"""
Copyright [2009-2025] EMBL-European Bioinformatics Institute
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


import polars as pl
import numpy as np
import requests
from gensim.models import Word2Vec
from sklearn.model_selection import train_test_split, KFold
from sklearn.ensemble import RandomForestClassifier
from skl2onnx import to_onnx
import onnxruntime as rt
import time
import pathlib
from datasets import load_dataset
from rnacentral_pipeline.rnacentral.genes.random_forest.preprocessing import compare_transcripts



url = "https://rest.ensembl.org/lookup/id/{0}?content-type=application/json"


def fetch_data(gene):
    r = requests.get(url.format(gene.split(".")[0]))
    data = r.json()

    return {
        "gene_start": data["start"],
        "gene_stop": data["end"],
        "gene_strand": data["strand"],
        "e_chromosome": data["seq_region_name"],
        "assembly": data["assembly_name"],
    }


def convert_csv_2_parquet(input_csv: str, output_parquet: str):
    """
    Convert CSV file to Parquet format using Polars library.

    Args:
        input_csv (str): Path to the input CSV file.
        output_parquet (str): Path to the output Parquet file.
    """
    df = pl.read_csv(input_csv)
    df.write_parquet(output_parquet)

def locate_genes(df, nearby_distance):
    """
    Loop on each gene, find any genes contained within another gene
    The search genes within 1kb from start or end
    These become the 'nearby' genes to provide the negative examples
    """
    candidates_df = []
    for gene in df.iter_rows(named=True):
        candidates = df.filter(
            pl.col("gene_start").is_between(
                gene["gene_start"] - nearby_distance,
                gene["gene_stop"] + nearby_distance,
            )
            | pl.col("gene_stop").is_between(
                gene["gene_start"] - nearby_distance,
                gene["gene_stop"] + nearby_distance,
            )
            & (pl.col("assembly") == gene["assembly"])
        )
        candidates_df.append(
            {
                "gene": gene["gene"],
                "candidates": [gene["gene"]],
                "labels": [1],
                "contained": [True],
            }
        )
        for c in candidates.iter_rows(named=True):
            if c["gene"] == gene["gene"]:
                continue
            candidates_df[-1]["candidates"].append(c["gene"])
            candidates_df[-1]["labels"].append(0)
            contained = (
                c["gene_start"] > gene["gene_start"]
                and c["gene_stop"] < gene["gene_stop"]
            )
            candidates_df[-1]["contained"].append(contained)

    candidates_df = pl.DataFrame(candidates_df)
    return candidates_df

def fetch_training_data(data: str, nearby_distance: int) -> pl.DataFrame:
    """
    
    """
    transcripts = pl.read_parquet(data)
    transcripts = transcripts.group_by(
        ["assembly_id", "gene", "region_name", "strand"], maintain_order=True
    ).agg(
        [
            pl.col("urs_taxid").first(),
            pl.col("chromosome").first(),
            pl.col("region_start").first(),
            pl.col("region_stop").first(),
            pl.col("exon_start"),
            pl.col("exon_stop"),
            pl.col("so_rna_type").first().alias("so_type"),
        ]
    )

    genes = transcripts.select(pl.col("gene").unique()).collect()
    if genes.height == 0:
        print("No genes found with selected parameters, try something else")
        exit()
    
    ## We don't store gene level coordinate data, so get it from ensembl
    genes = genes.with_columns(
        res=pl.col("gene").map_elements(fetch_data)
    ).unnest("res")
    genes = genes.filter(pl.col("gene_start").is_not_null())

    candidates = genes.group_by(["assembly", "e_chromosome", "gene_strand"]).map_groups(
        lambda x: locate_genes(x, nearby_distance)
    )

    return candidates


def build_training_features(candidates: pl.DataFrame, transcripts: pl.DataFrame, so_model_path:str) -> pl.DataFrame:
    so_model = Word2Vec.load(so_model_path)
    ## Features can now be built for these candidate genes
    features = []
    for row in candidates.iter_rows(named=True):
        ## The extract genes step has already guaranteed that genes and candidates are from the same assembly
        cand_genes = row["candidates"]
        cand_labels = row["labels"]
        row_transcripts = transcripts.filter(pl.col("gene") == row["gene"])

        for cand, label in zip(cand_genes, cand_labels):
            t_data = transcripts.filter(pl.col("gene") == cand)

            ## Mark up all transcripts in this gene as being in the gene
            features.extend(
                compare_transcripts(t_data, row_transcripts, so_model, label=label)
            )

    all_comparisons_features = pl.DataFrame(features).unique("comparison")

    return all_comparisons_features


def split_datasets(input_data: str, train_path: str, val_path:str, test_path:str, test_frac:float=0.2, val_frac:float=0.2, seed:int=1337, hub_repo:str=None):
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

def convert_model(model, output_path, check_data):
        ## Use the SKL model to predict some classes on our check data
        start_skl = time.time()
        classes = model.predict(check_data)
        end_skl = time.time()
        ## Convert to onnx
        model_onnx = to_onnx(model, check_data[:1])

        sess = rt.InferenceSession(output_path, providers=["CPUExecutionProvider"])
        input_name = sess.get_inputs()[0].name
        label_name = sess.get_outputs()[0].name
        start_onx = time.time()
        pred_onx = sess.run([label_name], {input_name: check_data})[0]
        end_onx = time.time()

        print(f"sklearn took {end_skl - start_skl} seconds")
        print(f"onnx took {end_onx - start_onx} seconds")
        if np.allclose(classes, pred_onx):
            with open(output_path, "wb") as f:
                f.write(model_onnx.SerializeToString())
            print(f"Model verification passed, saved to {output_path}")
        else:
            print("ONNX model predictions differ to sklearn, not saving")

def train(training_data, output_folder, folds, exclude, seed, basename, n_estimators):
    ## Load the data
    if not (training_data.endswith("pq") or training_data.endswith("parquet")):
        print("Convert your file to parquet format please")

    train = pl.read_parquet(training_data)

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

        model_output_path = pathlib.Path(output_folder) / f"{basename}{i}.onnx"
        names = train.select(pl.exclude(excludes)).columns
        feat_importances.append(
            {n: f for n, f in zip(names, model.feature_importances_)}
        )
        ## This will convert the model to ONNX and save it in the right place
        convert_model(model, model_output_path, X[test_idx])

    importances = pl.DataFrame(feat_importances)
    importances.write_parquet(pathlib.Path(output_folder) / "importances.parquet")



