# -*- coding: utf-8 -*-

"""
Copyright [2009-current] EMBL-European Bioinformatics Institute
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
import operator as op
import typing as ty

import polars as pl

LOGGER = logging.getLogger(__name__)

from pathlib import Path

from . import configuration, data, helpers

# from rnacentral_pipeline.databases import data


def as_expression(mapping):

    pass


def parse(data_path, lookup_path):
    """
    Process the jsonlines output from Rust into entries.

    The jsonlines output is already grouped by geneID and urs taxid so this
    should give us the transcript level linkage we're after without any further
    processing.
    """
    ## Step 0. Get a lookup from gene to URS taxid
    lookup = pl.scan_csv(lookup_path).select(["urs_taxid", "taxid", "external_id"])
    # This thing splits on | then expodes the resulting list and filters out empty IDs
    lookup.with_columns(pl.col("external_id").str.strip_chars().str.split("|")).explode(
        "external_id"
    ).with_columns(pl.col("external_id").str.split(",")).explode(
        "external_id"
    ).with_columns(
        external_id=pl.when(pl.col("external_id").str.strip_chars() == "None")
        .then(pl.lit(""))
        .otherwise(pl.col("external_id"))
    ).filter(
        pl.col("external_id").str.lengths() > 0
    ).collect()

    ## Step 1. Build a lookup between experiment name and configuration xml
    config_files = Path(data_path).glob("*configuration.xml")
    config_lookup = configuration.build_exp_config_lookup(config_files)

    ## Step 2. Split experiments based on type. Either differential or baseline
    differential_experiments = []
    baseline_experiments = []
    sdrf_data = []
    for exp_name, exp_config in config_lookup.items():
        if "differential" in exp_config["experimentType"]:
            analysis = exp_config["analytics"].get("array_design", None)
            print(exp_name, f"|{analysis}|", f"|{exp_config['experimentType']}|")
            if analysis is not None:
                data_filepath = Path(data_path) / f"{exp_name}_{analysis}-analytics.tsv"
            else:
                data_filepath = Path(data_path) / f"{exp_name}-analytics.tsv"

            if data_filepath.exists():
                # Now parse the differential experiment to a dataframe and return it here
                differential_experiments.append(
                    data.load_filter_differential(data_filepath, lookup, exp_name)
                )
            else:
                LOGGER.warning(
                    f"Could not find file for {data_filepath.name}, skipping the experiment"
                )
                continue
        else:
            data_filepath = Path(data_path) / f"{exp_name}-tpms.tsv"
            if data_filepath.exists():
                # Now parse the baseline experiment to a dataframe and return it here
                baseline_experiments.append(
                    data.load_filter_baseline(data_filepath, lookup, exp_name)
                )
            else:
                LOGGER.warning(
                    f"Could not find file for {data_filepath.name}, skipping the experiment"
                )
                continue

        ## Step 3. Now we have two big lists of URS_taxid - experiment name lookups
        ## We want to build the lookup between experiment and taxid, ontology terms

        sdrf_path = Path(data_path) / f"{exp_name}.condensed-sdrf.tsv"
        if sdrf_path.exists():
            exp_sdrf = data.load_sdrf_file(sdrf_path)
            sdrf_data.append(exp_sdrf)
        else:
            LOGGER.warning(
                f"Could not find an sdrf file for {exp_name}, skipping the experiment"
            )
            continue
    ## Step 4. Now we have loaded everything, we can stack all the dataframes and augment with sdrf
    ## info.

    differential_experiments.extend(baseline_experiments)
    # Concatenate all lazy frames with the expression data
    all_data = pl.concat(differential_experiments)
    # Concatenate all the experiment sdrf files
    all_sdrf = pl.concat(sdrf_data)

    # augment the urs -> experiment lookup with sdrf data
    all_data = all_data.join(
        all_sdrf.lazy(), on=["experiment_name", "taxid"], how="inner"
    )

    ## Step 5. Now we have the lookup from urs to experiment with all the extra data, we can join against the
    ## rest of the lookup data using only URS_taxid to get the right sequence and other info

    full_lookup = pl.scan_csv(lookup_path)

    all_data = all_data.join(full_lookup, on="urs_taxid", how="inner")

    return all_data
