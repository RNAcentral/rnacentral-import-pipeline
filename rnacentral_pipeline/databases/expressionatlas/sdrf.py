import logging
import os
from pathlib import Path
from typing import List, Union

import polars as pl


def parse_condensed_sdrf(path: Union[str, Path]) -> pl.DataFrame:
    """
    Parse a condensed SDRF file and return a polars DataFrame.

    A condensed SDRF file has 7 columns, but the last is often not delimited correctly,
    making it tricky to read with the polars default CSV reader.
    Therefore, manually parse the file into 6 series objects (one column seems to
    always be null) and construct a dataframe from them.

    Args:
        path: Path to the SDRF file

    Returns:
        A polars DataFrame containing the parsed SDRF data
    """
    logging.info(f"Loading sdrf data from {path}")

    # Read the file
    with open(path, "r") as file:
        content = file.read()

    # Parse lines and split by tabs
    part_parsed: List[List[str]] = [
        line.strip().split("\t") for line in content.splitlines()
    ]

    if len(part_parsed) == 0:
        raise ValueError(f"SDRF file {path} contained no data")
    # Initialize empty lists for each column
    exp_name_data = []
    assay_name_data = []
    feat_class_data = []
    feat_type_data = []
    feat_value_data = []
    ontology_data = []

    # Get the maximum number of columns
    max_columns = max(len(line) for line in part_parsed)

    # Parse according to the number of columns
    if max_columns == 7:
        for line in part_parsed:
            exp_name_data.append(line[0])
            assay_name_data.append(line[2])  # remember line[1] will be empty
            feat_class_data.append(line[3])
            feat_type_data.append(line[4])
            feat_value_data.append(line[5])
            if len(line) == 7:
                ontology_data.append(line[6])
            else:
                ontology_data.append(None)
    else:
        logging.warning(
            f"Unusual sdrf parsing with {max_columns} columns, not 7 for experiment {part_parsed[0][0]}"
        )
        for line in part_parsed:
            exp_name_data.append(line[0])
            assay_name_data.append(line[1])
            feat_class_data.append(line[2])
            feat_type_data.append(line[3])
            feat_value_data.append(line[4])

            if len(line) == 6:
                ontology_data.append(line[5])
            else:
                ontology_data.append(None)

    # Create the DataFrame with all series
    df = pl.DataFrame(
        {
            "exp_name": exp_name_data,
            "assay_name": assay_name_data,
            "feat_class": feat_class_data,
            "feat_type": feat_type_data,
            "feat_value": feat_value_data,
            "ontology": ontology_data,
        }
    )

    return df
