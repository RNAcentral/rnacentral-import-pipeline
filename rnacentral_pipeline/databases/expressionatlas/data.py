import pathlib

import polars as pl
import polars.selectors as cs


def load_filter_differential(
    file_path: pathlib.Path, lookup: pl.DataFrame, exp_name: str
) -> pl.LazyFrame:
    """
    Loads a differential experiment dataframe and filters it for
    abs(log2foldchange) > 1
    """
    significant_data = pl.scan_csv(file_path, separator="\t").filter(
        pl.any_horizontal(cs.contains("log2foldchange").abs().gt(1))
    )
    significant_data = significant_data.with_columns(experiment_name=pl.lit(exp_name))

    if not "GeneID" in significant_data.columns:
        significant_data = significant_data.rename({"Gene ID": "GeneID"})

    significant_data = significant_data.join(
        lookup, left_on="GeneID", right_on="external_id", how="inner"
    )
    return significant_data.select(["urs_taxid", "taxid", "experiment_name"])


def load_filter_baseline(
    file_path: pathlib.Path, lookup: pl.DataFrame, exp_name: str
) -> pl.LazyFrame:
    """
    Loads a baseline experiment dataframe and filters on
    """

    baseline_data = pl.scan_csv(file_path, separator="\t")
    baseline_data = baseline_data.with_columns(experiment_name=pl.lit(exp_name))
    baseline_data = baseline_data.with_columns(
        cs.starts_with("g")
        .str.split(",")
        .list.eval(pl.element().cast(pl.Float64))
        .list.median()
    ).filter(pl.any_horizontal(cs.starts_with("g") > 0.0))

    if not "GeneID" in baseline_data.columns:
        baseline_data = baseline_data.rename({"Gene ID": "GeneID"})

    baseline_data = baseline_data.join(
        lookup, left_on="GeneID", right_on="external_id", how="inner"
    )

    return baseline_data.select(["urs_taxid", "taxid", "experiment_name"])


def load_sdrf_file(file_path: pathlib.Path):
    """
    Read an SDRF file and extract the characteristics

    SDRF files are annoying because they don't have a stable number of columns so we probably have to parse them line-by line

    The columns are (probably)
    0 - The experiment name, same on every line
    1 - Usually empty (there's a double tab)
    2 - Assay Name
    3 - feature class
    4 - feature type
    5 - feature name
    6 - ontology
    """
    ontology_terms = {}
    taxid = None
    with open(file_path, "r") as sdrf:
        for line in sdrf:
            parts = line.split("\t")
            exp_name = parts[0]
            array_name = parts[2]
            if not array_name in ontology_terms.keys():
                ontology_terms[array_name] = []

            if parts[4] == "organism":
                taxid = int(parts[6].split("_")[-1])
            elif len(parts) > 5 and parts[-1].startswith("http"):
                term_id = parts[-1].split("/")[-1].replace("_", ":").strip()
                ontology_terms[array_name].append(term_id)

    ## For now, we will just flatten all the ontology terms in the whole experiment into a set
    ## This way we link ontology terms to an experiment rather than an array

    ontology_set = set()
    for array_ontologies in ontology_terms.values():
        ontology_set.update(array_ontologies)

    df = pl.DataFrame({"experiment_name": exp_name, "taxid": taxid})
    df = df.with_columns(ontology_terms=list(ontology_set))
    return df
