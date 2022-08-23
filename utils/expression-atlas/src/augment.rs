use crate::configuration::*;
/// This module will combine the config and dataframes to create a df we can extract all the
/// necessary information from.
use anyhow::Result;
use multimap::MultiMap;
use polars::frame::DataFrame;
use polars::prelude::*;

use log::{info, warn};

/// This function will 'augment' the experiment dataframe with information from the config
/// What that means is adding location and factor information as series to the right of the
/// dataframe. I do this based on the assay group id (the gN in the expression data) which I
/// map to the assay name. The assay name is then mapped to the location and factor data in the
/// sdrf file, where I try to separate factors and locations. This also allows grabbing the
/// taxid for the experiment
pub fn augment_differential_df(
    df: &mut DataFrame,
    config: &Config,
    sdrf: &DataFrame,
) -> Result<DataFrame> {
    // preprocess the config into a hash map for later convenience
    info!("Parsing config into lookup MultiMap");
    let mut contrast_lookup: MultiMap<String, (String, String)> = MultiMap::new();
    for analysis in &config.analytics {
        for cont in &analysis.contrasts.as_ref().unwrap().contrast {
            let test_group = cont.test_group.clone();

            let mut assay_names = analysis.assay_groups.assay_group.clone();
            assay_names.retain(|ag| ag.id == test_group);

            contrast_lookup.insert(
                test_group,
                (format!("{}.log2foldchange", cont.id), format!("{}.p-value", cont.id)),
            );
        }
    }

    // Set up some dataframes for the things we want
    let mut taxid_df = DataFrame::default();
    let mut localisation_df = DataFrame::default();
    let mut disease_df = DataFrame::default();
    let mut cell_type_df = DataFrame::default();

    let mut df_result = DataFrame::default();

    for analysis in &config.analytics {
        for ass_group in &analysis.assay_groups.assay_group {
            // Build the series of assay names for matching later
            let mut assay_names =
                Utf8ChunkedBuilder::new("assay_names", ass_group.assays.len(), 128);
            for ass_nm in &ass_group.assays {
                assay_names.append_value(ass_nm);
            }
            // Select out just the bits for this assay group from sdrf
            let assay_df = sdrf
                .clone()
                .lazy()
                .filter(col("assay_name").is_in(lit(assay_names.finish().into_series())))
                .with_column(lit(ass_group.id.as_str()).alias("group_id"));

            let mut df_inter = df.clone();

            df_inter = df_inter.lazy().with_column(lit(NULL).alias("group_id")).collect()?;

            // If there is localisation data, this will select it
            if localisation_df.height() == 0 {
                localisation_df = get_localisation_data(assay_df.clone())?;
            } else {
                localisation_df.vstack_mut(&get_localisation_data(assay_df.clone())?)?;
            }

            // If there is disease data, this will select it
            if disease_df.height() == 0 {
                disease_df = get_disease_data(assay_df.clone())?;
            } else {
                disease_df.vstack_mut(&get_disease_data(assay_df.clone())?)?;
            }

            // If there is cell type data, this will select it
            if cell_type_df.height() == 0 {
                cell_type_df = get_cell_type_data(assay_df.clone())?;
            } else {
                cell_type_df.vstack_mut(&get_cell_type_data(assay_df.clone())?)?;
            }

            // Get the taxonomy ontology reference. T
            if taxid_df.height() == 0 {
                taxid_df = get_taxonomy_data(assay_df.clone())?;
            } else {
                taxid_df.vstack_mut(&get_taxonomy_data(assay_df.clone())?)?;
            }

            if contrast_lookup.get(&ass_group.id) == None {
                continue;
            }

            let check_cols_vec = contrast_lookup.get_vec(&ass_group.id).unwrap();

            for check_cols in check_cols_vec.iter() {
                // Check columns are actually present
                if !&df_inter.get_column_names().contains(&check_cols.0.as_str())
                    || !&df_inter.get_column_names().contains(&check_cols.1.as_str())
                {
                    warn!("Couldn't find either {} or {} in the DataFrame columns, skipping this contrast", check_cols.0, check_cols.1);
                    continue;
                }

                df_inter = df_inter
                    .lazy()
                    .with_column(
                        when(col(&check_cols.0).and(col(&check_cols.1)))
                            .then(lit(ass_group.id.as_ref()))
                            .otherwise(col("group_id"))
                            .alias("group_id"),
                    )
                    .filter(col(&check_cols.0).and(col(&check_cols.1)))
                    .collect()?;

                if df_result.height() == 0 {
                    df_result = df_inter.clone();
                } else {
                    df_result.vstack_mut(&df_inter)?;
                }
            }
        }
    } // closes loop on analyses

    // Use a df to join on selectively
    df_result =
        join_augmentations(&df_result, &taxid_df, &localisation_df, &disease_df, &cell_type_df)?;

    Ok(df_result)
}

/// Augmentation function for baseline experiments
pub fn augment_baseline_df(
    df: &mut DataFrame,
    config: &Config,
    sdrf: &DataFrame,
) -> Result<DataFrame> {
    // Set up some dataframes for the things we want
    let mut taxid_df = DataFrame::default();
    let mut localisation_df = DataFrame::default();
    let mut disease_df = DataFrame::default();
    let mut cell_type_df = DataFrame::default();

    let mut df_inter = df.clone();
    df_inter = df_inter.lazy().with_column(lit(NULL).alias("group_id")).collect()?;

    let mut df_result = DataFrame::default();

    for analysis in &config.analytics {
        for assay_group in &analysis.assay_groups.assay_group {
            let assay_names = assay_group.assays.clone();
            let ass_group = assay_group.id.clone();

            let mut assay_names_series =
                Utf8ChunkedBuilder::new("assay_names", assay_names.len(), 128);
            for ass_nm in &assay_names {
                assay_names_series.append_value(ass_nm);
            }

            let assay_df = sdrf
                .clone()
                .lazy()
                .filter(col("assay_name").is_in(lit(assay_names_series.finish().into_series())))
                .with_column(lit(ass_group.as_str()).alias("group_id"));

            if localisation_df.height() == 0 {
                localisation_df = get_localisation_data(assay_df.clone())?;
            } else {
                localisation_df.vstack_mut(&get_localisation_data(assay_df.clone())?)?;
            }

            if disease_df.height() == 0 {
                disease_df = get_disease_data(assay_df.clone())?;
            } else {
                disease_df.vstack_mut(&get_disease_data(assay_df.clone())?)?;
            }

            if cell_type_df.height() == 0 {
                cell_type_df = get_cell_type_data(assay_df.clone())?;
            } else {
                cell_type_df.vstack_mut(&get_cell_type_data(assay_df.clone())?)?;
            }

            if taxid_df.height() == 0 {
                taxid_df = get_taxonomy_data(assay_df.clone())?;
            } else {
                taxid_df.vstack_mut(&get_taxonomy_data(assay_df.clone())?)?;
            }

            // Filter the experiment data according to how expresison atlas does it
            df_inter = df_inter
                .lazy()
                .with_column(
                    when(col(ass_group.as_str()))
                        .then(lit(ass_group.as_ref()))
                        .otherwise(col("group_id"))
                        .alias("group_id"),
                )
                .filter(col(ass_group.as_str()))
                .collect()?;

            if df_result.height() == 0 {
                df_result = df_inter.clone();
            } else {
                df_result.vstack_mut(&df_inter)?;
            }
        } // close loop on assay groups
    } //close loop on analyses

    df_result =
        join_augmentations(&df_result, &taxid_df, &localisation_df, &disease_df, &cell_type_df)?;

    Ok(df_result)
}

fn get_localisation_data(assay_df: LazyFrame) -> Result<DataFrame> {
    let localisation = assay_df
        .filter(
            ((col("feat_class").eq(lit("factor"))).or(col("feat_class").eq(lit("characteristic"))))
                .and(col("feat_type").eq(lit("organism part"))),
        )
        .select([col("group_id"), col("ontology").alias("location")])
        .first()
        .collect()?;

    Ok(localisation)
}

fn get_disease_data(assay_df: LazyFrame) -> Result<DataFrame> {
    let disease = assay_df
        .clone()
        .filter(
            ((col("feat_class").eq(lit("factor"))).or(col("feat_class").eq(lit("characteristic"))))
                .and(col("feat_type").eq(lit("disease"))),
        )
        .select([col("group_id"), col("ontology").alias("disease")])
        .first()
        .collect()?;

    Ok(disease)
}

fn get_cell_type_data(assay_df: LazyFrame) -> Result<DataFrame> {
    let cell_type = assay_df
        .clone()
        .filter(
            col("feat_class").eq(lit("characteristic")).and(col("feat_type").eq(lit("cell type"))),
        )
        .select([col("group_id"), col("ontology").alias("cell_type")])
        .first()
        .collect()?;

    Ok(cell_type)
}

fn get_taxonomy_data(assay_df: LazyFrame) -> Result<DataFrame> {
    let tax_data = assay_df
        .clone()
        // .lazy()
        .filter(col("feat_type").eq(lit("organism")))
        .select([col("group_id"), col("ontology").alias("taxonomy")])
        .first()
        .collect()?;

    Ok(tax_data)
}

fn join_augmentations(
    df_result_bare: &DataFrame,
    taxid_df: &DataFrame,
    localisation_df: &DataFrame,
    disease_df: &DataFrame,
    cell_type_df: &DataFrame,
) -> Result<DataFrame> {
    // Use a df to join on selectively
    let mut df_result = df_result_bare.clone();
    if taxid_df.height() > 0 {
        df_result = df_result.join(&taxid_df, ["group_id"], ["group_id"], JoinType::Inner, None)?;
    } else {
        df_result = df_result
            .lazy()
            .with_column(lit(NULL).cast(DataType::Utf8).alias("taxonomy"))
            .collect()?;
    }

    if localisation_df.height() > 0 {
        df_result =
            df_result.join(&localisation_df, ["group_id"], ["group_id"], JoinType::Inner, None)?;
    } else {
        df_result = df_result
            .lazy()
            .with_column(lit(NULL).cast(DataType::Utf8).alias("location"))
            .collect()?;
    }

    if disease_df.height() > 0 {
        df_result =
            df_result.join(&disease_df, ["group_id"], ["group_id"], JoinType::Inner, None)?;
    } else {
        df_result = df_result
            .lazy()
            .with_column(lit(NULL).cast(DataType::Utf8).alias("disease"))
            .collect()?;
    }

    if cell_type_df.height() > 0 {
        df_result =
            df_result.join(&cell_type_df, ["group_id"], ["group_id"], JoinType::Inner, None)?;
    } else {
        df_result = df_result
            .lazy()
            .with_column(lit(NULL).cast(DataType::Utf8).alias("cell_type"))
            .collect()?;
    }

    df_result = df_result.select([
        "GeneID",
        "Gene Name",
        "experiment",
        "taxonomy",
        "location",
        "disease",
        "cell_type",
    ])?;
    Ok(df_result)
}
