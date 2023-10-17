"""
Mgnify's data is not quite complete, so we load it up and make some modifications.
These include:
- Figuring out a taxid for each entry
- Adding a type based on the source model

"""
import logging
import re
from functools import lru_cache

from psycopg2 import connect

from rnacentral_pipeline.databases.helpers import phylogeny as phy

LOGGER = logging.getLogger(__name__)


@lru_cache
def get_inferred_species_taxid(inferred_lineage):
    species = re.findall(";s__(.*)", inferred_lineage)[0]
    species = re.sub(r"_\w+", "", species)
    if len(species) == 0:
        ## There is no species in the inferred phylogeny, so try next layer up
        return None
    taxid = None

    try:
        taxid = phy.taxid(species)
    except:
        LOGGER.warning(f"{species} doesn't have a taxid, trying something else")
        species_name = species.split()[0]
        env_sample_name = f"uncultured {species_name} sp."
        try:
            taxid = phy.taxid(env_sample_name)
        except:
            LOGGER.warning(
                f"Uncultured species {env_sample_name} doesn't have a taxid, trying something else"
            )

    return taxid


@lru_cache
def get_inferred_genus_taxid(inferred_lineage):
    genus = re.findall(";g__(.*);", inferred_lineage)[0]
    if len(genus) == 0:
        ## There is no genus in the inferred phylogeny, so try next layer up
        return None

    genus = genus.split("_")[0]
    taxid = None

    try:
        taxid = phy.taxid(f"uncultured {genus} sp.")
    except:
        LOGGER.warning(f"{genus} doesn't have a taxid, trying something else")

    return taxid


@lru_cache
def get_inferred_family_taxid(inferred_lineage):
    family = re.findall(";f__(.*);g", inferred_lineage)[0]
    if len(family) == 0:
        ## There is no family in the inferred phylogeny, so try next layer up
        return None
    taxid = None

    try:
        taxid = phy.taxid(f"{family} bacterium")
    except:
        LOGGER.warning(f"{family} doesn't have a taxid, trying something else")

    return taxid


@lru_cache
def get_inferred_order_taxid(inferred_lineage):
    order = re.findall(";o__(.*);f", inferred_lineage)[0]
    if len(order) == 0:
        ## There is no order in the inferred phylogeny, so try next layer up
        return None
    taxid = None

    try:
        taxid = phy.taxid(f"{order}")
    except:
        LOGGER.warning(f"{order} doesn't have a taxid, trying something else")

    return taxid


def get_so_type(conn_str):
    conn = connect(conn_str)
    cur = conn.cursor()
    cur.execute(
        "SELECT rfam_model_id, so_rna_type FROM rfam_models ORDER BY rfam_model_id ASC;"
    )
    lookup = cur.fetchall()
    cur.close()
    conn.close()

    so_type = {a: b for a, b in lookup}

    return so_type


def prepare_mgnify_data(data, conn_str):

    ## Get Rfam model ID to SO type lookup
    so_type = get_so_type(conn_str)

    ## Define fallback taxids of the general metagenome of the environment
    ## These are used if we can't do any better
    fallback = {
        "zebrafish fecal genome catalogue": 1331678,  # zebrafish metagenome - more accurate then generic fish fecal?
        "human gut genome catalogue": 408170,  # human gut metagenome
        "human oral genome catalogue": 447426,  # human oral metagenome
        "marine genome catalogue": 2994539,  # human skin metagenome
        "pig gut genome catalogue": 1510822,  # pig gut metagenome
        "non model fish gut genome catalogue": 1602388,  # non model fish gut metagenome
        "honeybee gut genome catalogue": 1202446,  # insect gut metagenome
        "cow rumen genome catalogue": 506599,  # bovine gut metagenome
        "chicken gut genome catalogue": 506600,  # mouse gut metagenome
    }

    prepared_data = []
    for entry in data["data"]:

        taxid = get_inferred_species_taxid(entry["inferredPhylogeny"])

        if taxid is None:
            taxid = get_inferred_genus_taxid(entry["inferredPhylogeny"])

        if taxid is None:
            taxid = get_inferred_family_taxid(entry["inferredPhylogeny"])

        if taxid is None:
            taxid = get_inferred_order_taxid(entry["inferredPhylogeny"])

        if taxid is None:
            LOGGER.warning("falling back to generic metagenome taxid")
            taxid = fallback[entry["additionalAnnotations"]["catalog_name"]]

        entry["taxonId"] = f"NCBITaxon:{taxid}"

        ## Add the type of the sequence
        entry["soTermId"] = so_type[entry["sourceModel"].split(":")[1]]

        prepared_data.append(entry)

    raw = {"data": prepared_data, "metaData": data["metaData"]}

    return raw
