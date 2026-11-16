"""
Search URS_taxids from RNAcentral that have GO annotations for a matching GO-CAM
"""

import logging
import time
import typing as ty

import httpx

from rnacentral_pipeline.databases.gocams.urls import GOCAMS_API

# import pdb
# pdb.set_trace()


def urs_taxid_2_curie(urs_taxid):
    if ":" in urs_taxid:
        logging.error(f"Must supply bare URS_taxid, not a cure. Got {urs_taxid}")
        raise ValueError(f"Must supply bare URS_taxid, not a cure. Got {urs_taxid}")
    else:
        curie = f"RNAcentral:{urs_taxid}"
    return curie


DEFAULT_RATE_LIMIT = 5.0  # requests per second


def find_gocam(
    urs_taxid_iter: ty.Iterable[str], rate_limit: float = DEFAULT_RATE_LIMIT
) -> ty.Iterable[tuple[str, str]]:
    """
    Search a list of URS_taxid for a valid GO-CAM and return the ID
    Must be URS_taxids, not curies - they are converted in the function
    Requests are spaced to at most `rate_limit` per second.

    returns a list of all searched URS_taxids and their goCAMs, one per line
    so if a URS_taxid is involved in multiple CAMs, it will have multiple lines
    """

    if rate_limit <= 0:
        raise ValueError(f"rate_limit must be positive, got {rate_limit}")
    min_interval = 1.0 / rate_limit
    last_request = None

    with httpx.Client(timeout=10) as client:
        for urs_taxid in urs_taxid_iter:
            curie = urs_taxid_2_curie(urs_taxid)
            if last_request is not None:
                wait = min_interval - (time.monotonic() - last_request)
                if wait > 0:
                    time.sleep(wait)
            last_request = time.monotonic()
            url = f"{GOCAMS_API}/gp/{curie}/models"
            try:
                # breakpoint()
                resp = client.get(url).raise_for_status()
                for m in resp.json():
                    yield urs_taxid, m["gocam"].split("/")[-1]

            except httpx.HTTPStatusError as e:
                if e.response.status_code == 404:
                    logging.debug(f"Failed to find a gocam for {url}")
                else:
                    raise e
