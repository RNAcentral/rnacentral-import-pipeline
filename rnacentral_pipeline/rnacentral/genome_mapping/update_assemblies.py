import logging
import re
import typing as ty

import psycopg2

from rnacentral_pipeline.rnacentral.genome_mapping.urls import FtpHost, ftp

LOGGER = logging.getLogger(__name__)

# Divisions we check — skip bacteria (directory is huge)
CHECKED_DIVISIONS = {
    FtpHost.ensembl,
    FtpHost.ensembl_metazoa,
    FtpHost.ensembl_plants,
    FtpHost.ensembl_protists,
    FtpHost.ensembl_fungi,
}


def _fasta_base_path(host: FtpHost) -> str:
    if host is FtpHost.ensembl:
        return "current_fasta"
    groups = {
        FtpHost.ensembl_metazoa: "metazoa",
        FtpHost.ensembl_plants: "plants",
        FtpHost.ensembl_protists: "protists",
        FtpHost.ensembl_fungi: "fungi",
    }
    return f"current/{groups[host]}/fasta"


def _fetch_species_listing(host: FtpHost) -> ty.Set[str]:
    with ftp(host.host) as conn:
        conn.cwd(f"pub/{_fasta_base_path(host)}")
        return set(conn.nlst())


def _fetch_assembly_id(host: FtpHost, species: str) -> ty.Optional[str]:
    path = f"pub/{_fasta_base_path(host)}/{species}/dna"
    with ftp(host.host) as conn:
        try:
            conn.cwd(path)
            files = conn.nlst()
        except Exception:
            return None
    upper = species[0].upper() + species[1:]
    pattern = re.compile(rf"{re.escape(upper)}\.(.+?)\.dna\.(toplevel|primary_assembly).*\.fa\.gz")
    for f in files:
        m = pattern.match(f)
        if m:
            return m.group(1)
    return None


def _species_prefix(ensembl_url: str) -> str:
    for sep in ("_gca", "_gcf"):
        if sep in ensembl_url:
            return ensembl_url.split(sep)[0]
    return ensembl_url


def update(db_url: str) -> None:
    available: ty.Dict[FtpHost, ty.Set[str]] = {}
    for host in CHECKED_DIVISIONS:
        LOGGER.info("Fetching species listing for %s", host.name)
        try:
            available[host] = _fetch_species_listing(host)
        except Exception as e:
            LOGGER.warning("Could not fetch listing for %s: %s", host.name, e)
            available[host] = set()

    all_species: ty.Dict[str, FtpHost] = {
        s: host for host, species_set in available.items() for s in species_set
    }

    conn = psycopg2.connect(db_url)
    try:
        with conn.cursor() as cur:
            cur.execute("SELECT ensembl_url, assembly_id, division FROM ensembl_assembly")
            rows = cur.fetchall()

        updates = []
        for ensembl_url, assembly_id, division in rows:
            if ensembl_url in all_species:
                continue

            try:
                host = FtpHost.from_string(division)
            except ValueError:
                continue
            if host not in CHECKED_DIVISIONS:
                continue

            prefix = _species_prefix(ensembl_url)
            search_pool = available.get(host, set())
            matches = [s for s in search_pool if s.startswith(prefix + "_") or s == prefix]

            if len(matches) == 1:
                new_url = matches[0]
                new_assembly_id = _fetch_assembly_id(host, new_url)
                if new_assembly_id:
                    LOGGER.info(
                        "Updating %s -> %s (assembly %s -> %s)",
                        ensembl_url, new_url, assembly_id, new_assembly_id,
                    )
                    updates.append((new_url, new_assembly_id, ensembl_url, assembly_id))
                else:
                    LOGGER.warning(
                        "Found updated URL %s for %s but could not determine assembly_id",
                        new_url, ensembl_url,
                    )
            elif len(matches) == 0:
                LOGGER.info("No Ensembl entry for %s (genuinely missing)", ensembl_url)
            else:
                LOGGER.warning(
                    "Multiple candidates for %s: %s — skipping", ensembl_url, matches
                )

        if updates:
            with conn.cursor() as cur:
                for new_url, new_assembly_id, old_url, old_assembly_id in updates:
                    cur.execute(
                        "SELECT 1 FROM ensembl_assembly WHERE assembly_id = %s",
                        (new_assembly_id,),
                    )
                    new_exists = cur.fetchone()

                    cur.execute(
                        "SELECT 1 FROM rnc_sequence_regions WHERE assembly_id = %s LIMIT 1",
                        (old_assembly_id,),
                    )
                    has_mappings = cur.fetchone()

                    if new_exists or has_mappings:
                        LOGGER.info(
                            "Skipping %s — %s",
                            old_url,
                            "new assembly already in DB" if new_exists else "has existing mappings",
                        )
                    else:
                        cur.execute(
                            "UPDATE ensembl_assembly SET ensembl_url = %s, assembly_id = %s"
                            " WHERE ensembl_url = %s",
                            (new_url, new_assembly_id, old_url),
                        )
            conn.commit()
            LOGGER.info("Updated %d assemblies", len(updates))
        else:
            LOGGER.info("No assembly updates needed")
    finally:
        conn.close()
