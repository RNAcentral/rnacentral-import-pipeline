import gzip
import re
from pathlib import Path

import polars as pl

urs_regex = re.compile(r"(URS[0-9A-Z]+_\d+)")


def get_all_written_ids_gff(ftp_path: Path) -> pl.DataFrame:
    """
    Load all the GFF/bed files in the FTP export and construct a big ol list
    of them to compare against the database
    """
    if isinstance(ftp_path, str):
        ftp_path = Path(ftp_path)

    all_written_ids = set()
    for gff_path in ftp_path.glob("*.gff3.gz"):
        print(gff_path)
        with gzip.open(gff_path, "r") as handle:
            written_ids = set()
            this_org_count = 0
            content = handle.read().decode("utf-8")
            lines = content.splitlines()
            for line in lines:
                if line.startswith("#"):
                    continue
                match = urs_regex.search(line)
                if match:
                    this_org_count += 1
                    written_ids.add(match.group(1))
        print(
            f"Filename: {gff_path.name}: Total coordinates: {this_org_count} Total unique ids: {len(written_ids)}"
        )
        all_written_ids.update(written_ids)
    print("Total written ids:", len(all_written_ids))
    return pl.DataFrame({"urs_taxid": list(all_written_ids)})


def get_all_written_ids_bed(ftp_path: Path) -> pl.DataFrame:
    """
    Load all the bed files in the FTP export and construct a big ol list
    of them to compare against the database
    """
    if isinstance(ftp_path, str):
        ftp_path = Path(ftp_path)

    all_written_ids = set()
    for gff_path in ftp_path.glob("*.bed.gz"):
        print(gff_path)
        with gzip.open(gff_path, "r") as handle:
            written_ids = set()
            this_org_count = 0
            content = handle.read().decode("utf-8")
            lines = content.splitlines()
            for line in lines:
                if line.startswith("#"):
                    continue
                match = urs_regex.search(line)
                if match:
                    this_org_count += 1
                    written_ids.add(match.group(1))
        print(
            f"Filename: {gff_path.name}: Total coordinates: {this_org_count} Total unique ids: {len(written_ids)}"
        )
        all_written_ids.update(written_ids)
    print("Total written ids:", len(all_written_ids))
    return pl.DataFrame({"urs_taxid": list(all_written_ids)})


def get_all_mapped_ids(db: str) -> pl.DataFrame:
    """
    Get all the URS_taxids from the sequence regions table and make a unique list

    """
    query = "SELECT id as urs_taxid FROM rnc_rna_precomputed WHERE (is_active AND has_coordinates)"

    all_db_ids = pl.read_database_uri(query, db).unique()

    return all_db_ids


def compare_gff(ftp_path: Path, db: str):
    """
    Compare the ids in the GFF files to the ids in the database
    """
    all_written_ids = get_all_written_ids_gff(ftp_path)
    exit()
    all_db_ids = get_all_mapped_ids(db)

    ids_missing_in_files = all_db_ids.join(all_written_ids, on="urs_taxid", how="anti")

    print(f"Found {len(all_written_ids)} ids in the files")
    print(f"Found {len(all_db_ids)} ids in the database")
    print(
        f"Found {len(ids_missing_in_files)} ids in the database that are not in the files"
    )


def compare_bed(ftp_path: Path, db: str):
    """
    Compare the ids in the GFF files to the ids in the database
    """
    all_written_ids = get_all_written_ids_bed(ftp_path)
    exit()
    all_db_ids = get_all_mapped_ids(db)

    ids_missing_in_files = all_db_ids.join(all_written_ids, on="urs_taxid", how="anti")

    print(f"Found {len(all_written_ids)} ids in the files")
    print(f"Found {len(all_db_ids)} ids in the database")
    print(
        f"Found {len(ids_missing_in_files)} ids in the database that are not in the files"
    )
