#!/usr/bin/env python3
"""
Scan a single RNA ID against Europe PMC and save results to the litscan DB.

This is a synchronous, self-contained port of the seek_references() function
from rnacentral-references (consumer/views/submit_job.py) for use as a SLURM job.

Usage:
    litscan-scan-job.py --job-id <job_id>

Environment variables:
    PSYCOPG_CONN    PostgreSQL connection URI
    LITSCAN_MODEL   Path to the scikit-learn SVC pipeline .pkl file
"""
import argparse
import datetime
import gzip
import logging
import os
import re
from contextlib import ExitStack
from functools import partial
from pathlib import Path

import joblib
import nltk
import polars as pl
import psycopg2
from lxml import etree as ET
from polars.exceptions import NoDataError

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def get_previous_searches(conn, job_ids):
    """
    PostgreSQL/psycopg2 optimized version using the ANY() operator.
    """
    query = "SELECT job_id, pmcid FROM litscan_result WHERE job_id = ANY(%s)"

    df = pl.read_database(
        query=query, connection=conn, execute_options={"parameters": (job_ids,)}
    )
    return df.cast({"job_id": pl.String, "pmcid": pl.String})


# ---------------------------------------------------------------------------
# Core scan logic  (ported from seek_references() in submit_job.py)
# ---------------------------------------------------------------------------

bad_tags = {
    "counts",
    "table-wrap",
    "table",
    "fig-group",
    "fig",
    "supplementary-material",
}


def extract_clean_text(element):
    if element is None:
        return ""
    # .itertext() grabs all text fragments inside the element and its children
    raw_text_pieces = element.itertext()

    # We strip out weird XML whitespace/newlines and ignore empty chunks,
    # then join everything back together with a single clean space.
    clean_text = " ".join(piece.strip() for piece in raw_text_pieces if piece.strip())

    return clean_text


def get_section_text(element, ignore_tags):
    """Extracts text, but stops at nested <sec> tags to prevent duplication."""
    if not isinstance(element.tag, str):
        return ""

    text_pieces = []
    if element.text:
        text_pieces.append(element.text)

    for child in element:
        if not isinstance(child.tag, str):
            continue

        # We add 'sec' to the ignore list so we don't recurse into sub-sections.
        # We also ignore 'title' so the section header isn't duplicated in the body text.
        if child.tag not in ignore_tags and child.tag not in ["sec", "title"]:
            text_pieces.append(get_section_text(child, ignore_tags))

        # ALWAYS grab the tail, even if we ignored the tag itself!
        if child.tail:
            text_pieces.append(child.tail)

    return " ".join(piece.strip() for piece in text_pieces if piece and piece.strip())


def parse_body_to_dict(body_elem, ignore_tags):
    sections_dict = {}
    if body_elem is None:
        return sections_dict

    # Find every <sec> tag anywhere in the body
    for i, sec in enumerate(body_elem.findall(".//sec")):

        # 1. Grab the Title
        title_tag = sec.find("title")
        if title_tag is not None:
            # Use itertext() just in case the title has bold/italic tags inside it
            section_name = "".join(title_tag.itertext()).strip()
            t = section_name.lower()
            if re.match(r".*intro.+", t):
                key = "intro"
            elif re.match(r".*results", t):
                key = "results"
            elif re.match(r".*discussion", t):
                key = "discussion"
            elif re.match(r".*conclusion.*", t):
                key = "conclusion"
            elif re.match(r".*method.+", t):
                key = "method"
            else:
                key = "other"
        else:
            key = "other"

        # Prevent key collisions (e.g., if there are two "Results" sub-sections)
        if key in sections_dict:
            section_name = f"{key}_{i+1}"

        # 2. Grab the Text (excluding nested sections)
        section_text = get_section_text(sec, ignore_tags)

        # 3. Add to dictionary if it actually contains text
        if section_text:
            sections_dict[key] = section_text

    return sections_dict


def parse_floats_to_dict(floats_elem):
    floats_dict = {"figures": {}, "tables": {}}
    if floats_elem is not None:
        # Set up a dictionary with sub-categories

        # --- 1. Process Figures ---
        for fig in floats_elem.findall(".//fig"):
            # .findtext() safely grabs the text, or returns None if the tag doesn't exist
            fig_label = fig.findtext("label") or "Unknown_Figure"

            caption_tag = fig.find("caption")
            if caption_tag is not None:
                caption_text = get_section_text(caption_tag, ignore_tags=[])
                floats_dict["figures"][fig_label] = nltk.sent_tokenize(caption_text)

        # --- 2. Process Tables ---
        for table in floats_elem.findall(".//table-wrap"):
            table_label = table.findtext("label") or "Unknown_Table"

            caption_tag = table.find("caption")
            if caption_tag is not None:
                # We use the exact same get_section_text function here!
                caption_text = get_section_text(caption_tag, ignore_tags=[])
                floats_dict["tables"][table_label] = nltk.sent_tokenize(caption_text)

    return floats_dict


def extract_article_base(elem, rna_pipeline):
    """Extract the job-agnostic fields of an article. Expensive work (text
    extraction, sentence tokenization, classifier) happens here once per article
    so it can be reused across all jobs scanning the same pmcid."""
    title = extract_clean_text(elem.find(".//article-title"))
    abstract = extract_clean_text(elem.find(".//abstract"))
    body = parse_body_to_dict(elem.find(".//body"), ignore_tags=bad_tags)
    floats = parse_floats_to_dict(elem)

    abstract_all_sentences = nltk.sent_tokenize(abstract) if abstract else []
    body_section_sentences = {
        key: nltk.sent_tokenize(section_text) for key, section_text in body.items()
    }

    doi_el = elem.find('.//article-id[@pub-id-type="doi"]')
    doi = doi_el.text if doi_el is not None else ""
    pmid_el = elem.find('.//article-id[@pub-id-type="pmid"]')
    pmid = pmid_el.text if pmid_el is not None else ""

    year = 0
    for item in elem.findall("./front/article-meta/pub-date"):
        if {"epub", "ppub", "pub"}.intersection(item.attrib.values()):
            year_el = item.find("year")
            if year_el is not None and year_el.text:
                year = int(year_el.text)

    journal = ""
    journal_el = elem.find("./front/journal-meta/journal-title-group/journal-title")
    if journal_el is None:
        journal_el = elem.find("./front/journal-meta/journal-title")
    if journal_el is not None:
        journal = journal_el.text or ""

    authors = []
    for auth in elem.findall(".//contrib-group//name"):
        surname = auth.findtext("surname", default="").strip()
        given = auth.findtext("given-names", default="").strip()
        if surname and given:
            authors.append(f"{surname}, {given}")
        elif surname or given:
            authors.append(surname + given)
    author = "; ".join(authors)

    article_type = (
        elem.attrib.get("article-type", "").strip().replace("-", " ").capitalize()
    )

    probability, rna_related = classify_abstract(abstract, rna_pipeline)

    return {
        "title": title,
        "title_lower": title.lower(),
        "abstract": abstract,
        "body": body,
        "floats": floats,
        "doi": doi,
        "pmid": pmid,
        "year": year,
        "journal": journal,
        "author": author,
        "rna_related": rna_related,
        "probability": probability,
        "abstract_all_sentences": abstract_all_sentences,
        "body_section_sentences": body_section_sentences,
        "type": article_type,
    }


def _job_regex(job_id):
    return (
        r"(^|\s|\(|\u201c|\u2018|\u201d|\;)"
        + re.escape(job_id.lower())
        + r"($|[\s.,:;?\u2019\u201d\u201c\"/)])"
    )


def finalize_for_job(base, pmcid, job_id, cite_count):
    """Apply the per-job regex to cached article text and produce a row dict."""
    regex = _job_regex(job_id)

    abstract_sentences = [
        s for s in base["abstract_all_sentences"] if re.search(regex, s.lower())
    ]
    id_in_abstract = len(abstract_sentences) > 0
    id_in_title = bool(re.search(regex, base["title_lower"]))

    body_sentences = []
    locations = []
    for key, section_sentences in base["body_section_sentences"].items():
        matched = [s for s in section_sentences if re.search(regex, s.lower())]
        body_sentences.extend(matched)
        locations.extend([key] * len(matched))
    id_in_body = len(body_sentences) > 0

    return {
        "pmcid": pmcid,
        "job_id": job_id,
        "cite_count": cite_count,
        "title": base["title"],
        "abstract": base["abstract"],
        "doi": base["doi"],
        "pmid": base["pmid"],
        "year": base["year"],
        "journal": base["journal"],
        "author": base["author"],
        "rna_related": base["rna_related"],
        "probability": base["probability"],
        "id_in_title": id_in_title,
        "id_in_abstract": id_in_abstract,
        "id_in_body": id_in_body,
        "abstract_sentences": abstract_sentences,
        "body_sentences": body_sentences,
        "locations": locations,
        "score": len(abstract_sentences) + len(body_sentences),
        "retracted": False,
        "type": base["type"],
    }


def classify_abstract(abstract_text, rna_pipeline):
    probability = rna_pipeline.predict_proba([abstract_text])[0][1]
    rna_related = probability >= 0.5
    probability = round(float(probability), 2)

    return probability, rna_related


def scan_shard(xml_file_path, pmcid_to_jobs, rna_pipeline):
    """Parse a single .xml.gz shard once and produce article rows for every
    (pmcid, job_id) pair that wants it.

    `pmcid_to_jobs` is a dict: pmcid -> list[(job_id, cite_count)].
    Returns the four polars DataFrames in the same shape scan_job used to return.
    """
    start = datetime.datetime.now()
    n_requested_pmcids = len(pmcid_to_jobs)
    n_requested_rows = sum(len(v) for v in pmcid_to_jobs.values())
    logger.info(
        "scan_shard start xml=%s pmcids=%d jobs=%d",
        xml_file_path,
        n_requested_pmcids,
        n_requested_rows,
    )
    if xml_file_path is None:
        empty = pl.DataFrame()
        return empty, empty, empty, empty

    pending = set(pmcid_to_jobs.keys())
    articles = []

    with gzip.open(xml_file_path, "rb") as xml_file:
        context = ET.iterparse(xml_file, events=("end",), tag="article")
        for _event, elem in context:
            id_tag = elem.find('.//article-id[@pub-id-type="pmcid"]')
            if id_tag is not None:
                current_id = id_tag.text
                if current_id in pending:
                    if elem.find(".//trans-title-group") is not None:
                        logger.debug(
                            "Skipping %s: trans-title-group present", current_id
                        )
                        pending.remove(current_id)
                    else:
                        base = extract_article_base(elem, rna_pipeline)
                        for job_id, cite_count in pmcid_to_jobs[current_id]:
                            articles.append(
                                finalize_for_job(base, current_id, job_id, cite_count)
                            )
                        pending.remove(current_id)
                    if not pending:
                        elem.clear()
                        break
            elem.clear()
            # keep memory flat on large shards
            while elem.getprevious() is not None:
                del elem.getparent()[0]

    elapsed = (datetime.datetime.now() - start).total_seconds()
    found_pmcids = n_requested_pmcids - len(pending)

    if not articles:
        logger.warning(
            "xml=%s found=0/%d pmcids elapsed=%.1fs — no articles extracted",
            xml_file_path,
            n_requested_pmcids,
            elapsed,
        )
        empty = pl.DataFrame()
        return empty, empty, empty, empty

    articles_df = pl.DataFrame(articles)

    results_csv = articles_df.select(
        ["pmcid", "job_id", "id_in_title", "id_in_abstract", "id_in_body"]
    )
    articles_csv = articles_df.select(
        [
            "pmcid",
            "title",
            "abstract",
            "author",
            "pmid",
            "doi",
            "year",
            "journal",
            "score",
            "cite_count",
            "retracted",
            "rna_related",
            "probability",
            "type",
        ]
    )
    abstract_sentences = (
        articles_df.select(["pmcid", "job_id", "abstract_sentences"])
        .explode("abstract_sentences")
        .filter(pl.col("abstract_sentences").is_not_null())
    )
    body_sentences = (
        articles_df.select(["pmcid", "job_id", "body_sentences", "locations"])
        .explode("body_sentences", "locations")
        .filter(pl.col("body_sentences").is_not_null())
    )

    logger.info(
        "xml=%s found=%d/%d pmcids rows=%d elapsed=%.1fs",
        xml_file_path,
        found_pmcids,
        n_requested_pmcids,
        len(articles),
        elapsed,
    )

    return results_csv, articles_csv, abstract_sentences, body_sentences


def add_xml_shard(pmcid: str, xml_lookup: dict[int, Path]) -> str:
    if not xml_lookup:
        raise ValueError("xml_lookup is empty — no XML shards available")
    keys = sorted(xml_lookup.keys())
    pmc_number = int(re.match("PMC(\d+)", pmcid).group(1))
    candidates = [x for x in keys if pmc_number > x]
    if not candidates:
        xml_to_use = keys[0]
        logger.warning(
            "pmcid=%s (num=%d) is below smallest shard key %d; falling back to it",
            pmcid,
            pmc_number,
            xml_to_use,
        )
    else:
        xml_to_use = max(candidates)
    return str(xml_lookup[xml_to_use])


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Extract RNA mentions from ePMC OA articles. Writes loadable CSV files in cwd"
    )
    parser.add_argument(
        "--search-results",
        required=True,
        help="Output of the search results. A csv with no header, 5 columns: job_id, last_search_date, hit_count, pmcid, cite_count",
    )
    parser.add_argument(
        "--xml-directory",
        required=True,
        help="Path to the directory containing the ePMC XML dumps",
    )
    args = parser.parse_args()

    conn_str = os.environ["PSYCOPG_CONN"]
    model_path = os.environ["LITSCAN_MODEL"]

    rna_pipeline = joblib.load(model_path)

    conn = psycopg2.connect(conn_str)
    try:
        ## Skip searching if we searched this article for the same job_id before
        search_results = pl.read_csv(
            args.search_results,
            has_header=False,
            new_columns=[
                "job_id",
                "last_search_date",
                "hit_count",
                "pmcid",
                "cite_count",
            ],
        ).cast({"hit_count": pl.Int64, "cite_count": pl.Int64})
        if search_results.filter(pl.col("hit_count") > 0).height == 0:
            logger.warning(
                "Search results report no hits for anything, creating empty files and moving on"
            )
            Path("litscan_results.csv").touch()
            Path("litscan_articles.csv").touch()
            Path("litscan_abstract_sentences.csv").touch()
            Path("litscan_body_sentences.csv").touch()
            hit_counts = search_results.select(["job_id", "hit_count"])
            status = search_results.select("job_id").with_columns(
                status=pl.lit("success")
            )

            hit_counts.write_csv("litscan_hit_counts.csv", include_header=False)
            status.write_csv("litscan_job_status.csv", include_header=False)
            return

        logger.info("Loaded %d search result rows", search_results.height)

        previous_searches = get_previous_searches(
            conn, search_results["job_id"].to_list()
        )
        new_searches = search_results.join(
            previous_searches, on=["pmcid", "job_id"], how="anti"
        )
        logger.info(
            "Filtered previously-searched rows: %d -> %d",
            search_results.height,
            new_searches.height,
        )

        ## Build lookup for pmcid -> xml path
        xml_directory = Path(args.xml_directory)
        if not xml_directory.is_dir():
            raise FileNotFoundError(
                f"XML directory does not exist or is not a directory: {xml_directory}"
            )
        xml_files = list(xml_directory.glob("*.xml.gz"))
        if not xml_files:
            xml_files = list(xml_directory.rglob("*.xml.gz"))
        logger.info("Found %d XML shards in %s", len(xml_files), xml_directory)
        if not xml_files:
            raise FileNotFoundError(
                f"No *.xml.gz files found under {xml_directory} (checked non-recursive and recursive)"
            )

        regex = r"PMC(\d+).*"
        xml_lookup = {}
        for f in xml_files:
            m = re.match(regex, f.name)
            if m is None:
                logger.warning("Skipping XML file with unexpected name: %s", f.name)
                continue
            xml_lookup[int(m.group(1))] = f
        if not xml_lookup:
            raise FileNotFoundError(
                f"No XML files matching {regex!r} found under {xml_directory}"
            )
        add_xml_shard_partial = partial(add_xml_shard, xml_lookup=xml_lookup)
        scan_jobs = new_searches.with_columns(
            xml_path=pl.col("pmcid").map_elements(
                add_xml_shard_partial, return_dtype=str
            )
        )

        shard_groups = scan_jobs.group_by("xml_path", maintain_order=True).agg(
            [pl.col("pmcid"), pl.col("job_id"), pl.col("cite_count")]
        )
        logger.info(
            "Grouped into %d shard scan tasks (from %d job/pmcid rows)",
            shard_groups.height,
            scan_jobs.height,
        )

        with ExitStack() as stack:
            results_csv_fh = stack.enter_context(
                Path("litscan_results.csv").open(mode="a")
            )
            articles_csv_fh = stack.enter_context(
                Path("litscan_articles.csv").open(mode="a")
            )
            abstract_csv_fh = stack.enter_context(
                Path("litscan_abstract_sentences.csv").open(mode="a")
            )
            body_csv_fh = stack.enter_context(
                Path("litscan_body_sentences.csv").open(mode="a")
            )

            for row in shard_groups.iter_rows(named=True):
                xml_file_path = row["xml_path"]
                pmcid_to_jobs: dict[str, list[tuple[str, int]]] = {}
                for pmcid, job_id, cite_count in zip(
                    row["pmcid"], row["job_id"], row["cite_count"]
                ):
                    pmcid_to_jobs.setdefault(pmcid, []).append((job_id, cite_count))

                results, articles, abstract, body = scan_shard(
                    xml_file_path, pmcid_to_jobs, rna_pipeline
                )

                results.write_csv(results_csv_fh, include_header=False)
                articles.write_csv(articles_csv_fh, include_header=False)
                abstract.write_csv(abstract_csv_fh, include_header=False)
                body.write_csv(body_csv_fh, include_header=False)

        hit_counts = scan_jobs.select(["job_id", "hit_count"]).unique()
        hit_counts.write_csv("litscan_hit_counts.csv", include_header=False)
        status = (
            scan_jobs.select("job_id").unique().with_columns(status=pl.lit("success"))
        )
        status.write_csv("litscan_job_status.csv", include_header=False)
        logger.info("Wrote hit_counts and job_status CSVs")
    except NoDataError:
        logger.warning(
            "No data in search result CSV, no papers to scan for this batch!"
        )
        logger.error("No papers to search and can't create empty files, bailing!")
        return 1

    finally:
        conn.close()


if __name__ == "__main__":
    main()
