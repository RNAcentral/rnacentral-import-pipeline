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


def extract_article(elem, rna_pipeline, regex):
    title = ""
    abstract = ""
    body = {}
    floats = {}

    title = extract_clean_text(elem.find(".//article-title"))
    abstract = extract_clean_text(elem.find(".//abstract"))
    body = parse_body_to_dict(elem.find(".//body"), ignore_tags=bad_tags)
    floats = parse_floats_to_dict(elem)

    ## Use regex to find stuff
    abstract_sentences = [
        s for s in nltk.sent_tokenize(abstract) if re.search(regex, s.lower())
    ]
    id_in_abstract = len(abstract_sentences) > 0
    id_in_title = bool(re.search(regex, title.lower()))

    body_sentences = []
    locations = []
    for key, section_text in body.items():
        section_sentences = [
            s for s in nltk.sent_tokenize(section_text) if re.search(regex, s.lower())
        ]
        body_sentences.extend(section_sentences)
        locations.extend([key] * len(section_sentences))

    id_in_body = len(body_sentences) > 0

    doi = elem.find('.//article-id[@pub-id-type="doi"]') or ""
    pmid = elem.find('.//article-id[@pub-id-type="pmid"]') or ""

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
        "abstract": nltk.sent_tokenize(abstract),
        "body": body,
        "floats": floats,
        "doi": doi,
        "pmid": pmid,
        "year": year,
        "journal": journal,
        "author": author,
        "rna_related": rna_related,
        "probability": probability,
        "id_in_title": id_in_title,
        "id_in_abstract": id_in_abstract,
        "id_in_body": id_in_body,
        "abstract_sentences": abstract_sentences,
        "body_sentences": body_sentences,
        "locations": locations,
        "score": len(abstract_sentences) + len(body_sentences),
        "retracted": False,
        "type": article_type,
    }


def classify_abstract(abstract_text, rna_pipeline):
    probability = rna_pipeline.predict_proba([abstract_text])[0][1]
    rna_related = probability >= 0.5
    probability = round(float(probability), 2)

    return probability, rna_related


def scan_job(job_id, pmcid_list, cite_counts, xml_file_path, rna_pipeline):
    start = datetime.datetime.now()
    n_requested = len(pmcid_list)
    logger.info(
        "scan_job start job_id=%s pmcids=%s xml=%s",
        job_id,
        pmcid_list,
        xml_file_path,
    )
    if xml_file_path is None:
        empty = pl.DataFrame()
        return empty, empty, empty, empty

    context = ET.iterparse(xml_file_path, events=("end",))
    cite_lookup = {pmcid: cc for pmcid, cc in zip(pmcid_list, cite_counts)}
    pmcid_list = set(pmcid_list)
    articles = []

    # word-boundary regex matching the job_id in article text
    regex = (
        r"(^|\s|\(|\u201c|\u2018|\u201d|\;)"
        + re.escape(job_id.lower())
        + r"($|[\s.,:;?\u2019\u201d\u201c\"/)])"
    )

    for event, elem in context:

        if elem.tag == "article":
            id_tag = elem.find('.//article-id[@pub-id-type="pmcid"]')
            if id_tag is not None:
                current_id = id_tag.text
                if current_id in pmcid_list:
                    logger.debug("Found %s for job_id=%s", current_id, job_id)
                    if elem.find(".//trans-title-group") is not None:
                        logger.debug(
                            "Skipping %s: trans-title-group present", current_id
                        )
                        continue
                    article = extract_article(elem, rna_pipeline, regex)
                    article["job_id"] = job_id
                    article["cite_count"] = cite_lookup[current_id]
                    articles.append(article)
                    pmcid_list.remove(current_id)
                    if len(pmcid_list) == 0:
                        break
            elem.clear()

    elapsed = (datetime.datetime.now() - start).total_seconds()
    found = n_requested - len(pmcid_list)

    if not articles:
        logger.warning(
            "job_id=%s found=0/%d elapsed=%.1fs — no articles extracted",
            job_id,
            n_requested,
            elapsed,
        )
        empty = pl.DataFrame()
        return empty, empty, empty, empty

    ## Now the articles list contains a load of dicts, which we can convert to a dataframe, then delete
    articles = pl.DataFrame(articles)

    ## Now slice the dataframe up to create the various CSVs that will be loaded

    results_csv = articles.select(
        ["pmcid", "job_id", "id_in_title", "id_in_abstract", "id_in_body"]
    )
    articles_csv = articles.select(
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
    abstract_sentences = articles.select(
        ["pmcid", "job_id", "abstract_sentences"]
    ).explode("abstract_sentences")
    body_sentences = articles.select(
        ["pmcid", "job_id", "body_sentences", "locations"]
    ).explode("body_sentences", "locations")

    logger.info(
        "job_id=%s found=%d/%d elapsed=%.1fs",
        job_id,
        found,
        n_requested,
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
        help="Output of the search results. A csv with no header, 3 columns: job_id, pmcid, cite count",
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
        xml_lookup = {int(re.match(regex, f.name).group(1)): f for f in xml_files}
        add_xml_shard_partial = partial(add_xml_shard, xml_lookup=xml_lookup)
        scan_jobs = new_searches.with_columns(
            xml_path=pl.col("pmcid").map_elements(
                add_xml_shard_partial, return_dtype=str
            )
        )

        scan_jobs_grouped = scan_jobs.group_by(
            "job_id", "xml_path", maintain_order=True
        ).agg([pl.col("pmcid"), pl.col("cite_count"), pl.first("hit_count")])
        logger.info(
            "Grouped into %d (job_id, xml_shard) scan tasks",
            scan_jobs_grouped.height,
        )

        ## scan_jobs_grouped is now 1 row -> 1 JobID + xml file + N PMCIDs
        ## Iterate over rows and accumulate results

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

            for row in scan_jobs_grouped.iter_rows(named=True):
                job_id = row["job_id"]
                xml_file_path = row["xml_path"]

                pmcid_list = row["pmcid"]
                cite_counts = row["cite_count"]
                results, articles, abstract, body = scan_job(
                    job_id, pmcid_list, cite_counts, xml_file_path, rna_pipeline
                )

                results.write_csv(results_csv_fh, include_header=False)
                articles.write_csv(articles_csv_fh, include_header=False)
                abstract.write_csv(abstract_csv_fh, include_header=False)
                body.write_csv(body_csv_fh, include_header=False)

        hit_counts = scan_jobs.select(["job_id", "hit_count"])
        hit_counts.write_csv("litscan_hit_counts.csv", include_header=False)
        status = scan_jobs.select("job_id").with_columns(status=pl.lit("success"))
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
