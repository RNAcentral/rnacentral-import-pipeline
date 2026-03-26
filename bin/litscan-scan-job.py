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
import sys
import time

import joblib
import nltk
import psycopg2
import requests
from xml.etree import ElementTree as ET
from xml.etree.ElementTree import ParseError

# Download NLTK tokenizer data if not already present
nltk.download("punkt", quiet=True)
nltk.download("punkt_tab", quiet=True)

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

EUROPE_PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"


# ---------------------------------------------------------------------------
# Text helpers  (from rnacentral-references/training/export_data.py)
# ---------------------------------------------------------------------------

def clean_text(text):
    text = text.lower()
    text = re.sub(r"<[^>]*>", " ", text)
    text = re.sub(r"\[.*?\]", "", text)
    text = re.sub(r"https?://\S+|www\.\S+", "", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


# ---------------------------------------------------------------------------
# XML helpers  (from rnacentral-references/consumer/views/submit_job.py)
# ---------------------------------------------------------------------------

_AVOID_TAGS = {
    "xref", "ext-link", "media", "caption", "monospace", "label",
    "disp-formula", "inline-formula", "inline-graphic", "def", "def-list",
    "def-item", "term", "funding-source", "award-id", "graphic",
    "alternatives", "tex-math", "sec-meta", "kwd-group", "kwd", "object-id",
    "{http://www.w3.org/1998/Math/MathML}math",
    "{http://www.w3.org/1998/Math/MathML}mrow",
    "{http://www.w3.org/1998/Math/MathML}mi",
    "{http://www.w3.org/1998/Math/MathML}mo",
    "{http://www.w3.org/1998/Math/MathML}msub",
    "{http://www.w3.org/1998/Math/MathML}mn",
    "{http://www.w3.org/1998/Math/MathML}msup",
    "{http://www.w3.org/1998/Math/MathML}mtext",
    "{http://www.w3.org/1998/Math/MathML}msubsup",
    "{http://www.w3.org/1998/Math/MathML}mover",
    "{http://www.w3.org/1998/Math/MathML}mstyle",
    "{http://www.w3.org/1998/Math/MathML}munderover",
    "{http://www.w3.org/1998/Math/MathML}mspace",
    "{http://www.w3.org/1998/Math/MathML}mfenced",
    "{http://www.w3.org/1998/Math/MathML}mpadded",
    "{http://www.w3.org/1998/Math/MathML}mfrac",
    "{http://www.w3.org/1998/Math/MathML}msqrt",
}


def get_text(sec):
    sec_sentences = [
        "".join(item.itertext())
        for item in sec.iter(tag="p")
        if item.text and item.tag not in _AVOID_TAGS
    ]
    sec_sentences = [" ".join(s.split()) for s in sec_sentences if len(s.split()) > 1]
    return " ".join(sec_sentences) if sec_sentences else ""


def get_sections(tree):
    sections = tree.findall("./body/sec")
    section_map = {}
    count = 0
    for sec in sections:
        sec_title = sec.find("title")
        try:
            t = sec_title.text.lower()
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
        except AttributeError:
            key = "other"
        section_map[key + str(count)] = sec
        count += 1
    return section_map


# ---------------------------------------------------------------------------
# Europe PMC API
# ---------------------------------------------------------------------------

def articles_list(job_id, query_filter, date, page="*"):
    search_date = (
        f" AND (FIRST_PDATE:[{date} TO {datetime.date.today().strftime('%Y-%m-%d')}])"
        if date else ""
    )
    qf = f" AND {query_filter}" if query_filter else ""
    query = (
        f'search?query=("{job_id}"{qf} AND IN_EPMC:Y AND OPEN_ACCESS:Y AND NOT SRC:PPR{search_date})'
        f"&sort_date:y&pageSize=500&cursorMark={page}"
    )
    try:
        response = requests.get(EUROPE_PMC + query, timeout=60)
        response.raise_for_status()
        articles = response.text
    except requests.exceptions.RequestException as e:
        logger.warning("Error fetching article list from Europe PMC: %s", e)
        return [], None

    try:
        root = ET.fromstring(articles)
    except ParseError:
        return [], None

    pmcid_list = [
        {"pmcid": item.find("pmcid").text, "cited_by": item.find("citedByCount").text}
        for item in root.findall("./resultList/result")
        if item.find("pmcid") is not None and item.find("citedByCount") is not None
    ]
    try:
        next_page = root.find("nextCursorMark").text
    except AttributeError:
        next_page = None

    return pmcid_list, next_page


# ---------------------------------------------------------------------------
# DB helpers (synchronous psycopg2)
# ---------------------------------------------------------------------------

def get_query_and_limit(cur, job_id):
    cur.execute("SELECT query, search_limit FROM litscan_job WHERE job_id = %s", (job_id,))
    row = cur.fetchone()
    return (row[0], row[1]) if row else (None, None)


def get_search_date(cur, job_id):
    cur.execute("SELECT finished FROM litscan_job WHERE job_id = %s", (job_id,))
    row = cur.fetchone()
    return row[0] if row else None


def set_job_status(cur, job_id, status):
    if status in ("success", "error"):
        cur.execute(
            "UPDATE litscan_job SET status = %s, finished = %s WHERE job_id = %s",
            (status, datetime.datetime.now(), job_id),
        )
    else:
        cur.execute(
            "UPDATE litscan_job SET status = %s WHERE job_id = %s",
            (status, job_id),
        )


def get_pmcid_in_result(cur, job_id):
    cur.execute("SELECT pmcid FROM litscan_result WHERE job_id = %s", (job_id,))
    return [row[0] for row in cur.fetchall()]


def get_pmcid(cur, pmcid):
    cur.execute("SELECT pmcid FROM litscan_article WHERE pmcid = %s", (pmcid,))
    row = cur.fetchone()
    return row[0] if row else None


def save_article(cur, article):
    cur.execute(
        """
        INSERT INTO litscan_article
            (pmcid, title, abstract, author, pmid, doi, year, journal,
             score, cited_by, retracted, rna_related, probability, type)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT (pmcid) DO NOTHING
        """,
        (
            article["pmcid"], article["title"], article.get("abstract", ""),
            article["author"], article["pmid"], article["doi"],
            article["year"], article["journal"], article["score"],
            article["cited_by"], article["retracted"],
            article["rna_related"], article["probability"],
            article.get("type", ""),
        ),
    )


def save_result(cur, result):
    cur.execute(
        """
        INSERT INTO litscan_result (pmcid, job_id, id_in_title, id_in_abstract, id_in_body)
        VALUES (%s, %s, %s, %s, %s)
        ON CONFLICT (pmcid, job_id) DO NOTHING
        RETURNING id
        """,
        (
            result["pmcid"], result["job_id"], result["id_in_title"],
            result["id_in_abstract"], result["id_in_body"],
        ),
    )
    row = cur.fetchone()
    return row[0] if row else None


def save_abstract_sentences(cur, sentences):
    cur.executemany(
        "INSERT INTO litscan_abstract_sentence (result_id, sentence) VALUES (%s, %s)",
        [(s["result_id"], s["sentence"]) for s in sentences],
    )


def save_body_sentences(cur, sentences):
    cur.executemany(
        "INSERT INTO litscan_body_sentence (result_id, sentence, location) VALUES (%s, %s, %s)",
        [(s["result_id"], s["sentence"], s["location"]) for s in sentences],
    )


def get_hit_count(cur, job_id):
    cur.execute("SELECT hit_count FROM litscan_job WHERE job_id = %s", (job_id,))
    row = cur.fetchone()
    return row[0] if (row and row[0] is not None) else 0


def save_hit_count(cur, job_id, hit_count):
    cur.execute(
        "UPDATE litscan_job SET hit_count = %s WHERE job_id = %s",
        (hit_count, job_id),
    )


# ---------------------------------------------------------------------------
# Core scan logic  (ported from seek_references() in submit_job.py)
# ---------------------------------------------------------------------------

def scan_job(job_id, conn, rna_pipeline):
    start = datetime.datetime.now()
    # word-boundary regex matching the job_id in article text
    regex = (
        r"(^|\s|\(|\u201c|\u2018|\u201d|\;)"
        + re.escape(job_id.lower())
        + r"($|[\s.,:;?\u2019\u201d\u201c\"/)])"
    )
    pmcid_list = []
    hit_count = 0

    with conn.cursor() as cur:
        set_job_status(cur, job_id.lower(), "started")
        conn.commit()

        last_search = get_search_date(cur, job_id.lower())
        date = last_search.date() if last_search else None

        query_filter, search_limit = get_query_and_limit(cur, job_id.lower())
        if query_filter:
            query_filter = query_filter.lower().replace(job_id.lower(), "")
        search_limit = search_limit if search_limit else 1_000_000

        # ---- collect article list from Europe PMC ----
        temp_list, next_page = articles_list(job_id, query_filter, date)
        for item in temp_list:
            if len(pmcid_list) < search_limit and item not in pmcid_list:
                pmcid_list.append(item)

        while len(pmcid_list) < search_limit and next_page:
            temp_list, next_page = articles_list(job_id, query_filter, date, next_page)
            for item in temp_list:
                if len(pmcid_list) < search_limit and item not in pmcid_list:
                    pmcid_list.append(item)

        # on re-scan, skip articles already processed for this job
        if date and pmcid_list:
            pmcid_in_db = get_pmcid_in_result(cur, job_id.lower())
            pmcid_list = [item for item in pmcid_list if item["pmcid"] not in pmcid_in_db]

        # ---- process each article ----
        for element in pmcid_list:
            time.sleep(0.6)  # respect Europe PMC rate limit (10 req/s)

            try:
                response = requests.get(
                    EUROPE_PMC + element["pmcid"] + "/fullTextXML", timeout=60
                )
                response.raise_for_status()
                get_article = response.text
            except requests.exceptions.RequestException as e:
                logger.warning("Error fetching %s: %s", element["pmcid"], e)
                get_article = None

            if not get_article:
                continue

            # extract text sections
            abstract_txt = re.search(r"<abstract(.*?)</abstract>", get_article, re.DOTALL)
            body_txt = re.search(r"<body(.*?)</body>", get_article, re.DOTALL)
            floats_txt = re.search(r"<floats-group(.*?)</floats-group>", get_article, re.DOTALL)

            if abstract_txt and body_txt and floats_txt:
                full_txt = abstract_txt[0] + body_txt[0] + floats_txt[0]
            elif abstract_txt and body_txt:
                full_txt = abstract_txt[0] + body_txt[0]
            elif body_txt:
                full_txt = body_txt[0]
            else:
                continue

            # quick check: job_id actually appears in the text
            full_txt_no_tags = re.sub(r"<[^>]*>", " ", full_txt)
            if not re.search(regex, full_txt_no_tags.lower()):
                continue

            # remove tables, figures, supplementary material before full parse
            full_txt = re.sub(
                r"(?is)<(counts|table-wrap|table|fig-group|fig|supplementary-material).*?>.*?(</\1>)",
                "",
                get_article,
            )

            try:
                article = ET.fromstring(full_txt)
            except ParseError as e:
                logger.warning("XML parse error for %s: %s", element["pmcid"], e)
                continue

            # skip non-English articles
            if article.find("./front/article-meta/title-group/trans-title-group") is not None:
                continue

            # title
            get_title = article.find("./front/article-meta/title-group/article-title")
            try:
                title = "".join(get_title.itertext()).strip()
            except AttributeError:
                continue

            result_response = {
                "id_in_title": job_id.lower() in title.lower(),
                "job_id": job_id.lower(),
                "pmcid": element["pmcid"],
            }

            # abstract
            abstract_types = [
                "teaser", "web-summary", "summary", "precis", "graphical", "author-highlights"
            ]
            get_abstract_tags = [
                item for item in article.findall(".//abstract")
                if not any(v in item.attrib.values() for v in abstract_types)
            ]
            abstract_text = [" ".join(item.itertext()) for item in get_abstract_tags]
            abstract = " ".join(abstract_text).replace(" .", ".").replace("  ", " ")

            abstract_sentences = [
                s for s in nltk.sent_tokenize(abstract)
                if re.search(regex, s.lower())
            ]
            result_response["id_in_abstract"] = bool(abstract_sentences)

            # body sentences
            sections = get_sections(article)
            body_sentences = {}
            for section_name, section in sections.items():
                item_text = get_text(section)
                tokenized = nltk.sent_tokenize(item_text)
                body_sentences[section_name] = []
                for idx, sentence in enumerate(tokenized):
                    if re.search(regex, sentence.lower()) and len(sentence.split()) > 3:
                        prev_s = tokenized[idx - 1] if idx > 0 else None
                        next_s = tokenized[idx + 1] if idx < len(tokenized) - 1 else None
                        if prev_s and next_s:
                            body_sentences[section_name].append(prev_s + " " + sentence + " " + next_s)
                        elif prev_s:
                            body_sentences[section_name].append(prev_s + " " + sentence)
                        elif next_s:
                            body_sentences[section_name].append(sentence + " " + next_s)
                        else:
                            body_sentences[section_name].append(sentence)

            if abstract_sentences and not any(body_sentences.values()):
                result_response["id_in_body"] = False
            elif not abstract_sentences and not any(body_sentences.values()):
                result_response["id_in_body"] = True
                body_sentences["other"] = [
                    "%s found in an image, table or supplementary material" % job_id
                ]
            else:
                result_response["id_in_body"] = True

            # save article metadata if not already in DB
            article_in_db = get_pmcid(cur, element["pmcid"])
            if not article_in_db:
                article_response = {
                    "pmcid": element["pmcid"],
                    "title": title,
                    "abstract": abstract,
                    "retracted": False,
                    "score": len(abstract_sentences) + len(body_sentences),
                }

                try:
                    cited_by = int(element.get("cited_by") or 0)
                except (ValueError, TypeError):
                    cited_by = 0
                article_response["cited_by"] = cited_by

                article_response["type"] = (
                    article.attrib.get("article-type", "").strip().replace("-", " ").capitalize()
                )

                # authors
                article_response["author"] = ""
                contrib_group = article.find("./front/article-meta/contrib-group")
                if contrib_group is not None:
                    authors = []
                    for auth in contrib_group.findall(".//name"):
                        surname_el = auth.find("surname")
                        given_el = auth.find("given-names")
                        surname = (surname_el.text or "") if surname_el is not None else ""
                        given = (given_el.text or "") if given_el is not None else ""
                        if surname and given:
                            authors.append(f"{surname}, {given}")
                        elif surname or given:
                            authors.append(surname + given)
                    article_response["author"] = "; ".join(authors)

                # pmid and doi
                article_response["doi"] = ""
                article_response["pmid"] = ""
                for item in article.findall("./front/article-meta/article-id"):
                    if item.attrib == {"pub-id-type": "doi"}:
                        article_response["doi"] = item.text or ""
                    elif item.attrib == {"pub-id-type": "pmid"}:
                        article_response["pmid"] = item.text or ""

                # year
                article_response["year"] = 0
                for item in article.findall("./front/article-meta/pub-date"):
                    if {"epub", "ppub", "pub"}.intersection(item.attrib.values()):
                        year_el = item.find("year")
                        if year_el is not None and year_el.text:
                            article_response["year"] = int(year_el.text)

                # journal
                article_response["journal"] = ""
                journal_el = article.find("./front/journal-meta/journal-title-group/journal-title")
                if journal_el is None:
                    journal_el = article.find("./front/journal-meta/journal-title")
                if journal_el is not None:
                    article_response["journal"] = journal_el.text or ""

                # ML classification
                cleaned = clean_text(abstract)
                relevance_label = rna_pipeline.predict([cleaned])[0]
                article_response["rna_related"] = bool(int(relevance_label))
                probability = rna_pipeline.predict_proba([cleaned])[0][1]
                article_response["probability"] = round(float(probability), 2)

                save_article(cur, article_response)

            # save result and sentences
            result_id = save_result(cur, result_response)
            conn.commit()

            if result_id:
                if abstract_sentences:
                    save_abstract_sentences(
                        cur,
                        [{"result_id": result_id, "sentence": s} for s in abstract_sentences],
                    )

                loc_prefix_map = {
                    "intro": "intro", "results": "results", "discussion": "discussion",
                    "conclusion": "conclusion", "method": "method",
                }
                body_to_save = []
                for loc, sentences in body_sentences.items():
                    location = next(
                        (v for k, v in loc_prefix_map.items() if loc.startswith(k)),
                        "other",
                    )
                    for s in sentences:
                        body_to_save.append({"result_id": result_id, "sentence": s, "location": location})
                if body_to_save:
                    save_body_sentences(cur, body_to_save)

                conn.commit()

            hit_count += 1

        # ---- finalise job ----
        if date:
            hit_count += get_hit_count(cur, job_id.lower())
        save_hit_count(cur, job_id.lower(), hit_count)
        set_job_status(cur, job_id.lower(), "success")
        conn.commit()

    elapsed = (datetime.datetime.now() - start).total_seconds()
    print(f"job_id={job_id} hit_count={hit_count} elapsed={elapsed:.1f}s")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Scan an RNA ID against Europe PMC and store results in the litscan DB"
    )
    parser.add_argument("--job-id", required=True, help="The RNA ID / job_id to scan")
    args = parser.parse_args()

    conn_str = os.environ["PSYCOPG_CONN"]
    model_path = os.environ["LITSCAN_MODEL"]

    rna_pipeline = joblib.load(model_path)

    conn = psycopg2.connect(conn_str)
    try:
        scan_job(args.job_id, conn, rna_pipeline)
    except Exception:
        try:
            with conn.cursor() as cur:
                set_job_status(cur, args.job_id.lower(), "error")
                conn.commit()
        except Exception:
            pass
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    main()
