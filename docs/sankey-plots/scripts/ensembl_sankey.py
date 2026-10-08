#!/usr/bin/env python3
"""
Sankey of how Ensembl sequences flow into RNAcentral, for talks and papers.

Ensembl counts come from the public Ensembl MySQL servers, RNAcentral counts
from the database.

Writes ensembl-sankey.png to docs/sankey-plots.
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import psycopg2
import pymysql
from ena_sankey import AQUA, BLUE, GAP, GREY, MUTED, ORANGE, TOP, Sankey

OUT = Path(__file__).parent.parent

# Vertebrates, then Plants, Fungi, Metazoa and Protists (Bacteria is skipped).
SERVERS = (("ensembldb.ensembl.org", 3306), ("mysql-eg-publicsql.ebi.ac.uk", 4157))

NONCODING = {"snoncoding", "lnoncoding", "mnoncoding"}

# GENCODE (47) is left out: it repeats Ensembl's human and mouse transcripts.
QUERY = """
with ens as (
  select urs, taxid, count(*) as records
  from xref
  where dbid in (25, 31, 34, 35, 36) and deleted = 'N'
  group by urs, taxid
)
select
  case
    when q.urs is null then 'No QC result'
    when q.has_issue then 'Flagged by QC'
    else 'Passed QC'
  end,
  case
    when q.possible_contamination then 'Possible contamination'
    when q.possible_orf or q.possible_orf_stopfree or q.possible_orf_tcode then 'Possible ORF'
    when q.missing_rfam_match then 'Missing Rfam match'
    when q.incomplete_sequence then 'Incomplete sequence'
    when q.from_repetitive_region then 'Repetitive region'
  end,
  coalesce(p.rna_type, 'other'),
  count(*),
  sum(e.records)::bigint
from ens e
left join qa_status q on q.urs = e.urs and q.taxid = e.taxid
left join rnc_rna_precomputed p on p.urs = e.urs and p.taxid = e.taxid
group by 1, 2, 3
"""


def ensembl_counts():
    """Transcripts per Ensembl biotype group, over every latest-release core db."""
    counts = {}
    for host, port in SERVERS:
        conn = pymysql.connect(host=host, port=port, user="anonymous")
        with conn, conn.cursor() as cur:
            cur.execute("show databases like '%\\_core\\_%'")
            dbs = [r[0] for r in cur.fetchall()]
            release = max(int(d.split("_core_")[1].split("_")[0]) for d in dbs)
            for db in dbs:
                if f"_core_{release}_" not in db:
                    continue
                cur.execute(
                    f"select 1 from {db}.meta where meta_key = 'species.division'"
                    " and meta_value = 'EnsemblBacteria'"
                )
                if cur.fetchone():
                    continue
                cur.execute(
                    f"select b.biotype_group, count(*) from {db}.transcript t"
                    f" left join {db}.biotype b"
                    " on b.name = t.biotype and b.object_type = 'transcript'"
                    " group by 1"
                )
                for group, n in cur.fetchall():
                    counts[group] = counts.get(group, 0) + n
    return counts


def fetch():
    with psycopg2.connect(os.environ["PGDATABASE"]) as conn, conn.cursor() as cur:
        cur.execute(QUERY)
        rows = cur.fetchall()

    status, reasons, types, records = {}, {}, {}, 0
    for state, reason, rna_type, seqs, recs in rows:
        records += recs
        status[state] = status.get(state, 0) + seqs
        if state == "Flagged by QC":
            reasons[reason] = reasons.get(reason, 0) + seqs
        elif state == "Passed QC":
            name = rna_type.replace("_", " ")
            types[name] = types.get(name, 0) + seqs

    groups = ensembl_counts()
    top = sorted(types.items(), key=lambda kv: -kv[1])[:4]
    other = sum(types.values()) - sum(v for _, v in top)
    return {
        "coding": groups.get("coding", 0),
        "noncoding": sum(v for k, v in groups.items() if k in NONCODING),
        "other": sum(v for k, v in groups.items() if k not in NONCODING | {"coding"}),
        "records": records,
        "status": status,
        "reasons": sorted(reasons.items(), key=lambda kv: -kv[1]),
        "types": top + [("other types", other)],
    }


def draw(c):
    status = c["status"]
    in_rnac = sum(status.values())
    total = c["coding"] + c["other"] + c["noncoding"]
    skipped = c["noncoding"] - c["records"]
    merged = c["records"] - in_rnac

    fig, ax = plt.subplots(figsize=(16, 7))
    ax.set(xlim=(0, 1), ylim=(-0.15, 1))
    ax.axis("off")
    sk = Sankey(ax)

    # Each stage gets its own scale; the zoom wedges show the magnification.
    s1 = (TOP - 2 * GAP) / total
    sk.column(0.02, [("All Ensembl transcripts", total, GREY, "above")], s1)
    sk.column(
        0.13,
        [
            ("Protein-coding", c["coding"], GREY, "right"),
            ("Pseudogenes and other", c["other"], GREY, "right"),
            ("Non-coding", c["noncoding"], BLUE, "below"),
        ],
        s1,
    )
    sk.link("All Ensembl transcripts", "Protein-coding", c["coding"], s1, GREY)
    sk.link("All Ensembl transcripts", "Pseudogenes and other", c["other"], s1, GREY)
    sk.link("All Ensembl transcripts", "Non-coding", c["noncoding"], s1, BLUE)

    s2 = (TOP - 2 * GAP) / c["noncoding"]
    sk.column(
        0.30, [("Ensembl non-coding transcripts", c["noncoding"], BLUE, "above")], s2
    )
    sk.column(
        0.41,
        [
            ("Not imported", skipped, GREY, "right"),
            ("Identical copies merged", merged, GREY, "right"),
            ("Unique sequences", in_rnac, BLUE, "below"),
        ],
        s2,
    )
    sk.link("Ensembl non-coding transcripts", "Not imported", skipped, s2, GREY)
    sk.link(
        "Ensembl non-coding transcripts", "Identical copies merged", merged, s2, GREY
    )
    sk.link("Ensembl non-coding transcripts", "Unique sequences", in_rnac, s2, BLUE)
    sk.zoom("Non-coding", "Ensembl non-coding transcripts")

    states = [
        (k, status[k])
        for k in ("Passed QC", "Flagged by QC", "No QC result")
        if status.get(k)
    ]
    colours = {"Passed QC": AQUA, "Flagged by QC": ORANGE, "No QC result": GREY}
    leaves = [(name, v, AQUA) for name, v in c["types"]] + [
        (name, v, ORANGE) for name, v in c["reasons"]
    ]
    s3 = (TOP - (len(leaves) - 1) * GAP) / in_rnac
    sk.column(0.58, [("Unique sequences in RNAcentral", in_rnac, BLUE, "above")], s3)
    sk.column(
        0.70, [(k, v, colours[k], "right") for k, v in states], s3, top=TOP - 0.5 * GAP
    )
    sk.column(0.84, [(n, v, col, "right") for n, v, col in leaves], s3)
    sk.zoom("Unique sequences", "Unique sequences in RNAcentral")
    for k, v in states:
        sk.link("Unique sequences in RNAcentral", k, v, s3, colours[k])
    for name, v in c["types"]:
        sk.link("Passed QC", name, v, s3, AQUA)
    for name, v in c["reasons"]:
        sk.link("Flagged by QC", name, v, s3, ORANGE)

    ax.text(
        0.02,
        -0.12,
        "Transcripts are from every Ensembl and Ensembl Plants, Fungi, Metazoa and "
        "Protists genome, grouped by Ensembl's biotype groups. Not imported covers "
        "pseudogene-like and filtered records, species imported from their own "
        "databases (mouse strains, worm, fly, yeasts), and release differences. "
        "GENCODE is left out because it repeats Ensembl's human and mouse transcripts. "
        "Identical sequences from one species become one entry. Flagged sequences stay "
        "in RNAcentral with a warning, each counted once under its most serious flag.",
        fontsize=8.5,
        color=MUTED,
        va="top",
        wrap=True,
    )
    fig.savefig(OUT / "ensembl-sankey.png", dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    draw(fetch())
