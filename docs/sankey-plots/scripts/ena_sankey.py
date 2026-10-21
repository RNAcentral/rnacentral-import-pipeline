#!/usr/bin/env python3
"""
Sankey of how ENA sequences flow into RNAcentral, for talks and papers.

ENA counts come from the ENA portal API, RNAcentral counts from the database.

Writes ena-sankey.png to docs/sankey-plots.
"""

import os
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import psycopg2
from matplotlib.patches import PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MPath

OUT = Path(__file__).parent.parent

QUERY = """
with ena as (
  select urs, taxid, count(*) as records
  from xref
  where dbid = 1 and deleted = 'N'
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
from ena e
left join qa_status q on q.urs = e.urs and q.taxid = e.taxid
left join rnc_rna_precomputed p on p.urs = e.urs and p.taxid = e.taxid
group by 1, 2, 3
"""

GREY = "#b4b2a9"
BLUE = "#2a78d6"
AQUA = "#1baf7a"
ORANGE = "#eb6834"
INK = "#0b0b0b"
MUTED = "#52514e"

TOP = 0.88
GAP = 0.025
WIDTH = 0.010


def ena_count(result):
    url = f"https://www.ebi.ac.uk/ena/portal/api/count?result={result}"
    with urllib.request.urlopen(url, timeout=600) as response:
        return int(response.read().split()[-1])


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

    top = sorted(types.items(), key=lambda kv: -kv[1])[:4]
    other = sum(types.values()) - sum(v for _, v in top)
    return {
        "coding": ena_count("coding"),
        "noncoding": ena_count("noncoding"),
        "records": records,
        "status": status,
        "reasons": sorted(reasons.items(), key=lambda kv: -kv[1]),
        "types": top + [("other types", other)],
    }


def human(n):
    for size, unit in ((1e9, "B"), (1e6, "M"), (1e3, "k")):
        if n >= size:
            return f"{n / size:.3g} {unit}"
    return str(n)


class Sankey:
    def __init__(self, ax):
        self.ax = ax
        self.nodes = {}

    def column(self, x, items, scale, top=TOP):
        """Stack (name, value, colour, label side) nodes downward from top."""
        y, label_y = top, 1.0
        for name, value, colour, side in items:
            h = value * scale
            self.nodes[name] = [x, y - h, y, y, y]  # x, bottom, top, out, in
            self.ax.add_patch(Rectangle((x, y - h), WIDTH, h, color=colour, lw=0))
            if side == "above":
                self.label(name, value, x, TOP + 0.05)
            elif side == "below":
                self.label(name, value, x, y - h - 0.04)
            else:
                # Small neighbouring nodes would otherwise stack their labels.
                label_y = min(y - h / 2, label_y - 0.065)
                self.label(name, value, x + WIDTH + 0.005, label_y)
            y -= h + GAP

    def label(self, name, value, x, y):
        box = dict(facecolor="white", edgecolor="none", pad=1, alpha=0.85)
        self.ax.text(
            x, y, name, va="bottom", fontsize=10, color=INK, fontweight="bold", bbox=box
        )
        self.ax.text(
            x, y - 0.004, human(value), va="top", fontsize=10, color=MUTED, bbox=box
        )

    def link(self, src, dst, value, scale, colour):
        s, d = self.nodes[src], self.nodes[dst]
        h = value * scale
        x0, x1 = s[0] + WIDTH, d[0]
        a0, a1 = s[3], s[3] - h
        b0, b1 = d[4], d[4] - h
        s[3], d[4] = a1, b1
        mid = (x0 + x1) / 2
        path = MPath(
            [
                (x0, a0),
                (mid, a0),
                (mid, b0),
                (x1, b0),
                (x1, b1),
                (mid, b1),
                (mid, a1),
                (x0, a1),
                (x0, a0),
            ],
            [
                MPath.MOVETO,
                MPath.CURVE4,
                MPath.CURVE4,
                MPath.CURVE4,
                MPath.LINETO,
                MPath.CURVE4,
                MPath.CURVE4,
                MPath.CURVE4,
                MPath.CLOSEPOLY,
            ],
        )
        self.ax.add_patch(PathPatch(path, color=colour, alpha=0.35, lw=0))

    def zoom(self, src, dst):
        """Wedge from a thin node to the full-height node that magnifies it."""
        s, d = self.nodes[src], self.nodes[dst]
        corners = [
            (s[0] + WIDTH, s[2]),
            (d[0], d[2]),
            (d[0], d[1]),
            (s[0] + WIDTH, s[1]),
        ]
        self.ax.add_patch(Polygon(corners, color="#e9e8e4", lw=0, zorder=0))
        for (xa, ya), (xb, yb) in ((corners[0], corners[1]), (corners[3], corners[2])):
            self.ax.plot([xa, xb], [ya, yb], ls=(0, (3, 3)), lw=0.8, color=MUTED)


def draw(c):
    status = c["status"]
    in_rnac = sum(status.values())
    total = c["coding"] + c["noncoding"]
    filtered = c["noncoding"] - c["records"]
    merged = c["records"] - in_rnac

    fig, ax = plt.subplots(figsize=(16, 7))
    ax.set(xlim=(0, 1), ylim=(-0.15, 1))
    ax.axis("off")
    sk = Sankey(ax)

    # Each stage gets its own scale; the zoom wedges show the magnification.
    s1 = (TOP - GAP) / total
    sk.column(0.02, [("All ENA annotated sequences", total, GREY, "above")], s1)
    sk.column(
        0.13,
        [
            ("Protein-coding", c["coding"], GREY, "right"),
            ("Non-coding", c["noncoding"], BLUE, "below"),
        ],
        s1,
    )
    sk.link("All ENA annotated sequences", "Protein-coding", c["coding"], s1, GREY)
    sk.link("All ENA annotated sequences", "Non-coding", c["noncoding"], s1, BLUE)

    s2 = (TOP - 2 * GAP) / c["noncoding"]
    sk.column(0.30, [("ENA non-coding sequences", c["noncoding"], BLUE, "above")], s2)
    sk.column(
        0.41,
        [
            ("Removed by import filters", filtered, GREY, "right"),
            ("Identical copies merged", merged, GREY, "right"),
            ("Unique sequences", in_rnac, BLUE, "below"),
        ],
        s2,
    )
    sk.link("ENA non-coding sequences", "Removed by import filters", filtered, s2, GREY)
    sk.link("ENA non-coding sequences", "Identical copies merged", merged, s2, GREY)
    sk.link("ENA non-coding sequences", "Unique sequences", in_rnac, s2, BLUE)
    sk.zoom("Non-coding", "ENA non-coding sequences")

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
        -0.15,
        "Import filters drop sequences under 10 nt or with more than 10% N, mRNAs, "
        "protein-like and pseudogene records, and metagenomic small-subunit rRNA that "
        "fails ribotyper. Identical sequences from one species become one entry. "
        "Flagged sequences stay in RNAcentral with a warning, each counted once under "
        "its most serious flag. ENA counts are live and RNAcentral counts are from the "
        "current release, so the removed figure is approximate.",
        fontsize=8.5,
        color=MUTED,
        wrap=True,
    )
    fig.savefig(OUT / "ena-sankey.png", dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    draw(fetch())
