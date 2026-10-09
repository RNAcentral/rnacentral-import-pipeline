"""
Sankey plots of how each database's records flow into RNAcentral, one PNG per
database plus one for every source together, for the database pages.

docs/sankey-plots holds the hand-run ENA and Ensembl versions that also show the
provider's own totals.
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import psycopg2
from matplotlib.patches import PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MPath

LOGGER = logging.getLogger(__name__)

# One pass over xref for every database; dbid 0 is all of them, where a sequence
# several databases share is one unique sequence.
QUERY = """
with per_db as materialized (
  select dbid, urs, taxid, count(*) as records
  from xref
  where deleted = 'N'
  group by dbid, urs, taxid
),
imported as (
  select * from per_db
  union all
  select 0, urs, taxid, sum(records) from per_db group by urs, taxid
)
select
  e.dbid,
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
from imported e
left join qa_status q on q.urs = e.urs and q.taxid = e.taxid
left join rnc_rna_precomputed p on p.urs = e.urs and p.taxid = e.taxid
group by 1, 2, 3, 4
"""

# Config keys that are not rnc_database.descr once case and underscores are ignored.
ALIASES = {"PDB": "PDBE"}

TRACKING_SQL = """
CREATE TABLE IF NOT EXISTS rnacen.pipeline_tracking_sankey (
    database text PRIMARY KEY,
    last_run timestamptz NOT NULL
)
"""

# Loaded since its plot was last drawn, or never drawn.
DUE_QUERY = """
SELECT d.descr
FROM rnc_database d
LEFT JOIN rnacen.pipeline_tracking_sankey t ON t.database = d.descr
WHERE t.last_run IS NULL
   OR t.last_run < (SELECT max(end_time) FROM release_stats WHERE dbid = d.id)
"""

RECORD_SQL = """
INSERT INTO rnacen.pipeline_tracking_sankey (database, last_run)
SELECT unnest(%s::text[]), now()
ON CONFLICT (database) DO UPDATE SET last_run = excluded.last_run
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


def summarise(rows):
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
        self.lowest = 0

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
        self.lowest = min(self.lowest, y)
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


def draw(c, path):
    status = c["status"]
    in_rnac = sum(status.values())
    merged = c["records"] - in_rnac

    fig, ax = plt.subplots(figsize=(13, 7))
    ax.set(xlim=(0, 0.8), ylim=(-0.15, 1))
    ax.axis("off")
    sk = Sankey(ax)

    # Each stage gets its own scale; the zoom wedge shows the magnification.
    s1 = (TOP - GAP) / c["records"]
    sk.column(0.02, [("Imported records", c["records"], BLUE, "above")], s1)
    split = [("Identical copies merged", merged, GREY, "right")] if merged else []
    sk.column(0.13, split + [("Unique sequences", in_rnac, BLUE, "below")], s1)
    if merged:
        sk.link("Imported records", "Identical copies merged", merged, s1, GREY)
    sk.link("Imported records", "Unique sequences", in_rnac, s1, BLUE)

    states = [
        (k, status[k])
        for k in ("Passed QC", "Flagged by QC", "No QC result")
        if status.get(k)
    ]
    colours = {"Passed QC": AQUA, "Flagged by QC": ORANGE, "No QC result": GREY}
    types = [(name, v) for name, v in c["types"] if v]
    leaves = [(name, v, AQUA) for name, v in types] + [
        (name, v, ORANGE) for name, v in c["reasons"]
    ]
    s2 = (TOP - (len(leaves) - 1) * GAP) / in_rnac
    sk.column(0.30, [("Unique sequences in RNAcentral", in_rnac, BLUE, "above")], s2)
    sk.column(
        0.42, [(k, v, colours[k], "right") for k, v in states], s2, top=TOP - 0.5 * GAP
    )
    sk.column(0.56, [(n, v, col, "right") for n, v, col in leaves], s2)
    sk.zoom("Unique sequences", "Unique sequences in RNAcentral")
    for k, v in states:
        sk.link("Unique sequences in RNAcentral", k, v, s2, colours[k])
    for name, v in types:
        sk.link("Passed QC", name, v, s2, AQUA)
    for name, v in c["reasons"]:
        sk.link("Flagged by QC", name, v, s2, ORANGE)

    ax.text(
        0.02,
        # Labels of tiny leaves can stack below the plot.
        min(-0.12, sk.lowest - 0.08),
        "Identical sequences from one species become one entry. Flagged sequences "
        "stay in RNAcentral with a warning, each counted once under its most serious "
        "flag.",
        fontsize=8.5,
        color=MUTED,
        va="top",
        wrap=True,
    )
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def key(name):
    name = name.upper()
    return ALIASES.get(name, name).replace("_", "")


def write_all(db_url, out: Path, names=()):
    """
    One plot per database in names (config keys or rnc_database.descr), or else per
    database loaded since its plot was last drawn, plus one for all of them. Only
    the latter is tracked: plots drawn by hand are not the published ones.
    """
    with psycopg2.connect(db_url) as conn, conn.cursor() as cur:
        cur.execute("select id, descr from rnc_database")
        descrs = dict(cur.fetchall())
        if not names:
            cur.execute(TRACKING_SQL)
            cur.execute(DUE_QUERY)
            due = [descr for (descr,) in cur.fetchall()]
        cur.execute(QUERY)
        rows = cur.fetchall()

        wanted = {key(n) for n in names or due}
        if unknown := wanted - {key(d) for d in descrs.values()}:
            LOGGER.warning("No such database: %s", ", ".join(sorted(unknown)))

        by_db = {}
        for dbid, *row in rows:
            by_db.setdefault(dbid, []).append(row)
        drawn = []
        for dbid, db_rows in by_db.items():
            descr = descrs.get(dbid, "all")
            if dbid and key(descr) not in wanted:
                continue
            draw(summarise(db_rows), out / f"{descr.lower()}.png")
            if dbid:
                drawn.append(descr)

        if not names:
            cur.execute(RECORD_SQL, (drawn,))
