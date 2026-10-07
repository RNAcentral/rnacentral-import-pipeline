# -*- coding: utf-8 -*-

"""
Forward merge of gene calls across releases, based on shared members.

A previous gene continues as a next-release gene when most of the previous
gene's members (transcript region names) are still together in that gene. This
replaces matching on gene extents, which lost identity whenever a short gene's
ends moved, followed a stale first-exon-only extent, and could never notice a
split or a merge.

Rules:

1. Members whose URS has been retired are dropped from the previous genes first;
   a gene left with no members is retired.
2. Every (previous, next) pair is scored by the number of members they share.
   Pairs are taken strongest first, each gene used at most once, and a pair is
   accepted only if it keeps at least `min_keep` of the *previous* gene's
   members. A small fragment of a large gene therefore cannot take its
   identity; it becomes a new gene that records `split_from`.
3. The heir keeps the previous name and internal name. Its extent is recomputed
   from its members' full spans, and its version increases only if the
   membership changed.
4. History is kept in one list column, `history`, of {event, internal_name, name}
   (the previous gene's stable key and its `public_name.version`). `merged_from`:
   other previous genes that contributed at least `link_share` of their members
   to this gene; `split_from`: for a gene with no heir, the previous genes it was
   carved out of.
5. Previous genes with no heir are retired (not in the output).
"""

import collections
import re
import typing as ty

EXON = re.compile(r"(\d+)-(\d+)")


def region_urs_taxid(region_name: str) -> str:
    return region_name.split("@", 1)[0]


def region_urs(region_name: str) -> str:
    """The bare URS, as listed in the release's inactive-ids file."""
    return region_urs_taxid(region_name).split("_", 1)[0]


def region_span(region_name: str) -> ty.Tuple[int, int]:
    """Full span (first exon start, last exon stop) of 'URS@chrom/s-e,s-e:strand'."""
    locus = region_name.rsplit(":", 1)[0].rsplit("/", 1)[-1]
    exons = [(int(a), int(b)) for a, b in EXON.findall(locus)]
    return min(a for a, _ in exons), max(b for _, b in exons)


def true_extent(members: ty.Iterable[str]) -> ty.Tuple[int, int]:
    spans = [region_span(m) for m in members]
    return min(a for a, _ in spans), max(b for _, b in spans)


def split_name(name: str) -> ty.Tuple[str, int]:
    base, _, version = name.rpartition(".")
    return base, int(version)


def merge(
    previous: ty.List[dict],
    following: ty.List[dict],
    retired_urs: ty.AbstractSet[str] = frozenset(),
    min_keep: float = 0.5,
    link_share: float = 0.2,
) -> ty.Tuple[ty.List[dict], ty.Counter]:
    """
    Merge `previous` (the stored active genes) with `following` (the new calls).

    Both are lists of dicts with at least: name (`base.version`), internal_name,
    assembly_id, members (list of region names), first_release, last_release.
    Returns (genes, event counts); each returned gene also has: version,
    start, stop, history. `name` of the result has the new
    version appended and `last_release` is NOT changed here (the caller stamps it).
    """
    kept = []
    events: ty.Counter = collections.Counter()
    for gene in previous:
        members = frozenset(
            m for m in gene["members"] if region_urs(m) not in retired_urs
        )
        if members:
            kept.append(
                {
                    **gene,
                    "members": members,
                    "stored_members": frozenset(gene["members"]),
                }
            )
        else:
            events["retired"] += 1
    nxt = [{**g, "members": frozenset(g["members"])} for g in following]

    by_member: ty.Dict[str, ty.List[int]] = collections.defaultdict(list)
    for j, gene in enumerate(nxt):
        for member in gene["members"]:
            by_member[member].append(j)
    shared: ty.Dict[ty.Tuple[int, int], int] = collections.defaultdict(int)
    for i, gene in enumerate(kept):
        for member in gene["members"]:
            for j in by_member.get(member, ()):
                if kept[i]["assembly_id"] == nxt[j]["assembly_id"]:
                    shared[i, j] += 1

    def strength(item):
        (i, j), count = item
        union = len(kept[i]["members"] | nxt[j]["members"])
        return (-count, -count / union, kept[i]["internal_name"], j)

    heir_of: ty.Dict[int, int] = {}  # previous index -> next index
    prev_of: ty.Dict[int, int] = {}  # next index -> previous index
    for (i, j), count in sorted(shared.items(), key=strength):
        if i in heir_of or j in prev_of:
            continue
        if count / len(kept[i]["members"]) >= min_keep:
            heir_of[i] = j
            prev_of[j] = i

    successors: ty.Dict[int, ty.Set[int]] = collections.defaultdict(set)
    feeders: ty.Dict[int, ty.Set[int]] = collections.defaultdict(set)
    for (i, j), count in shared.items():
        if count / len(kept[i]["members"]) >= link_share:
            successors[i].add(j)
            feeders[j].add(i)
    ref = lambda event, i: {
        "event": event,
        "internal_name": kept[i]["internal_name"],
        "name": kept[i]["name"],
    }

    result = []
    for j, gene in enumerate(nxt):
        start, stop = true_extent(gene["members"])
        i = prev_of.get(j)
        others = sorted(k for k in feeders[j] if k != i)
        if i is None:
            base, _ = split_name(gene["name"])
            split_from = [k for k in others if len(successors[k]) > 1]
            merged_from = [k for k in others if len(successors[k]) <= 1]
            events["new"] += 1
            events["new_from_split"] += bool(split_from)
            version, first_release = 1, gene["first_release"]
            internal_name = gene["internal_name"]
        else:
            old = kept[i]
            base, old_version = split_name(old["name"])
            changed = (
                old["stored_members"] != gene["members"]
            )  # vs what is stored, not after retirement
            split_from, merged_from = [], others
            version = old_version + int(changed)
            first_release, internal_name = old["first_release"], old["internal_name"]
            events["changed" if changed else "continued"] += 1
            events["split"] += len(successors[i]) > 1
        events["merged_in"] += bool(merged_from)
        history = [ref("split_from", k) for k in split_from] + [
            ref("merged_from", k) for k in merged_from
        ]
        result.append(
            {
                "name": f"{base}.{version}",
                "internal_name": internal_name,
                "members": sorted(gene["members"]),
                "start": start,
                "stop": stop,
                "strand": gene["strand"],
                "chromosome": gene["chromosome"],
                "assembly_id": gene["assembly_id"],
                "version": version,
                "first_release": first_release,
                "last_release": gene["last_release"],
                "history": history,
            }
        )
    events["retired"] += sum(1 for i in range(len(kept)) if i not in heir_of)
    return result, events


def read_retired_urs(path) -> ty.Set[str]:
    """URS identifiers from the first column of a CSV / one-per-line file."""
    if path is None:
        return set()
    with open(path) as handle:
        cells = (line.split(",")[0].strip().strip('"') for line in handle)
        return {c.split("_", 1)[0] for c in cells if c.startswith("URS")}
