"""The results, as chromosomix opens them: a store of every species' files.

``plexembler browser`` writes ``chromosomix.store.json`` at the top of the
results (docs/plan.md, M6): chromosomix's store (its Contract 1, version 1),
a row per species' genome and Hi-C map, each with its tracks (genes, RNA-seq
coverage, nanopore depth) placed, every file named by its path in the results.
Laid out as cards, one per species, it opens as it is in chromosomix: its app
with the results folder connected, or chromosomix.app with ``?store=<url>``
once the results are on a bucket (``--base``).

Each file is given as ``SPECIES=PATH``, its path relative to the results, and
known by its name: a ``.mcool`` is the species' map, a FASTA its genome, a
GFF3, GTF or bigWig one of its tracks.
"""

from __future__ import annotations

import json
import logging
import re
import time
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

STORE_KIND = "chromosomix.store"
STORE_VERSION = 1
STORE_NAME = "chromosomix.store.json"

_GENOME = re.compile(r"\.(fa|fasta|fna)(\.gz)?$")
_TRACK = re.compile(r"\.(gff3|gff|gtf)(\.gz)?$|\.(bigwig|bw)$", re.IGNORECASE)


@dataclass
class Species:
    """A species' files: its genome, its map, its tracks."""

    name: str
    pool: str = ""
    genome: str = ""
    map: str = ""
    tracks: list[str] = field(default_factory=list)


def track_label(path: str, species: str) -> str:
    """A track's label, from its file's name: ``genes``, ``rna-ypd-1 forward``, …"""
    name = path.rsplit("/", 1)[-1]
    if _GENOME.search(name) is None and re.search(r"\.(gff3|gff|gtf)(\.gz)?$", name):
        return "genes"
    stem = re.sub(r"\.(bigwig|bw)$", "", name, flags=re.IGNORECASE)
    if stem.startswith(f"{species}."):
        stem = stem[len(species) + 1 :]
    return stem.replace(".", " ")


def gather(pools: dict[str, str], files: list[tuple[str, str]]) -> list[Species]:
    """Each species' files, sorted by what they are.

    Raises
    ------
    ValueError
        When a file names a species no ``--species`` gives, or is none of a
        map, a genome or a track.
    """
    found = {name: Species(name, pool) for name, pool in pools.items()}
    for species, path in files:
        if species not in found:
            raise ValueError(f"{path}: {species} is not one of the species given")
        entry = found[species]
        if path.endswith(".mcool"):
            entry.map = path
        elif _GENOME.search(path):
            entry.genome = path
        elif _TRACK.search(path):
            entry.tracks.append(path)
        else:
            raise ValueError(f"{path}: neither a map, a genome nor a track")
    for entry in found.values():
        entry.tracks.sort(key=lambda p: (track_label(p, entry.name) != "genes", p))
    return sorted(found.values(), key=lambda s: s.name)


def store(
    species: list[Species], name: str, base: str | None = None, description: str = ""
) -> dict:
    """The store: a row per species' genome and map, its tracks on each."""
    keyed: list[tuple[tuple[str, str], dict]] = []
    for entry in species:
        tracks = [
            {"label": track_label(track, entry.name), "url": track}
            for track in entry.tracks
        ]
        rows = [
            (entry.genome, "genome", "Genome"),
            (entry.map, "map", "Hi-C map"),
        ]
        for uri, kind, kind_label in rows:
            if not uri:
                continue
            fields = {
                "assembly": entry.name,
                "species": entry.name,
                "type": kind_label,
                "tracks": str(len(tracks)),
            }
            if entry.pool:
                fields["pool"] = entry.pool
            row = {
                "label": f"{entry.name} {kind_label.lower()}",
                "uri": uri,
                "kind": kind,
                "fields": fields,
            }
            if tracks:
                row["tracks"] = tracks
            keyed.append(((entry.name, str(row["label"])), row))
    # Ordered by assembly, then label, as chromosomix's store wants them.
    entries = [row for _, row in sorted(keyed, key=lambda pair: pair[0])]
    document: dict = {
        "kind": STORE_KIND,
        "version": STORE_VERSION,
        "name": name,
        "generatedAt": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if base:
        document["base"] = base if base.endswith("/") else base + "/"
    if description:
        document["description"] = description
    document.update(
        {
            "columns": [
                {"key": "species", "label": "Species"},
                {"key": "pool", "label": "Pool"},
                {"key": "type", "label": "Type"},
                {"key": "tracks", "label": "Tracks"},
            ],
            "layout": "cards",
            "group": "species",
            "headline": ["pool", "type"],
            "entries": entries,
        }
    )
    return document


def write(document: dict, outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / STORE_NAME
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    logger.info(
        "%s: %d rows, %d species",
        path,
        len(document["entries"]),
        len({row["fields"]["species"] for row in document["entries"]}),
    )
    return path
