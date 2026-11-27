"""A species' closest known genome: sourmash gather's best match, and its sequence.

``plexembler closest`` reads the matches of ``sourmash gather`` between a
species' genome and a database of genomes (one of sourmash's prepared
databases of NCBI's genomes, docs/usage.md), and writes them best first: the
first is the species' closest known genome. With the database's taxonomy, each
match has its lineage.

The first match's sequence is then downloaded from NCBI Datasets, by the
assembly accession its name starts with (``GCA_000146045.2 Saccharomyces
cerevisiae S288C``): the reference of the species' dotplot
(``plexembler dotplot``). The download is the step's only use of the network,
and it may fail: the node has no network, NCBI does not know the accession, or
the database's names are not NCBI's. The command then warns, writes its table
without the sequence and ends without an error: the species has no dotplot,
and the rest of the run goes on.
"""

from __future__ import annotations

import csv
import gzip
import http.client
import io
import logging
import re
import shutil
import ssl
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
import zlib
from dataclasses import dataclass, replace
from pathlib import Path
from typing import IO

from . import __version__
from .rename import open_text

logger = logging.getLogger(__name__)

#: NCBI Datasets' API, version 2, which the closest genome is downloaded from.
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"

#: An NCBI assembly accession: GenBank's (GCA_) or RefSeq's (GCF_), and its version.
ACCESSION = re.compile(r"\bGC[AF]_\d{9}(?:\.\d+)?\b")

#: The ranks of a sourmash taxonomy, in order; a file has some of them.
RANKS = (
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
)

#: What became of the closest genome's sequence: the table's ``sequence`` column.
DOWNLOADED = "downloaded"
OFFLINE = "offline"
NOT_FOUND = "not_found"
FAILED = "failed"
NOT_NCBI = "not_ncbi"
NOT_DOWNLOADED = "not_downloaded"

#: HTTP answers worth asking again: too many requests, and the server's errors.
_TRANSIENT = frozenset({408, 429, 500, 502, 503, 504})

#: What a connection that fails on the way raises: the network's errors, not a
#: disk's, which are OSErrors too but no reason to say the node is offline.
_NETWORK_ERRORS = (
    TimeoutError,
    ConnectionError,
    ssl.SSLError,
    http.client.HTTPException,
)

#: How much of NCBI's answer is read at once, in bytes.
_CHUNK = 1 << 20


@dataclass(frozen=True)
class Match:
    """A genome of the database, as sourmash gather matched it to a species'."""

    rank: int
    name: str
    accession: str | None
    ani: float | None
    species_containment: float | None
    match_containment: float | None
    shared_bp: int
    lineage: str | None = None

    @property
    def title(self) -> str:
        """The genome's name without its accession: its organism and strain."""
        if self.accession and self.name.startswith(self.accession):
            return self.name[len(self.accession) :].strip() or self.name
        return self.name


@dataclass(frozen=True)
class Fetched:
    """What became of a download: its status, and why, in words."""

    status: str
    detail: str


def _percent(value: str | None) -> float | None:
    """A fraction of gather's CSV as a percentage; None when it is empty."""
    if value in (None, ""):
        return None
    try:
        return float(value) * 100
    except ValueError:
        return None


def read_gather(path: Path) -> list[Match]:
    """sourmash gather's matches, best first.

    An empty file is no match: gather writes one with ``--create-empty-results``
    when no genome of the database shares enough of the species' k-mers.

    Raises
    ------
    ValueError
        When the file has no column ``name`` or ``intersect_bp``.
    """
    with open_text(path) as handle:
        text = handle.read()
    if not text.strip():
        return []
    reader = csv.DictReader(io.StringIO(text))
    fields = set(reader.fieldnames or [])
    name_column = "name" if "name" in fields else "match_name"
    missing = {name_column, "intersect_bp"} - fields
    if missing:
        raise ValueError(f"{path}: no column {', '.join(sorted(missing))}")
    rows = list(reader)
    if "gather_result_rank" in fields:
        rows.sort(key=lambda row: int(row["gather_result_rank"] or 0))
    matches = []
    for rank, row in enumerate(rows, start=1):
        name = row[name_column].strip()
        found = ACCESSION.search(name)
        matches.append(
            Match(
                rank=rank,
                name=name,
                accession=found.group(0) if found else None,
                ani=_percent(row.get("max_containment_ani")),
                species_containment=_percent(row.get("f_orig_query")),
                match_containment=_percent(row.get("f_match")),
                shared_bp=int(float(row["intersect_bp"] or 0)),
            )
        )
    return matches


def _unversioned(accession: str) -> str:
    """An accession without its version: sourmash's taxonomy matches them so."""
    return accession.split(".", 1)[0]


def read_lineages(path: Path, accessions: set[str]) -> dict[str, str]:
    """The lineages of the given accessions, from a sourmash taxonomy (CSV).

    The file names its genomes in an ``ident`` column (or ``identifiers``, or
    ``accession``, as ``sourmash tax`` reads them) and gives their lineages in
    the ranks' columns (``superkingdom`` … ``strain``) or in a ``lineage``
    column. Accessions are matched without their version.

    Returns
    -------
    dict
        Each accession found, versioned as asked, and its lineage, its ranks
        joined by ``;``.

    Raises
    ------
    ValueError
        When the file names no genome or gives no lineage.
    """
    wanted = {_unversioned(accession): accession for accession in accessions}
    found: dict[str, str] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        ident = next(
            (f for f in ("ident", "identifiers", "accession") if f in fields), None
        )
        ranks = [rank for rank in RANKS if rank in fields]
        if ident is None or not (ranks or "lineage" in fields):
            raise ValueError(
                f"{path}: not a sourmash taxonomy: no column ident, identifiers or "
                "accession, or no ranks"
            )
        for row in reader:
            key = _unversioned((row.get(ident) or "").strip())
            if key not in wanted:
                continue
            if ranks:
                values = [(row.get(rank) or "").strip() for rank in ranks]
            else:
                values = [
                    part.strip() for part in (row.get("lineage") or "").split(";")
                ]
            found[wanted[key]] = ";".join(value for value in values if value)
    return found


def download_url(api: str, accession: str) -> str:
    """The address of an assembly's genome sequence, in NCBI Datasets' API."""
    return (
        f"{api.rstrip('/')}/genome/accession/{accession}/download"
        "?include_annotation_type=GENOME_FASTA"
    )


def unpack(archive: IO[bytes], accession: str, output: Path) -> Fetched:
    """Write the genome of NCBI Datasets' package to output, gzipped.

    NCBI answers an accession it does not know with a broken package, or one
    without a sequence: both are ``not_found``.
    """
    try:
        bundle = zipfile.ZipFile(archive)
    except zipfile.BadZipFile:
        return Fetched(NOT_FOUND, f"NCBI sent no genome for {accession}")
    with bundle:
        members = sorted(name for name in bundle.namelist() if name.endswith(".fna"))
        if not members:
            return Fetched(NOT_FOUND, f"NCBI's package for {accession} has no sequence")
        partial = output.with_name(f".{output.name}.partial")
        try:
            with (
                bundle.open(members[0]) as source,
                gzip.open(partial, "wb", compresslevel=6) as sink,
            ):
                shutil.copyfileobj(source, sink)
        except (zipfile.BadZipFile, EOFError, zlib.error) as error:
            partial.unlink(missing_ok=True)
            return Fetched(
                FAILED, f"NCBI's package for {accession} is broken ({error})"
            )
    partial.replace(output)
    return Fetched(
        DOWNLOADED, f"{accession}, from NCBI ({members[0].rsplit('/', 1)[-1]})"
    )


def _receive(
    request: urllib.request.Request, archive: IO[bytes], timeout: float
) -> None:
    """Write NCBI's answer to archive. The network's errors are raised as they
    come; a disk's are raised too, but are not the network's (see fetch_genome)."""
    with urllib.request.urlopen(request, timeout=timeout) as response:
        while chunk := response.read(_CHUNK):
            archive.write(chunk)


def fetch_genome(
    accession: str,
    output: Path,
    api: str = DATASETS,
    timeout: float = 60.0,
    attempts: int = 3,
    wait: float = 2.0,
) -> Fetched:
    """Download an assembly's genome sequence from NCBI Datasets to output.

    Never raises for the network: a node without network, a timeout, or an
    answer that is not the genome each give a status, which the caller logs.

    Parameters
    ----------
    accession
        The assembly's accession (``GCA_000146045.2``).
    output
        The genome, gzipped FASTA, written only when it was downloaded whole.
    api
        NCBI Datasets' API.
    timeout
        How long to wait for NCBI, in seconds, at each step of each attempt.
    attempts
        How many times to try, when NCBI cannot be reached or is busy.
    wait
        How long to wait before the second attempt, in seconds; twice as long
        before the third, and so on.

    Returns
    -------
    Fetched
        ``downloaded``; ``offline`` when NCBI could not be reached;
        ``not_found`` when it does not know the accession; ``failed`` for
        any other answer.
    """
    request = urllib.request.Request(
        download_url(api, accession),
        headers={
            "User-Agent": f"plexembler/{__version__}",
            "Accept": "application/zip",
        },
    )
    result = Fetched(FAILED, "no attempt")
    for attempt in range(1, attempts + 1):
        # Beside the output, not in a node's /tmp, which a genome may not fit in.
        with tempfile.TemporaryFile(dir=output.parent) as archive:
            try:
                _receive(request, archive, timeout)
            except urllib.error.HTTPError as error:
                if error.code in (400, 404):
                    return Fetched(
                        NOT_FOUND,
                        f"NCBI has no assembly {accession} (HTTP {error.code})",
                    )
                result = Fetched(
                    FAILED, f"NCBI answered HTTP {error.code} ({error.reason})"
                )
                if error.code not in _TRANSIENT:
                    return result
            except urllib.error.URLError as error:
                result = Fetched(OFFLINE, f"NCBI could not be reached ({error.reason})")
            except _NETWORK_ERRORS as error:
                # A timeout, or a connection lost on the way.
                result = Fetched(OFFLINE, f"NCBI could not be reached ({error!r})")
            else:
                archive.seek(0)
                return unpack(archive, accession, output)
        if attempt < attempts:
            logger.debug("Attempt %d of %d: %s", attempt, attempts, result.detail)
            time.sleep(wait * 2 ** (attempt - 1))
    return result


def _number(value: float | None) -> str:
    return f"{value:.2f}" if value is not None else "NA"


def write_outputs(
    matches: list[Match], fetched: Fetched | None, outdir: Path, prefix: str
) -> dict[str, Path]:
    """Write the species' matches and MultiQC's table.

    Returns
    -------
    dict
        The paths written: ``table`` and ``multiqc``.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "table": outdir / f"{prefix}.closest.tsv",
        "multiqc": outdir / f"{prefix}.closest_mqc.tsv",
    }
    columns = (
        "rank",
        "accession",
        "name",
        "ani",
        "species_containment",
        "match_containment",
        "shared_bp",
        "lineage",
        "sequence",
    )
    with paths["table"].open("w", encoding="utf-8") as sink:
        sink.write("\t".join(columns) + "\n")
        for match in matches:
            sequence = fetched.status if match.rank == 1 and fetched else "NA"
            cells = [
                str(match.rank),
                match.accession or "NA",
                match.name,
                _number(match.ani),
                _number(match.species_containment),
                _number(match.match_containment),
                str(match.shared_bp),
                match.lineage or "NA",
                sequence,
            ]
            sink.write("\t".join(cells) + "\n")

    # One table of every species: each species' file has the same section ID.
    lines = [
        "# id: 'plexembler_closest'",
        "# section_name: 'Closest known genome'",
        "# description: 'Each species'' closest genome in the sourmash database: "
        "its ANI, estimated from the k-mers they share, and whether its sequence "
        "was downloaded from NCBI for the dotplot (plexembler closest).'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'plexembler_closest_table'",
        "#     namespace: 'plexembler'",
        "Sample\tClosest genome\tAccession\tANI (%)\tSpecies k-mers in it (%)"
        "\tIts k-mers in the species (%)\tLineage\tSequence",
    ]
    if matches:
        best = matches[0]
        cells = [
            prefix,
            best.title,
            best.accession or "NA",
            _number(best.ani),
            _number(best.species_containment),
            _number(best.match_containment),
            best.lineage or "NA",
            fetched.status if fetched else "NA",
        ]
    else:
        cells = [prefix, "no match", "NA", "NA", "NA", "NA", "NA", "NA"]
    lines.append("\t".join(cells))
    paths["multiqc"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return paths


def closest(
    gather: Path,
    prefix: str,
    outdir: Path,
    lineages: Path | None = None,
    top: int = 5,
    download: bool = True,
    api: str = DATASETS,
    timeout: float = 60.0,
    attempts: int = 3,
) -> tuple[list[Match], Fetched | None]:
    """A species' closest genomes, and the first one's sequence when it can be had.

    Writes ``<prefix>.closest.tsv``, ``<prefix>.closest_mqc.tsv`` and, when
    the sequence was downloaded, ``<prefix>.closest.fa.gz`` in outdir. A
    sequence that cannot be had is a warning, never an error.

    Raises
    ------
    ValueError
        When gather's file or the taxonomy cannot be read.
    """
    matches = read_gather(gather)[:top]
    if lineages is not None and matches:
        found = read_lineages(lineages, {m.accession for m in matches if m.accession})
        matches = [
            replace(match, lineage=found.get(match.accession or ""))
            for match in matches
        ]
    outdir.mkdir(parents=True, exist_ok=True)
    fetched: Fetched | None = None
    if not matches:
        logger.warning(
            "%s: no genome of the database shares enough of its k-mers: the "
            "species is far from every genome there, and has no dotplot",
            prefix,
        )
    else:
        best = matches[0]
        logger.info(
            "%s: its closest genome is %s, at %s %% ANI",
            prefix,
            best.name,
            _number(best.ani),
        )
        if best.lineage:
            logger.info("%s: %s", prefix, best.lineage)
        if best.accession is None:
            fetched = Fetched(NOT_NCBI, f"{best.name!r} has no NCBI assembly accession")
        elif not download:
            fetched = Fetched(NOT_DOWNLOADED, "--no-download")
        else:
            fetched = fetch_genome(
                best.accession,
                outdir / f"{prefix}.closest.fa.gz",
                api=api,
                timeout=timeout,
                attempts=attempts,
            )
        if fetched.status == DOWNLOADED:
            logger.info("%s: downloaded %s", prefix, fetched.detail)
        elif fetched.status == NOT_DOWNLOADED:
            logger.info("%s: its sequence not downloaded (--no-download)", prefix)
        else:
            logger.warning(
                "%s: its closest genome's sequence was not downloaded: %s. "
                "%s has no dotplot; the rest of the run goes on",
                prefix,
                fetched.detail,
                prefix,
            )
    write_outputs(matches, fetched, outdir, prefix)
    return matches, fetched
