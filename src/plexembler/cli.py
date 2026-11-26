"""plexembler's command line: the commands the pipeline's steps run.

Each command that writes outputs takes ``-o/--output`` (or ``-o/--outdir``),
refuses to overwrite without ``--force``, names its run with ``-p/--prefix``
(the species, or the pool) and writes its log beside its outputs
(CONVENTIONS.md § 10, § 11). docs/cli.md documents every command and option.
"""

from __future__ import annotations

from pathlib import Path

import rich_click as click

from . import __version__
from ._log import logger, run_log
from .rename import SPECIES, rename_contigs, write_names

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


def _refuse_existing(paths: list[Path], force: bool) -> None:
    """Fail when an output exists, unless ``--force`` was given."""
    existing = [str(path) for path in paths if path.exists()]
    if existing and not force:
        raise click.ClickException(
            f"{', '.join(existing)} exists already: pass --force to overwrite it"
        )


def _species(_context: click.Context, _parameter: click.Parameter, value: str) -> str:
    """Check a species' ID as the option is read: a usage error otherwise."""
    from .rename import check_species

    try:
        return check_species(value)
    except ValueError as error:
        raise click.BadParameter(str(error)) from error


def _pool(_context: click.Context, _parameter: click.Parameter, value: str) -> str:
    """Check a pool's ID as the option is read: a usage error otherwise."""
    if not SPECIES.match(value):
        raise click.BadParameter(
            f"{value!r} is not a pool ID: a letter, then letters, digits and hyphens"
        )
    return value


#: The options every command that writes outputs has (CONVENTIONS.md § 10).
_FORCE = click.option(
    "--force", is_flag=True, help="Overwrite outputs that exist already."
)
_VERBOSE = click.option(
    "-v", "--verbose", is_flag=True, help="Debug messages, in the log and on stderr."
)
_QUIET = click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Only warnings and errors on stderr; the log keeps everything.",
)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__, "--version", prog_name="plexembler")
@click.command_panel("Contigs", commands=["rename", "select", "screen"])
@click.command_panel("Reads", commands=["digest"])
@click.command_panel("Pools", commands=["ani", "split"])
def cli() -> None:
    """The commands plexembler's pipeline runs on its data."""


@cli.command(context_settings=CONTEXT_SETTINGS)
@click.argument("fasta", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_species,
    help="The species' ID, which every new name starts with: a letter, then letters, digits and hyphens.",
)
@click.option(
    "-o",
    "--output",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="The renamed contigs, as FASTA; gzipped when the name ends in .gz.",
)
@click.option(
    "-n",
    "--names",
    type=click.Path(dir_okay=False, path_type=Path),
    help="The table of new and original names (default: the output's name, with .names.tsv for its extensions).",
)
@click.option("--force", is_flag=True, help="Overwrite outputs that exist already.")
@click.option(
    "-v", "--verbose", is_flag=True, help="Debug messages, in the log and on stderr."
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Only warnings and errors on stderr; the log keeps everything.",
)
def rename(
    fasta: Path,
    prefix: str,
    output: Path,
    names: Path | None,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Name every contig after its species: <species>_<contig>."""
    if names is None:
        stem = output.name.removesuffix(".gz")
        for extension in (".fasta", ".fa", ".fna"):
            stem = stem.removesuffix(extension)
        names = output.with_name(f"{stem}.names.tsv")
    _refuse_existing([output, names], force)
    parameters = {
        "fasta": fasta,
        "prefix": prefix,
        "output": output,
        "names": names,
        "force": force,
    }
    with run_log(output.parent, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            renamed = rename_contigs(fasta, output, prefix)
        except ValueError as error:
            raise click.ClickException(f"{fasta}: {error}") from error
        write_names(renamed, names)


@cli.command(name="select", context_settings=CONTEXT_SETTINGS)
@click.argument(
    "candidates", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "-s",
    "--step",
    type=click.Choice(["assembly", "polish"]),
    default="assembly",
    show_default=True,
    help="The step whose candidates these are, which sets the rules.",
)
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_species,
    help="The species' ID, which names the outputs and the log.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The folder of the kept sequence, the table of candidates and the log.",
)
@click.option(
    "-k",
    "--keep",
    help="A candidate to keep whatever its numbers: the samplesheet's choice.",
)
@click.option(
    "--busco-tolerance",
    type=float,
    default=0.5,
    show_default=True,
    help="Single-copy BUSCOs within this many points of the best are a tie.",
)
@click.option(
    "--qv-tolerance",
    type=float,
    default=1.0,
    show_default=True,
    help="QVs within this many units of the best are a tie.",
)
@_FORCE
@_VERBOSE
@_QUIET
def select_command(
    candidates: Path,
    step: str,
    prefix: str,
    outdir: Path,
    keep: str | None,
    busco_tolerance: float,
    qv_tolerance: float,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Keep one candidate of a step, by its checks, and say why."""
    from .select import output_paths, select

    _refuse_existing(list(output_paths(outdir, prefix, step).values()), force)
    parameters = {
        "candidates": candidates,
        "step": step,
        "prefix": prefix,
        "outdir": outdir,
        "keep": keep,
        "busco_tolerance": busco_tolerance,
        "qv_tolerance": qv_tolerance,
        "force": force,
    }
    with run_log(outdir, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            select(
                candidates, outdir, prefix, step, keep, busco_tolerance, qv_tolerance
            )
        except (OSError, ValueError) as error:
            raise click.ClickException(str(error)) from error


@cli.command(name="ani", context_settings=CONTEXT_SETTINGS)
@click.argument("skani", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_pool,
    help="The pool's ID, which names the outputs and the log.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The folder of the matrix, the pairs and the log.",
)
@click.option(
    "-s",
    "--species",
    multiple=True,
    help="A species of the pool, listed even when skani paired it with none; once per species.",
)
@click.option(
    "--warn-above",
    type=float,
    default=95.0,
    show_default=True,
    help="The identity, in percent, above which a pair is too close: a warning.",
)
@click.option(
    "--genus-above",
    type=float,
    default=85.0,
    show_default=True,
    help="The identity, in percent, from which a pair is of the same genus.",
)
@_FORCE
@_VERBOSE
@_QUIET
def ani_command(
    skani: Path,
    prefix: str,
    outdir: Path,
    species: tuple[str, ...],
    warn_above: float,
    genus_above: float,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """The identity between the species of a pool, from skani's pairs."""
    from .ani import ani_matrix, read_skani

    names = ("ani.tsv", "ani_pairs.tsv", "ani_mqc.tsv")
    _refuse_existing([outdir / f"{prefix}.{name}" for name in names], force)
    parameters = {
        "skani": skani,
        "prefix": prefix,
        "outdir": outdir,
        "species": ",".join(species),
        "warn_above": warn_above,
        "genus_above": genus_above,
        "force": force,
    }
    with run_log(outdir, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            pairs = read_skani(skani)
        except ValueError as error:
            raise click.ClickException(str(error)) from error
        ani_matrix(prefix, pairs, list(species), outdir, warn_above, genus_above)


@cli.command(name="split", context_settings=CONTEXT_SETTINGS)
@click.argument("sam", type=click.Path(allow_dash=True, dir_okay=False, path_type=Path))
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_pool,
    help="The pool's ID, which names the tables and the log.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The folder of each species' FASTQ files, the tables and the log.",
)
@click.option(
    "-l",
    "--library",
    required=True,
    help="The library's ID, which names each species' FASTQ files.",
)
@click.option(
    "-t",
    "--type",
    "kind",
    type=click.Choice(["wgs", "hic"]),
    required=True,
    help="The library's type: wgs, or hic, which drops duplicates, counts the contacts between contigs and estimates epsilon.",
)
@click.option(
    "--min-mapq",
    type=click.IntRange(0, 60),
    default=20,
    show_default=True,
    help="A mate's lowest mapping quality: below, it is ambiguous.",
)
@click.option(
    "--min-identity",
    type=click.FloatRange(0, 1),
    default=0.98,
    show_default=True,
    help="A mate's lowest identity to its contig: below, it is divergent, a read of another species.",
)
@click.option(
    "--holdout",
    type=click.FloatRange(0, 1),
    default=0.1,
    show_default=True,
    help="The share of the library's pairs kept apart, for the QV guard.",
)
@click.option(
    "--dedup/--no-dedup",
    default=True,
    show_default=True,
    help="Drop a Hi-C library's duplicate pairs.",
)
@click.option(
    "--compress-level",
    type=click.IntRange(1, 9),
    default=1,
    show_default=True,
    help="The gzip level of the FASTQ files.",
)
@_FORCE
@_VERBOSE
@_QUIET
def split_command(
    sam: Path,
    prefix: str,
    outdir: Path,
    library: str,
    kind: str,
    min_mapq: int,
    min_identity: float,
    holdout: float,
    dedup: bool,
    compress_level: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Split a pooled library between its pool's species, from its alignments."""
    from .split import open_sam, split_pairs, summary, write_contacts, write_outputs

    _refuse_existing([outdir / f"{prefix}.{library}.split.json"], force)
    parameters = {
        "sam": sam,
        "prefix": prefix,
        "outdir": outdir,
        "library": library,
        "type": kind,
        "min_mapq": min_mapq,
        "min_identity": min_identity,
        "holdout": holdout,
        "dedup": dedup,
        "compress_level": compress_level,
        "force": force,
    }
    with run_log(outdir, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            with open_sam(sam) as handle:
                species, tally = split_pairs(
                    handle,
                    outdir,
                    library,
                    kind,
                    min_mapq,
                    min_identity,
                    holdout,
                    dedup,
                    compress_level,
                )
        except (OSError, ValueError) as error:
            raise click.ClickException(str(error)) from error
        numbers = summary(species, tally, kind)
        write_outputs(numbers, outdir, prefix, library)
        if kind == "hic":
            write_contacts(tally, outdir, prefix, library)
        for name, row in numbers["species"].items():
            logger.info(
                "%s: %d pairs (%.1f %%)", name, row["pairs"], 100 * row["share"]
            )
        logger.info(
            "%s: %d pairs; %d ambiguous, %d divergent, %d between species, %d unmapped",
            library,
            numbers["pairs"],
            numbers["ambiguous"],
            numbers["divergent"],
            numbers["between_species"],
            numbers["unmapped"],
        )
        if numbers.get("epsilon") is not None:
            logger.info("Random ligation rate (epsilon): %.4f", numbers["epsilon"])


@cli.command(name="digest", context_settings=CONTEXT_SETTINGS)
@click.argument(
    "fastqs",
    nargs=-1,
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "-e",
    "--enzymes",
    required=True,
    help="The library's restriction enzymes, comma-separated (DpnII,HinfI), or arima.",
)
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_species,
    help="The species' ID, which names the outputs and the log.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The folder of the pieces, their numbers and the log.",
)
@click.option(
    "--min-piece",
    type=click.IntRange(1),
    default=40,
    show_default=True,
    help="The shortest piece kept, in bases.",
)
@_FORCE
@_VERBOSE
@_QUIET
def digest_command(
    fastqs: tuple[Path, ...],
    enzymes: str,
    prefix: str,
    outdir: Path,
    min_piece: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Cut Hi-C reads at their ligation junctions, into single reads to polish with."""
    from .digest import digest, parse_enzymes, write_outputs

    output = outdir / f"{prefix}.digested.fastq.gz"
    _refuse_existing([output, outdir / f"{prefix}.digest.json"], force)
    try:
        names = parse_enzymes(enzymes)
    except ValueError as error:
        raise click.BadParameter(str(error), param_hint="--enzymes") from error
    parameters = {
        "fastqs": " ".join(str(f) for f in fastqs),
        "enzymes": ",".join(names),
        "prefix": prefix,
        "outdir": outdir,
        "min_piece": min_piece,
        "force": force,
    }
    with run_log(outdir, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            counts = digest(list(fastqs), output, names, min_piece)
        except ValueError as error:
            raise click.ClickException(str(error)) from error
        write_outputs(counts, outdir, prefix)
        reads = sum(counts.reads.values())
        logger.info(
            "%d reads: %d without a junction, %d with one, %d with two or more",
            reads,
            counts.reads[0],
            counts.reads[1],
            counts.reads[2],
        )
        logger.info(
            "%d pieces of %d bases or more, %d bases",
            counts.pieces,
            min_piece,
            counts.bases,
        )


def _species_paf(
    _context: click.Context, _parameter: click.Parameter, values: tuple[str, ...]
) -> dict[str, Path]:
    """``SPECIES=PAF`` pairs, as the screen's reads take them."""
    found = {}
    for value in values:
        species, sep, path = value.partition("=")
        if not sep or not species or not path:
            raise click.BadParameter(f"{value!r} is not SPECIES=PAF")
        if not Path(path).exists():
            raise click.BadParameter(f"{path} does not exist")
        found[species] = Path(path)
    return found


@cli.command(name="screen", context_settings=CONTEXT_SETTINGS)
@click.argument("contigs", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "-r",
    "--reads",
    multiple=True,
    required=True,
    callback=_species_paf,
    help="A species' nanopore reads aligned to the pool's contigs, SPECIES=PAF; once per species.",
)
@click.option(
    "-i",
    "--identity",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The pool's contigs aligned to each other (minimap2 -x asm10 -X), PAF.",
)
@click.option(
    "-c",
    "--contacts",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="A Hi-C library's contacts between contigs (plexembler split); once per library.",
)
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_pool,
    help="The pool's ID, which names the table and the log.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The folder of each species' screened contigs, the table and the log.",
)
@click.option(
    "--organelle-depth",
    type=float,
    default=3.0,
    show_default=True,
    help="An organelle's depth, in times its species' contigs'.",
)
@click.option(
    "--organelle-gc",
    type=float,
    default=0.3,
    show_default=True,
    help="An organelle's highest GC share.",
)
@click.option(
    "--duplicate-share",
    type=float,
    default=0.9,
    show_default=True,
    help="The share of a contig found in another species' contigs that makes it theirs.",
)
@click.option(
    "--own-depth",
    type=float,
    default=0.2,
    show_default=True,
    help="Its own species' reads' depth, in times their usual depth, under which they are thin on it.",
)
@click.option(
    "--other-depth",
    type=float,
    default=0.5,
    show_default=True,
    help="Another species' reads' depth, in times their usual depth, at which they cover it.",
)
@click.option(
    "--hic-ratio",
    type=float,
    default=2.0,
    show_default=True,
    help="How many times its Hi-C contacts with another species outnumber its own species'.",
)
@click.option(
    "--min-contacts",
    type=int,
    default=20,
    show_default=True,
    help="The fewest Hi-C contacts with another species that move a contig.",
)
@_FORCE
@_VERBOSE
@_QUIET
def screen_command(
    contigs: Path,
    reads: dict[str, Path],
    identity: Path | None,
    contacts: tuple[Path, ...],
    prefix: str,
    outdir: Path,
    organelle_depth: float,
    organelle_gc: float,
    duplicate_share: float,
    own_depth: float,
    other_depth: float,
    hic_ratio: float,
    min_contacts: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Keep, move or drop each contig of a pool, and name its organelles."""
    from .screen import Thresholds, screen

    _refuse_existing([outdir / f"{prefix}.screen.tsv"], force)
    thresholds = Thresholds(
        organelle_depth=organelle_depth,
        organelle_gc=organelle_gc,
        duplicate_share=duplicate_share,
        own_depth=own_depth,
        other_depth=other_depth,
        hic_ratio=hic_ratio,
        min_contacts=min_contacts,
    )
    parameters = {
        "contigs": contigs,
        "reads": ",".join(f"{k}={v}" for k, v in reads.items()),
        "identity": identity,
        "contacts": ",".join(str(c) for c in contacts),
        "prefix": prefix,
        "outdir": outdir,
        **{key: getattr(thresholds, key) for key in thresholds.__dataclass_fields__},
        "force": force,
    }
    with run_log(outdir, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            found = screen(
                contigs, reads, identity, list(contacts), outdir, prefix, thresholds
            )
        except (OSError, ValueError, KeyError) as error:
            raise click.ClickException(str(error)) from error
        verdicts = {}
        for contig in found.values():
            verdicts[contig.verdict] = verdicts.get(contig.verdict, 0) + 1
        logger.info(
            "%d contigs: %s",
            len(found),
            ", ".join(f"{n} {v}" for v, n in sorted(verdicts.items())),
        )
