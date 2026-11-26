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
from ._log import run_log
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
@click.command_panel("Contigs", commands=["rename", "select"])
@click.command_panel("Pools", commands=["ani"])
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
    type=click.Choice(["assembly"]),
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
