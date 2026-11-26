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
from .rename import rename_contigs, write_names

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


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__, "--version", prog_name="plexembler")
@click.command_panel("Contigs", commands=["rename"])
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
