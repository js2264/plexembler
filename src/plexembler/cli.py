"""plexembler's command line: the pipeline, and the commands its steps run.

``plexembler run`` runs the pipeline with Nextflow. Each command that writes
outputs takes ``-o/--output`` (or ``-o/--outdir``), refuses to overwrite
without ``--force``, names its run with ``-p/--prefix`` (the species, or the
pool) and writes its log beside its outputs (CONVENTIONS.md § 10, § 11).
docs/cli.md documents every command and option.
"""

from __future__ import annotations

import shlex
from pathlib import Path

import rich_click as click

from . import __version__
from ._log import logger, run_log, stderr_log
from .rename import SPECIES, rename_contigs, write_names
from .run import REPOSITORY, split_arguments

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
@click.command_panel("Pipeline", commands=["run"])
@click.command_panel("Contigs", commands=["rename", "select", "screen"])
@click.command_panel(
    "Genomes", commands=["finalize", "junctions", "genes", "expression"]
)
@click.command_panel("Reports", commands=["report", "browser"])
@click.command_panel("Reads", commands=["digest"])
@click.command_panel("Pools", commands=["ani", "split"])
def cli() -> None:
    """plexembler's pipeline, and the commands it runs on its data."""


#: What `plexembler run` says when it finds no Nextflow to run.
NO_NEXTFLOW = (
    "Nextflow was not found, beside plexembler or on PATH: install plexembler "
    "with its dependencies (uvx plexembler run …, or uv tool install plexembler), "
    "which bring Nextflow's launcher, or install Nextflow 26.04.0 or newer "
    "(https://js2264.github.io/plexembler/installation/)"
)


class _NextflowCommand(click.RichCommand):
    """A command whose arguments, but its own options, go to Nextflow as they are.

    click would read Nextflow's single-dash options (``-profile``,
    ``-with-docker``) as bundles of its short options, ``-h`` among them: the
    command's own options are taken out by their exact words, and the rest is
    handed to click after ``--``, untouched.
    """

    def parse_args(self, ctx: click.Context, args: list[str]) -> list[str]:
        if not args:
            return super().parse_args(ctx, args)
        flags: set[str] = set()
        valued: set[str] = set()
        for param in self.get_params(ctx):
            if isinstance(param, click.Option):
                names = [*param.opts, *param.secondary_opts]
                (flags if param.is_flag else valued).update(names)
        try:
            own, nextflow_args = split_arguments(args, flags, valued)
        except ValueError as error:
            raise click.UsageError(str(error), ctx=ctx) from error
        return super().parse_args(ctx, [*own, "--", *nextflow_args])


@cli.command(
    name="run",
    cls=_NextflowCommand,
    no_args_is_help=True,
    context_settings=CONTEXT_SETTINGS,
    epilog=(
        "Every other argument goes to Nextflow as it is: its options (-profile, "
        "-resume, -c, -r) and the pipeline's parameters (--input, --libraries, "
        "--outdir). After --, every argument does: plexembler run -- --help "
        "lists the pipeline's parameters."
    ),
)
@click.argument("nextflow_args", nargs=-1, type=click.UNPROCESSED)
@click.option(
    "--pipeline",
    help=(
        "The pipeline: a folder (a clone), run as it is, or a repository, which "
        f"Nextflow pulls, run at v{__version__} unless -r names a revision "
        f"[default: the folder or the commit plexembler was installed from, "
        f"else {REPOSITORY}]."
    ),
)
@click.option(
    "--dry-run", is_flag=True, help="Print the Nextflow command, and run nothing."
)
def run_command(
    nextflow_args: tuple[str, ...], pipeline: str | None, dry_run: bool
) -> None:
    """Run the pipeline with Nextflow, at this version of plexembler."""
    from .run import find_nextflow, launch, nextflow_command, pipeline_to_run

    with stderr_log():
        name, revision = pipeline_to_run(pipeline, nextflow_args)
        nextflow = find_nextflow()
        if nextflow is None and not dry_run:
            raise click.ClickException(NO_NEXTFLOW)
        command = nextflow_command(
            nextflow or "nextflow", name, revision, nextflow_args
        )
        if dry_run:
            if nextflow is None:
                logger.warning(NO_NEXTFLOW)
            click.echo(shlex.join(command))
            return
        logger.info("Running %s", shlex.join(command))
        try:
            launch(command)
        except OSError as error:
            raise click.ClickException(f"{command[0]}: {error}") from error


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
    help="The folder of each species' FASTQ (or SAM) files, the tables and the log.",
)
@click.option(
    "-l",
    "--library",
    required=True,
    help="The library's ID, which names each species' FASTQ (or SAM) files.",
)
@click.option(
    "-t",
    "--type",
    "kind",
    type=click.Choice(["wgs", "hic", "rnaseq"]),
    required=True,
    help="The library's type: wgs; hic, which drops duplicates, counts the contacts "
    "between contigs and estimates epsilon; or rnaseq, STAR's alignments, written "
    "as each species' SAM.",
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
    default=2.0,
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
@click.option(
    "--hic-own-depth",
    type=float,
    default=0.5,
    show_default=True,
    help="Its own species' depth on a contig, as a share of their depth, under "
    "which Hi-C contacts can move it.",
)
@click.option(
    "--depth-ratio",
    type=float,
    default=2.0,
    show_default=True,
    help="Another species' depth on a duplicate, over its own species', each against "
    "their depth, from which the duplicate goes.",
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
    hic_own_depth: float,
    depth_ratio: float,
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
        hic_own_depth=hic_own_depth,
        depth_ratio=depth_ratio,
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


@cli.command(name="finalize", context_settings=CONTEXT_SETTINGS)
@click.argument("fasta", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_species,
    help="The species' ID, which starts every name and names the outputs.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The folder of the genome, its sizes, names, AGP and the log.",
)
@click.option(
    "-a",
    "--agp",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The scaffolder's AGP, written again with the new names.",
)
@click.option(
    "-n",
    "--n-chromosomes",
    type=click.IntRange(1),
    help="The species' number of nuclear chromosomes: its longest scaffolds, named chromosomes.",
)
@click.option(
    "--min-chromosome",
    type=click.IntRange(1),
    default=100_000,
    show_default=True,
    help="Without --n-chromosomes, the shortest scaffold named a chromosome, in bases.",
)
@_FORCE
@_VERBOSE
@_QUIET
def finalize_command(
    fasta: Path,
    prefix: str,
    outdir: Path,
    agp: Path | None,
    n_chromosomes: int | None,
    min_chromosome: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Name a species' scaffolds: its chromosomes by size, its organelle, the unplaced."""
    from .finalize import finalize

    _refuse_existing(
        [outdir / f"{prefix}.fa.gz", outdir / f"{prefix}.chrom.sizes"], force
    )
    parameters = {
        "fasta": fasta,
        "prefix": prefix,
        "outdir": outdir,
        "agp": agp,
        "n_chromosomes": n_chromosomes,
        "min_chromosome": min_chromosome,
        "force": force,
    }
    with run_log(outdir, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            finalize(fasta, prefix, outdir, agp, n_chromosomes, min_chromosome)
        except ValueError as error:
            raise click.ClickException(str(error)) from error


@cli.command(name="junctions", context_settings=CONTEXT_SETTINGS)
@click.argument("pairs", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "-a",
    "--agp",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The genome's AGP, in the names of the pairs' chromosomes.",
)
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_species,
    help="The species' ID, which names the outputs.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The folder of the table of joins, its MultiQC table and the log.",
)
@click.option(
    "-w",
    "--window",
    type=click.IntRange(1),
    default=20_000,
    show_default=True,
    help="How far from a join, on each side, a pair's mates may lie, in bases.",
)
@click.option(
    "--min-ratio",
    type=click.FloatRange(0, min_open=True),
    default=0.2,
    show_default=True,
    help="The pairs across a join over those expected, below which it is unsupported.",
)
@click.option(
    "--min-expected",
    type=click.FloatRange(0),
    default=10.0,
    show_default=True,
    help="The pairs a join is expected to have, below which it is too thin to judge.",
)
@_FORCE
@_VERBOSE
@_QUIET
def junctions_command(
    pairs: Path,
    agp: Path,
    prefix: str,
    outdir: Path,
    window: int,
    min_ratio: float,
    min_expected: float,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Score each join of a genome's scaffolds by the Hi-C pairs across it."""
    from .junctions import score, write_outputs

    _refuse_existing([outdir / f"{prefix}.junctions.tsv"], force)
    parameters = {
        "pairs": pairs,
        "agp": agp,
        "prefix": prefix,
        "outdir": outdir,
        "window": window,
        "min_ratio": min_ratio,
        "min_expected": min_expected,
        "force": force,
    }
    with run_log(outdir, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            found = score(agp, pairs, window, min_ratio, min_expected)
        except (OSError, ValueError) as error:
            raise click.ClickException(str(error)) from error
        write_outputs(found, outdir, prefix)
        for junction in found:
            if junction.verdict == "unsupported":
                logger.warning(
                    "%s: the join of %s and %s at %d: %d pairs across, %.0f expected",
                    junction.scaffold,
                    junction.left.name,
                    junction.right.name,
                    junction.left.end,
                    junction.observed,
                    junction.expected,
                )


@cli.command(name="genes", context_settings=CONTEXT_SETTINGS)
@click.argument("genome", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("models", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_species,
    help="The species' ID, which starts every gene's name and names the outputs.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The folder of the annotation and the log.",
)
@click.option(
    "-g",
    "--genetic-code",
    type=click.Choice([str(c) for c in (1, 3, 4, 11, 12, 26)]),
    default="1",
    show_default=True,
    help="NCBI's translation table of the species' genes: 12 for the CUG-Ser clade, 26 for CUG-Ala.",
)
@click.option(
    "-s",
    "--source",
    default="plexembler",
    show_default=True,
    help="The GFF3 and GTF files' second column: what predicted the genes.",
)
@_FORCE
@_VERBOSE
@_QUIET
def genes_command(
    genome: Path,
    models: Path,
    prefix: str,
    outdir: Path,
    genetic_code: str,
    source: str,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Name a species' gene models; write its GFF3, GTF, proteins and CDS."""
    from .genes import annotate

    _refuse_existing([outdir / f"{prefix}.gff3.gz"], force)
    parameters = {
        "genome": genome,
        "models": models,
        "prefix": prefix,
        "outdir": outdir,
        "genetic_code": genetic_code,
        "source": source,
        "force": force,
    }
    with run_log(outdir, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            annotate(genome, models, prefix, outdir, int(genetic_code), source)
        except (OSError, ValueError, KeyError) as error:
            raise click.ClickException(str(error)) from error


def _labels(
    _context: click.Context, _parameter: click.Parameter, values: tuple[str, ...]
) -> dict[str, str]:
    """``LIBRARY=LABEL`` pairs, as a mapping."""
    labels = {}
    for value in values:
        library, separator, label = value.partition("=")
        if not separator or not library or not label:
            raise click.BadParameter(f"{value!r} is not LIBRARY=LABEL")
        labels[library] = label
    return labels


@cli.command(name="expression", context_settings=CONTEXT_SETTINGS)
@click.argument("counts", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_species,
    help="The species' ID, which starts its alignments' names and names the outputs.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The folder of the tables and the log.",
)
@click.option(
    "-l",
    "--label",
    "labels",
    multiple=True,
    callback=_labels,
    help="A library's column name, LIBRARY=LABEL (rna-ypd-1=YPD_1); the library's ID otherwise.",
)
@click.option(
    "--min-count",
    type=click.IntRange(1),
    default=10,
    show_default=True,
    help="The fewest pairs that make a gene detected, in MultiQC's table.",
)
@_FORCE
@_VERBOSE
@_QUIET
def expression_command(
    counts: Path,
    prefix: str,
    outdir: Path,
    labels: dict[str, str],
    min_count: int,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Write a species' counts and TPM, genes by libraries, from featureCounts' table."""
    from .expression import read_featurecounts, write_outputs

    _refuse_existing([outdir / f"{prefix}.counts.tsv"], force)
    parameters = {
        "counts": counts,
        "prefix": prefix,
        "outdir": outdir,
        "labels": ",".join(f"{k}={v}" for k, v in labels.items()),
        "min_count": min_count,
        "force": force,
    }
    with run_log(outdir, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            table = read_featurecounts(counts, prefix, labels)
        except (OSError, ValueError) as error:
            raise click.ClickException(str(error)) from error
        write_outputs(table, outdir, prefix, min_count)


@cli.command(name="report", context_settings=CONTEXT_SETTINGS)
@click.argument(
    "files",
    nargs=-1,
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "-p",
    "--prefix",
    required=True,
    callback=_species,
    help="The species' ID, which its files' names start with and which names the report.",
)
@click.option("--pool", help="The species' pool, which starts its pool's files' names.")
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The folder of the report and the log.",
)
@_FORCE
@_VERBOSE
@_QUIET
def report_command(
    files: tuple[Path, ...],
    prefix: str,
    pool: str | None,
    outdir: Path,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Write a species' report from the files the pipeline wrote for it."""
    from .report import write

    _refuse_existing([outdir / f"{prefix}.report.txt"], force)
    parameters = {
        "files": ",".join(str(f) for f in files),
        "prefix": prefix,
        "pool": pool,
        "outdir": outdir,
        "force": force,
    }
    with run_log(outdir, prefix, parameters, verbose=verbose, quiet=quiet):
        try:
            write(prefix, pool, list(files), outdir)
        except (OSError, ValueError, KeyError) as error:
            raise click.ClickException(str(error)) from error


def _pairs(
    _context: click.Context, parameter: click.Parameter, values: tuple[str, ...]
) -> list[tuple[str, str]]:
    """``KEY=VALUE`` pairs, in order."""
    pairs = []
    for value in values:
        key, separator, rest = value.partition("=")
        if not separator or not key or not rest:
            raise click.BadParameter(f"{value!r} is not {parameter.metavar}")
        pairs.append((key, rest))
    return pairs


@cli.command(name="browser", context_settings=CONTEXT_SETTINGS)
@click.option(
    "-s",
    "--species",
    "species",
    multiple=True,
    required=True,
    metavar="SPECIES=POOL",
    callback=_pairs,
    help="A species and its pool, SPECIES=POOL; once per species.",
)
@click.option(
    "-f",
    "--file",
    "files",
    multiple=True,
    metavar="SPECIES=PATH",
    callback=_pairs,
    help="A species' map, genome or track, by its path in the results; once per file.",
)
@click.option(
    "-o",
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="The top of the results, where the store goes, and the log's folder.",
)
@click.option(
    "-n",
    "--name",
    default="plexembler",
    show_default=True,
    help="The store's name, as chromosomix's picker shows it.",
)
@click.option(
    "--base",
    help="Where the results are served (https://…/results/), when they are on a bucket.",
)
@_FORCE
@_VERBOSE
@_QUIET
def browser_command(
    species: list[tuple[str, str]],
    files: list[tuple[str, str]],
    outdir: Path,
    name: str,
    base: str | None,
    force: bool,
    verbose: bool,
    quiet: bool,
) -> None:
    """Write chromosomix's store of every species' map, genome and tracks."""
    from .browser import STORE_NAME, gather, store, write

    _refuse_existing([outdir / STORE_NAME], force)
    if base is not None and not base.startswith(("https://", "http://")):
        raise click.BadParameter(f"{base!r} is not an http(s) URL", param_hint="--base")
    parameters = {
        "species": ",".join(f"{k}={v}" for k, v in species),
        "files": len(files),
        "outdir": outdir,
        "name": name,
        "base": base,
        "force": force,
    }
    with run_log(outdir, "browser", parameters, verbose=verbose, quiet=quiet):
        try:
            found = gather(dict(species), files)
        except ValueError as error:
            raise click.ClickException(str(error)) from error
        write(store(found, name, base), outdir)
