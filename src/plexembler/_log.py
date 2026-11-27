"""The log of a command, beside its outputs (CONVENTIONS.md § 11).

Every command that writes outputs writes ``<sample>.plexembler_<YYYYMMDDHHMMSS>.log``
in the directory of its outputs, headed as hicstuff's logs are: the version,
the date, the command line, one line per parameter, then ``---``. Then one line
per record, and a last line that says what the run wrote and how long it took.
"""

from __future__ import annotations

import logging
import sys
import time
from collections.abc import Generator, Mapping
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path

from . import __version__

#: The logger of every module of the package; only the command line configures it.
LOGGER_NAME = "plexembler"
#: The format of a record, in the file and on stderr.
RECORD_FORMAT = "%(asctime)s :: %(levelname)s :: %(message)s"
#: The format of a record's date.
DATE_FORMAT = "%Y-%m-%d,%H:%M:%S"

logger = logging.getLogger(LOGGER_NAME)


def log_file_name(outdir: Path, sample: str | None, started: datetime) -> Path:
    """The log's path: ``<sample>.plexembler_<YYYYMMDDHHMMSS>.log`` in ``outdir``.

    Parameters
    ----------
    outdir
        The directory of the command's outputs.
    sample
        The run's name (its ``-p/--prefix``), or None for a command without one.
    started
        When the command started, in local time.

    Returns
    -------
    Path
        The log's path; ``plexembler_<YYYYMMDDHHMMSS>.log`` without a sample.
    """
    stamp = started.strftime("%Y%m%d%H%M%S")
    name = f"{sample}.plexembler_{stamp}.log" if sample else f"plexembler_{stamp}.log"
    return Path(outdir) / name


def log_header(
    started: datetime, command: str, parameters: Mapping[str, object]
) -> str:
    """The log's header, as hicstuff writes it.

    Parameters
    ----------
    started
        When the command started, in local time.
    command
        The command line, as it was typed.
    parameters
        Every parameter of the command, as given or by default.

    Returns
    -------
    str
        The header's lines, the last one ``---``.
    """
    lines = [
        f"## plexembler: v{__version__} log file",
        f"## date: {started:%Y-%m-%d %H:%M:%S}",
        f"## cmd: {command}",
    ]
    lines += [f"## {key}: {value}" for key, value in parameters.items()]
    lines.append("---")
    return "\n".join(lines) + "\n"


def elapsed(seconds: float) -> str:
    """A duration as the log's last line gives it: ``0h 17m 31s``."""
    minutes, secs = divmod(round(seconds), 60)
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes}m {secs}s"


@contextmanager
def stderr_log() -> Generator[None]:
    """Log a command's run to stderr alone, INFO and above, for as long as it runs.

    For a command that writes no log file (CONVENTIONS.md § 11), such as
    ``plexembler run``, whose run Nextflow logs.
    """
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter(RECORD_FORMAT, datefmt=DATE_FORMAT))
    handler.setLevel(logging.INFO)
    previous_level = logger.level
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    try:
        yield
    finally:
        logger.removeHandler(handler)
        handler.close()
        logger.setLevel(previous_level)


@contextmanager
def run_log(
    outdir: Path,
    sample: str | None,
    parameters: Mapping[str, object],
    *,
    verbose: bool = False,
    quiet: bool = False,
) -> Generator[Path]:
    """Log a command's run to stderr and to its log file, for as long as it runs.

    The file always gets INFO and above; ``verbose`` adds DEBUG, to the file and
    to stderr; ``quiet`` limits stderr to warnings, never the file.

    Parameters
    ----------
    outdir
        The directory of the command's outputs, where the log goes.
    sample
        The run's name, which names the log.
    parameters
        Every parameter of the command, for the header.
    verbose, quiet
        The command's ``-v/--verbose`` and ``-q/--quiet``.

    Yields
    ------
    Path
        The log's path.
    """
    started = datetime.now()
    clock = time.monotonic()
    path = log_file_name(outdir, sample, started)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        log_header(started, " ".join(sys.argv), parameters), encoding="utf-8"
    )

    formatter = logging.Formatter(RECORD_FORMAT, datefmt=DATE_FORMAT)
    to_file = logging.FileHandler(path, mode="a", encoding="utf-8")
    to_file.setLevel(logging.DEBUG if verbose else logging.INFO)
    to_stderr = logging.StreamHandler(sys.stderr)
    to_stderr.setLevel(
        logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO
    )
    previous_level = logger.level
    for handler in (to_file, to_stderr):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)
    try:
        yield path
        logger.info("Done after %s", elapsed(time.monotonic() - clock))
    finally:
        for handler in (to_file, to_stderr):
            logger.removeHandler(handler)
            handler.close()
        logger.setLevel(previous_level)
