"""The pipeline, run with Nextflow: ``plexembler run``.

``plexembler run ARGS…`` runs ``nextflow run js2264/plexembler -r v<version>
ARGS…``: the pipeline of the package's own version, every argument but
plexembler's own options going to Nextflow as it is (docs/cli.md). With
Nextflow beside it, as pixi's script installs both (docs/run.py), the package
is all a machine needs to run the pipeline, but a container engine.

plexembler then hands its process over to Nextflow (``exec``): Nextflow's exit
status, signals and logs are the run's, as when Nextflow is launched by hand.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import sys
import sysconfig
from collections.abc import Collection, Sequence
from importlib.metadata import distribution
from pathlib import Path
from typing import Any, NoReturn
from urllib.parse import urlparse
from urllib.request import url2pathname

from . import __version__

#: The pipeline's repository, as Nextflow names it.
REPOSITORY = "js2264/plexembler"
#: Nextflow's options that name the revision to run: given one, plexembler adds none.
REVISION_OPTIONS = frozenset({"-r", "-revision"})
#: A GitHub repository's URL, as pip records an install from git:
#: https://github.com/<owner>/<name>.git, ssh://git@github.com/<owner>/<name>.
_GITHUB = re.compile(r"github\.com[:/](?P<owner>[^/]+)/(?P<name>[^/]+?)(?:\.git)?/?$")


def split_arguments(
    arguments: Sequence[str], flags: Collection[str], valued: Collection[str]
) -> tuple[list[str], list[str]]:
    """plexembler's own options, and the arguments for Nextflow, as they are.

    An argument is plexembler's when it is one of its options, word for word
    (``--pipeline=…`` too), and the word after an option that takes a value is
    that value; every other argument is Nextflow's, in its order. After ``--``,
    every argument is Nextflow's.

    Parameters
    ----------
    arguments
        The command line after ``plexembler run``.
    flags
        plexembler's options that take no value (``--dry-run``, ``-h``).
    valued
        plexembler's options that take one (``--pipeline``).

    Returns
    -------
    tuple of list of str
        plexembler's arguments, and Nextflow's.

    Raises
    ------
    ValueError
        When an option that takes a value ends the command line.
    """
    own: list[str] = []
    nextflow: list[str] = []
    words = iter(arguments)
    for word in words:
        if word == "--":
            nextflow.extend(words)
            break
        if word in flags:
            own.append(word)
        elif word in valued:
            value = next(words, "--")
            if value == "--":
                raise ValueError(f"Option '{word}' requires an argument.")
            own += [word, value]
        elif word.split("=", 1)[0] in valued:
            own.append(word)
        else:
            nextflow.append(word)
    return own, nextflow


def installed_from() -> dict[str, Any] | None:
    """Where the package was installed from, as pip, uv and pixi record it (PEP 610).

    Returns
    -------
    dict or None
        The installation's ``direct_url.json``: a folder (``dir_info``), a git
        repository (``vcs_info``) or an archive; None for a package installed
        from an index (PyPI, bioconda) or with an unreadable record.
    """
    text = distribution("plexembler").read_text("direct_url.json")
    if not text:
        return None
    try:
        origin = json.loads(text)
    except json.JSONDecodeError:
        return None
    return origin if isinstance(origin, dict) else None


def default_pipeline() -> tuple[str, str | None]:
    """The pipeline of this installation of plexembler, and its revision.

    - Installed from a folder that holds the pipeline (``uv sync``, ``pip
      install -e .``, pixi from a git repository): that folder, as it is.
    - Installed from a git repository by pip or uv: that repository, at the
      commit installed.
    - Otherwise, a release (PyPI, bioconda, or a folder that is gone, such as a
      conda package's build): ``js2264/plexembler`` at its tag, ``v<version>``.

    Returns
    -------
    tuple of str and str or None
        The pipeline, a folder or a repository, and its revision: None for a
        folder.
    """
    origin = installed_from() or {}
    url = str(origin.get("url", ""))
    if "dir_info" in origin and url.startswith("file:"):
        folder = Path(url2pathname(urlparse(url).path))
        if (folder / "main.nf").is_file():
            return str(folder), None
    vcs = origin.get("vcs_info") or {}
    if vcs.get("vcs") == "git" and vcs.get("commit_id"):
        github = _GITHUB.search(url)
        repository = f"{github['owner']}/{github['name']}" if github else url
        return repository, str(vcs["commit_id"])
    return REPOSITORY, f"v{__version__}"


def pipeline_to_run(
    pipeline: str | None, nextflow_args: Sequence[str]
) -> tuple[str, str | None]:
    """The pipeline to run, and the revision plexembler gives it with ``-r``.

    Parameters
    ----------
    pipeline
        ``--pipeline``: a folder, run as it is, or a repository, run at the
        package's version; None for the installation's own
        (:func:`default_pipeline`).
    nextflow_args
        The arguments for Nextflow: a revision among them (``-r``) is the one
        run, of the repository.

    Returns
    -------
    tuple of str and str or None
        The pipeline, and the revision to add: None for a folder, or when the
        arguments name one.
    """
    revision_given = any(argument in REVISION_OPTIONS for argument in nextflow_args)
    if pipeline is None:
        pipeline, revision = default_pipeline()
        if revision is None and revision_given:
            # A revision is a repository's: a folder has none.
            pipeline = REPOSITORY
    elif Path(pipeline).is_dir():
        revision = None
    else:
        revision = f"v{__version__}"
    return pipeline, None if revision_given else revision


def find_nextflow() -> str | None:
    """Nextflow's executable: beside plexembler's Python, else on PATH.

    pixi and conda install Nextflow in plexembler's own environment, whose
    folder of executables need not be on PATH (an environment that is not
    activated).

    Returns
    -------
    str or None
        Nextflow's path; None when it is nowhere.
    """
    beside = shutil.which("nextflow", path=sysconfig.get_path("scripts"))
    return beside or shutil.which("nextflow")


def nextflow_command(
    nextflow: str, pipeline: str, revision: str | None, nextflow_args: Sequence[str]
) -> list[str]:
    """The command run: ``nextflow run <pipeline> [-r <revision>] <arguments>``."""
    revision_args = ["-r", revision] if revision is not None else []
    return [nextflow, "run", pipeline, *revision_args, *nextflow_args]


def environment(nextflow: str) -> dict[str, str]:
    """The environment Nextflow runs in: plexembler's own.

    When Nextflow is the one of plexembler's environment, that environment's
    executables come first on PATH, so that Nextflow finds its Java there too.
    """
    env = dict(os.environ)
    folder = Path(nextflow).parent
    if folder == Path(sysconfig.get_path("scripts")):
        env["PATH"] = os.pathsep.join(filter(None, [str(folder), env.get("PATH")]))
    return env


def launch(command: Sequence[str]) -> NoReturn:
    """Hand the process over to Nextflow (``exec``): it never comes back.

    Raises
    ------
    OSError
        When the command cannot be run.
    """
    sys.stdout.flush()
    sys.stderr.flush()
    os.execve(command[0], list(command), environment(command[0]))
