# /// script
# requires-python = ">=3.10"
# # The release, from PyPI: bump-my-version writes it (pyproject.toml).
# dependencies = ["plexembler==0.0.0"]
#
# [tool.pixi.workspace]
# channels = ["conda-forge", "bioconda"]
# # No conda package below is a PyPI one too: no mapping to fetch (prefix.dev).
# conda-pypi-map = false
#
# [tool.pixi.dependencies]
# # The pipeline's runner, manifest.nextflowVersion; it brings Java.
# nextflow = ">=26.04.0"
# ///
"""plexembler run, with nothing installed but pixi (docs/installation.md).

pixi makes an environment of plexembler, from PyPI, Nextflow, from bioconda,
and Java, from conda-forge, the first time, and keeps it:

    pixi run --script https://js2264.github.io/plexembler/run.py -- \\
        -profile docker --input species.csv --libraries libraries.csv --outdir results

Every argument after ``--`` is ``plexembler run``'s (docs/cli.md): the
pipeline of the release this file names. Its copy at a release's tag names
that release: https://raw.githubusercontent.com/js2264/plexembler/vX.Y.Z/docs/run.py
"""

import sys

from plexembler.cli import cli

cli(["run", *sys.argv[1:]], prog_name="plexembler")
