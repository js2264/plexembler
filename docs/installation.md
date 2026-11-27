# Installation

plexembler is a [Nextflow](https://www.nextflow.io) pipeline: Nextflow fetches
it, and every tool it runs comes in a container or a conda environment. It
runs on Linux, and on macOS for small tests.

## With pixi

[pixi](https://pixi.sh) runs the pipeline from a script of a few lines,
[`run.py`](run.py), that names what it needs: plexembler, from PyPI, and
Nextflow, from bioconda, which brings Java. pixi makes their environment the
first time and keeps it: with a container engine for the pipeline's tools
(Docker, Apptainer, Singularity, Podman), pixi is all a machine needs.

```bash
curl -fsSL https://pixi.sh/install.sh | sh         # pixi itself, once
pixi run --script https://js2264.github.io/plexembler/run.py -- \
    -profile test,docker --outdir results
```

Every argument after `--` is `plexembler run`'s ([command line](cli.md#run)),
which runs the pipeline of the latest release: `nextflow run js2264/plexembler
-r v<version>`, with its arguments. The `--` is needed: pixi would read
`-profile` as an option of its own.

- **A release:** the script of its tag runs that release.

    ```bash
    pixi run --script https://raw.githubusercontent.com/js2264/plexembler/v0.1.0/docs/run.py -- \
        -profile docker --input species.csv --libraries libraries.csv --outdir results
    ```

- **pixi 0.80.0 or newer** runs it: a PEP 723 script, with conda packages.
- **pixi fetches the script at each run**, and keeps its environment. Where
  the network is not to be relied on, as on a cluster's nodes, download it
  once and run the file: `curl -O https://js2264.github.io/plexembler/run.py`,
  then `pixi run --script run.py -- …`.

!!! note "To come"

    The script takes plexembler from PyPI, where it is published from its
    first release on ([the plan](plan.md#t8-plexembler-run-and-pixi)).

## The pipeline

You need:

- **Nextflow 26.04.0 or newer**, and Java 17 or newer;
- **a container engine** (Docker, Apptainer, Singularity, Podman), or conda.

=== "Nextflow's installer"

    ```bash
    curl -s https://get.nextflow.io | bash    # writes ./nextflow
    mv nextflow ~/bin/                        # anywhere on PATH
    ```

=== "conda"

    ```bash
    conda install -c conda-forge -c bioconda nextflow
    ```

Then fetch the pipeline, and run its test:

```bash
nextflow pull js2264/plexembler
nextflow run js2264/plexembler -profile test,docker --outdir results
```

`nextflow run js2264/plexembler -r <version>` runs a release; without `-r`,
the latest commit of the default branch.

## The Python package

The commands the pipeline's steps run (`plexembler <command>`, the
[command line](cli.md)) are a Python package, which the pipeline's containers
hold, and which runs the pipeline itself (`plexembler run`). You need it
yourself only to work on plexembler ([contributing](contributing.md)), to run
a command by hand, or to run the pipeline with `plexembler run` and a Nextflow
of your own on `PATH` (pixi installs both, above):

=== "uv"

    ```bash
    uv tool install git+https://github.com/js2264/plexembler
    ```

=== "pip"

    ```bash
    pip install git+https://github.com/js2264/plexembler
    ```

`plexembler --version` says which version is installed. It needs Python 3.10
or newer.
