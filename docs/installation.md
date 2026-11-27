# Installation

plexembler is a [Nextflow](https://www.nextflow.io) pipeline: Nextflow fetches
it, and every tool it runs comes in a container or a conda environment. It
runs on Linux, and on macOS for small tests.

## With uv

[uv](https://docs.astral.sh/uv/) runs plexembler from PyPI, with its
dependencies, Nextflow's launcher among them (PyPI's `nextflow`), which
fetches Nextflow at its first run. With Java, and a container engine for the
pipeline's tools (Docker, Apptainer, Singularity, Podman) or micromamba for
their conda environments ([without a container
engine](usage.md#without-a-container-engine)), uv is all a machine needs:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh    # uv itself, once
uvx plexembler run -profile test,docker --outdir results
```

`plexembler run` runs the pipeline of the release that uv installed, `nextflow
run js2264/plexembler -r v<version>`, with every other argument as it is
([command line](cli.md#run)).

- **A release:** `uvx plexembler@0.1.0 run …` runs that release.
- **Installed for good:** `uv tool install plexembler`, then `plexembler run
  …`; `uv tool upgrade plexembler` takes the latest release.
- **Java 17 or newer** is the machine's: `java` on `PATH`, or `JAVA_HOME`
  (`NXF_JAVA_HOME` to choose another). On a cluster, often `module load java`;
  elsewhere, the system's OpenJDK package, Homebrew's `openjdk`, or
  [SDKMAN!](https://sdkman.io).
- **Nextflow is fetched once**, into `~/.nextflow`, from nextflow.io: the
  version of the launcher, PyPI's `nextflow`, which uv took with plexembler.

!!! note "To come"

    uvx takes plexembler from PyPI, where it is published from its first
    release on ([the plan](plan.md#t8-plexembler-run-and-uv)).

## The pipeline

You need:

- **Nextflow 26.04.0 or newer**, and Java 17 or newer;
- **a container engine** (Docker, Apptainer, Singularity, Podman), or conda,
  mamba or micromamba ([without a container
  engine](usage.md#without-a-container-engine)).

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
a command by hand, or to run the pipeline with `plexembler run`, which brings
Nextflow's launcher with it (above):

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
