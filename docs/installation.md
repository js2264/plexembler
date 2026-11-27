# Installation

plexembler is a [Nextflow](https://www.nextflow.io) pipeline: Nextflow fetches
it, and every tool it runs comes in a container or a conda environment. It
runs on Linux, and on macOS for small tests.

## With pixi

[pixi](https://pixi.sh) installs plexembler with Nextflow and Java, from
conda-forge and bioconda, in an environment of their own: with a container
engine for the pipeline's tools (Docker, Apptainer, Singularity, Podman), it
is all a machine needs.

```bash
curl -fsSL https://pixi.sh/install.sh | sh         # pixi itself, once
pixi global install -c conda-forge -c bioconda \
    --git https://github.com/js2264/plexembler.git --tag v0.1.0
plexembler run -profile test,docker --outdir results
```

`plexembler run` is `nextflow run js2264/plexembler -r v0.1.0`, the pipeline
of the release installed: every argument after it is Nextflow's or the
pipeline's ([command line](cli.md#run)). Another release is installed the same
way, with its tag; `--branch main` installs the latest commit, whose pipeline
`plexembler run` runs. pixi builds the package from the repository's
`pixi.toml`, with its `pixi-build` preview: pixi 0.85.0 or newer.

!!! note "To come: nothing to install"

    Once plexembler is on bioconda ([the plan](plan.md#t8-plexembler-run-and-pixi)),
    `pixi exec` runs it from a temporary environment that pixi keeps, with
    nothing installed:

    ```bash
    pixi exec -c conda-forge -c bioconda plexembler run -profile docker \
        --input species.csv --libraries libraries.csv --outdir results
    ```

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
