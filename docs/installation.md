# Installation

plexembler is a [Nextflow](https://www.nextflow.io) pipeline: Nextflow fetches
it, and every tool it runs comes in a container or a conda environment. It
runs on Linux, and on macOS for small tests.

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
hold. You need it yourself only to work on plexembler
([contributing](contributing.md)), or to run a command by hand:

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
