# Quickstart

## The test

Two species of a few dozen kilobases, and a pool of three libraries made from
both (`tests/data/`, written by `tests/data/make_test_data.py`):

=== "Nextflow"

    ```bash
    nextflow run js2264/plexembler -profile test,docker --outdir results
    ```

=== "pixi"

    ```bash
    pixi run --script https://js2264.github.io/plexembler/run.py -- \
        -profile test,docker --outdir results
    ```

    The same run, of the latest release, with nothing installed but pixi
    ([installation](installation.md#with-pixi)).

Open `results/multiqc/multiqc_report.html`.

## Your pool

1. Write the **species** samplesheet: one row per nanopore run of a species.

    ```csv title="species.csv"
    species,pool,ont,genome_size
    Scer,poolA,ont/Scer.fastq.gz,12m
    Kmar,poolA,ont/Kmar.fastq.gz,11m
    ```

2. Write the **pooled libraries** samplesheet: one row per pair of FASTQ files.

    ```csv title="libraries.csv"
    library,pool,type,fastq_1,fastq_2,enzymes
    hic,poolA,hic,hic_R1.fastq.gz,hic_R2.fastq.gz,"DpnII,HinfI"
    wgs,poolA,wgs,wgs_R1.fastq.gz,wgs_R2.fastq.gz,
    ```

    A relative path is read from the samplesheet's folder.

3. Run it, with the profile of your machine ([on a workstation or a
   cluster](guides/compute.md)):

    === "Nextflow"

        ```bash
        nextflow run js2264/plexembler -profile docker \
            --input species.csv --libraries libraries.csv --outdir results
        ```

    === "pixi"

        ```bash
        pixi run --script https://js2264.github.io/plexembler/run.py -- -profile docker \
            --input species.csv --libraries libraries.csv --outdir results
        ```

[Usage](usage.md) describes every column, and [outputs](output.md) every
result.
