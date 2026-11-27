# The image of plexembler (CONVENTIONS.md § 15), which release.yml
# builds from the released wheel (dist/*.whl): it holds what PyPI holds, and
# the programs pip cannot install. Dependabot updates the base.
ARG BASE=ghcr.io/mamba-org/micromamba:2.10.0
FROM ${BASE}

# samtools: the pipeline's steps stream alignments into plexembler's commands
# (`samtools view | plexembler split -`, docs/plan.md, T1). procps: ps, which
# Nextflow runs in every step's container to measure it.
RUN micromamba install -y -n base -c conda-forge -c bioconda python=3.14 pip procps-ng samtools=1.22.1 \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER dist/*.whl /tmp/dist/
RUN micromamba run -n base pip install --no-cache-dir /tmp/dist/*.whl \
    && rm -rf /tmp/dist

# The environment on PATH without the entrypoint too: Nextflow runs a step's
# script with bash, and Apptainer and Singularity skip the entrypoint.
ENV PATH=/opt/conda/bin:$PATH

WORKDIR /home/mambauser
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
