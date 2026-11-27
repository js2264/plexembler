# Contributing

Bug reports, fixes, features and corrections to these pages are all welcome,
on [GitHub](https://github.com/js2264/plexembler). plexembler follows
the project's [conventions](https://github.com/koszullab/.github/blob/main/CONVENTIONS.md).

The repository holds two things that share a version: the **pipeline**
(Nextflow, on nf-core's template: `main.nf`, `workflows/`, `subworkflows/`,
`modules/`, `conf/`, `assets/`) and the **Python package** (`src/plexembler/`),
whose commands the pipeline's own steps run. [The plan](plan.md) says what
comes next.

## Setting up

The Python environment is [uv](https://docs.astral.sh/uv/)'s; `uv.lock` pins
every dependency and is committed. The pipeline needs
[Nextflow](https://www.nextflow.io) (26.04.0 or newer, with Java 17 or newer),
[nf-test](https://www.nf-test.com), [nf-core/tools](https://nf-co.re/docs/nf-core-tools)
and a container engine (Docker locally, as in CI).

```bash
git clone https://github.com/js2264/plexembler.git
cd plexembler
uv sync                       # the environment: plexembler, editable, and the dev group
uv run pre-commit install     # the hooks: ruff and `nextflow lint`, on every commit
uv tool install nf-core       # nf-core/tools: modules, lint, schema
```

## Testing

### The Python package

```bash
uv run pytest                                   # the tests, with coverage
uv run ruff check . && uv run ruff format .     # lint and format
uv run ty check src                             # types
```

- The command line is tested in-process, with click's `CliRunner`; a
  subprocess is for third-party programs only, and a test that needs a
  program that is not installed is skipped.
- A known bug is an `xfail` test naming its issue in `reason=`.
- Coverage stays at 80 %, overall and on the lines a pull request changes.
  Before pushing:

    ```bash
    uv run pytest
    uv run diff-cover coverage.xml --compare-branch=origin/main --fail-under=80
    ```

### The pipeline

```bash
nextflow lint .                                       # the scripts and configs (strict syntax)
nf-test test --profile +docker                        # every test: pipeline, subworkflows, modules
nf-test test --profile +docker --tag ont_reads        # one subworkflow's
nextflow run . -profile test,docker --outdir results  # the test profile, by hand
nextflow run . -profile test -stub --outdir results   # the plumbing alone, without tools
nf-core pipelines lint                                # nf-core's checks (.nf-core.yml lists the exceptions)
```

- **Every local module and subworkflow has its nf-test**, under `tests/`
  beside its `main.nf`, with a `-stub` test that needs no tool. The pipeline's
  own tests are `tests/default.nf.test`.
- **Tests assert what a step must produce**; snapshots (`.nf.test.snap`) come
  once a step's outputs are stable, recorded with `nf-test test --update-snapshot`.
- **The test data** are written by `tests/data/make_test_data.py`, beside it:
  two species of a few dozen kilobases and their pool. Run it again after
  changing it, and raise its `VERSION` when what it writes changes. Larger
  data go to a release asset (`test-data-v<N>`), never into git
  (CONVENTIONS.md § 13).

### The conda package

`pixi.toml` is plexembler's conda package, which pixi builds from the
repository ([installation](installation.md#with-pixi)): the Python package,
with Nextflow from bioconda. `pipeline.yml`'s `pixi` job builds it and runs the
test profile through it, as a machine with pixi alone would; with
[pixi](https://pixi.sh) 0.85.0 or newer:

```bash
pixi build --output-dir conda-dist          # the package, from the checkout
pixi exec -c conda-forge -c bioconda -s conda-dist/plexembler-*.conda -- \
    plexembler run --pipeline . -profile test -preview --outdir results
```

- **Its run-dependencies** are `pyproject.toml`'s dependencies, as conda-forge
  names them, and the Nextflow of `nextflow.config`'s
  `manifest.nextflowVersion`: `tests/test_run.py` fails when they differ.
- **pixi keeps a package by its file name**, which every build of a version
  shares, `--force-reinstall` or not: run a new build with a cache of its own,
  `PIXI_CACHE_DIR=$(mktemp -d) pixi exec …`.

### nf-core's modules and subworkflows

```bash
nf-core modules install <tool>          # a module of nf-core/modules, into modules/nf-core/
nf-core modules update --all            # newer versions
nf-core modules patch <tool>            # record a change to one, as a patch
nf-core subworkflows install <name>
```

They are never edited by hand: `modules.json` records each one's version, and
`nf-core pipelines lint` checks them. A tool without an nf-core module gets a
local one, `modules/local/<tool>/main.nf`, in nf-core's form (a `stub:` block,
its versions in the `versions` topic), and goes upstream to nf-core/modules
when it is useful beyond plexembler.

## Documenting

The site is built from `docs/` by [zensical](https://zensical.org), in
[chromokit](https://github.com/js2264/chromokit)'s look, and published at
<https://js2264.github.io/plexembler/>.

```bash
uv run --only-group docs zensical serve     # http://localhost:8000
```

- **Pages:** the Markdown files of `docs/`, listed in order by `nav` in
  `zensical.toml`. A new page goes in both.
- **The command line:** `docs/cli.md`, written by hand. `tests/test_docs.py`
  fails when it misses a command or an option, or documents one that does not
  exist.
- **The Python API:** `::: plexembler.<module>` directives in `docs/api/`,
  which mkdocstrings expands from the docstrings (numpy style).
- **The look:** chromokit's files (`overrides/`, `docs/stylesheets/chromokit.css`,
  `docs/javascripts/chromokit.js`, listed in `chromokit.json`) are never
  edited; the site's own styles go in `docs/stylesheets/site.css`.

    ```bash
    npx github:js2264/chromokit sync     # a newer chromokit
    npx github:js2264/chromokit check    # nothing of chromokit's edited
    ```

## Changing things

- **Every commit is a [Conventional Commit](https://www.conventionalcommits.org/en/v1.0.0/):**
  `type(scope): summary`, and a body that says why.
- **A change a user can notice** gets its entry under `## Unreleased` in
  `docs/changelog.md`, ending with a link to its pull request.
- **A pull request is reviewed as a whole, tested at its head, and squash-merged.**
- **Nothing a user calls is renamed:** commands, options, output files,
  modules. An alias is added instead.

## Releasing

From an up-to-date `main` with green CI:

```bash
uv run bump-my-version bump minor     # or patch; never major: plexembler stays at 0.x (plan)
git push --follow-tags
```

bump-my-version sets the version (`pyproject.toml` and `uv.lock`, the
pipeline's `manifest.version` in `nextflow.config`, `zensical.toml`'s
`release`, `CITATION.cff`), turns `## Unreleased` into the version's section,
commits `chore(release): vX.Y.Z` and tags `vX.Y.Z`. The tag starts
`release.yml`: CI and the pipeline's tests on the tag, the sdist and the wheel
checked, the GitHub release from the changelog, and the image,
`ghcr.io/js2264/plexembler:X.Y.Z` and `latest`, which the pipeline's own steps
run in. PyPI's trusted publishing waits for the repository variable
`PUBLISH_PYPI` to be `true` (docs/plan.md, T1). The pipeline is released by
the same tag: `nextflow run js2264/plexembler -r vX.Y.Z`, which `plexembler
run` of version X.Y.Z runs. A minor release that breaks something says so
under *Breaking changes* in its changelog section, and comes with
`docs/upgrading.md` ([the plan's versions](plan.md#versions-0x-for-good)).

### bioconda

The project's conda packages are bioconda's (CONVENTIONS.md § 0). plexembler's
recipe is to be submitted to [bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
once, after the first release ([the plan](plan.md#t8-plexembler-run-and-pixi)):
`recipes/plexembler/meta.yaml`, its source the release's archive. bioconda's
bot then follows each release; when the dependencies change (`pyproject.toml`,
`pixi.toml`), the recipe's are changed by hand, the same day. Its `sha256` is
the archive's (`curl -sL <url> | sha256sum`); with a local `path:` as its
source, `conda-build` builds and tests it.

```yaml title="recipes/plexembler/meta.yaml"
{% set name = "plexembler" %}
{% set version = "0.1.0" %}

package:
  name: {{ name }}
  version: {{ version }}

source:
  url: https://github.com/js2264/plexembler/archive/v{{ version }}.tar.gz
  sha256: <the archive's>

build:
  number: 0
  noarch: python
  script: {{ PYTHON }} -m pip install . --no-deps --no-build-isolation --no-cache-dir -vvv
  run_exports:
    # 0.x for good: a minor release may break.
    - {{ pin_subpackage(name, max_pin="x.x") }}

requirements:
  host:
    - python >=3.10
    - pip
    - hatchling >=1.27
  run:
    - python >=3.10
    - click >=8.1
    - rich-click >=1.9
    - jinja2 >=3.1
    - pyyaml >=6.0
    - nextflow >=26.04.0

test:
  imports:
    - plexembler
  commands:
    - plexembler --help
    - plexembler run --dry-run -profile test --outdir results

about:
  home: https://github.com/js2264/plexembler
  license: MIT
  license_family: MIT
  license_file: LICENSE
  summary: Genomes, gene annotations, Hi-C maps and expression tables for every species of a pool.
  description: |
    Reference genomes, gene annotations, Hi-C maps and expression tables for
    every species of a pool, from per-species nanopore reads and pooled Hi-C,
    WGS and RNA-seq: a Nextflow pipeline, run by `plexembler run`.
  doc_url: https://js2264.github.io/plexembler/
  dev_url: https://github.com/js2264/plexembler

extra:
  recipe-maintainers:
    - js2264
```

The project's template (chromokit's) brings its changes with `uvx copier update`.
