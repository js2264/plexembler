# Working on plexembler

The instructions for coding agents, here and nowhere else (CONVENTIONS.md
§ 21): Copilot and Claude Code read this file. Do not add a `CLAUDE.md`, which
would switch Claude Code off reading it; a personal one holding the line
`@AGENTS.md` is ignored by git.

Reference genomes, gene annotations, Hi-C maps and expression tables for every
species of a pool, from per-species nanopore reads and pooled Hi-C, WGS and
RNA-seq. Read [README.md](README.md), [docs/concepts.md](docs/concepts.md)
and [the implementation plan](docs/plan.md) first; the project's conventions
are [koszullab/.github's CONVENTIONS.md](https://github.com/koszullab/.github/blob/main/CONVENTIONS.md),
adapted here to a Nextflow pipeline (deviations are written where they are).

## Commands

```bash
uv sync                                         # the environment
uv run pytest                                   # the tests, with coverage
uv run ruff check . && uv run ruff format .     # lint and format
uv run ty check src                             # types
uv run pre-commit install                       # the hooks
uv run --only-group docs zensical serve         # the documentation, at http://localhost:8000
uv run bump-my-version bump minor               # a release (docs/contributing.md)
```

The pipeline (Nextflow ≥ 26.04.0, nf-test, nf-core/tools, Docker):

```bash
nextflow lint .                                         # strict syntax, every script and config
nf-test test --profile +docker                          # pipeline, subworkflow and module tests
nextflow run . -profile test -stub --outdir results     # the plumbing, without tools
nf-core pipelines lint                                  # nf-core's checks (.nf-core.yml: exceptions)
nf-core pipelines schema docs nextflow_schema.json --output docs/parameters.md --force
```

pixi's script, with the checkout's plexembler for PyPI's (pixi ≥ 0.80.0;
docs/contributing.md):

```bash
sed "s|\"plexembler==[^\"]*\"|\"plexembler @ file://$PWD\"|" docs/run.py > /tmp/run.py
pixi run --script /tmp/run.py -- --pipeline . -profile test -preview --outdir results
```

## Layout

```text
main.nf, workflows/          the pipeline: its entry, and the workflow of the steps
subworkflows/local/          plexembler's subworkflows, one folder each, with tests/
modules/local/               plexembler's modules, nf-core's form, each with tests/
modules/nf-core/, subworkflows/nf-core/   installed by nf-core/tools, never edited
conf/                        base (resources), modules (args, publishing), test profiles
assets/                      the samplesheets' schemas, MultiQC's config
nextflow.config, nextflow_schema.json     parameters, profiles, manifest; their schema
src/plexembler/              the Python package: the commands the pipeline's steps run
src/plexembler/_report/      chromokit's report kit: never edited (chromokit sync replaces it)
tests/                       pytest (test_*.py) and the pipeline's nf-test (default.nf.test)
tests/data/                  the test data, and make_test_data.py, which writes them
benchmarks/                  the benchmark pool's generator and checks (benchmarks/README.md)
docs/                        the site's pages; docs/plan.md is the implementation plan,
                             docs/run.py pixi's script, which the site serves
```

## Rules

- **A pull request is reviewed as a whole, and tested at its head**: the
  maintainer squash-merges it. A review says how the whole change affects the
  pipeline and the package: behaviour, outputs, parameters, samplesheets,
  command line, docs. This applies to the brief given to an independent
  reviewer too.
- **A change a user can notice carries its changelog entry**, under
  `## Unreleased` in `docs/changelog.md`.
- **Every commit is a Conventional Commit:** `type(scope): summary`, and a body
  that says why. Scopes: `pipeline`, `cli`, `docs`, `ci`, or a step's name.
- **Coverage stays at 80 %**, overall and on the lines a pull request changes:
  `uv run diff-cover coverage.xml --compare-branch=origin/main --fail-under=80`.
- **The command line is tested in-process**, with click's `CliRunner`; a
  subprocess is for third-party programs only.
- **Diagnostics go through logging** (`logging.getLogger("plexembler")`),
  never `print`.
- **Nothing a user calls is renamed** (parameters, samplesheet columns,
  commands, options, output paths): an alias is added instead.
- **chromokit's files are never edited**: `overrides/`,
  `docs/stylesheets/chromokit.css`, `docs/javascripts/chromokit.js` and
  `src/plexembler/_report/`, which `chromokit.json` lists and `npx
  github:js2264/chromokit check` verifies. The template's files take their
  changes from `uvx copier update`.
- **nf-core's components are never edited**: `nf-core modules patch` records a
  change; a tool without an nf-core module gets a local one in nf-core's form.
- **Every species' contigs are named `<species>_<contig>`** from the assembly
  on: the pooled reads are split between the species by that prefix
  (docs/concepts.md). A species' ID has no underscore.
- **The version is the package's** (`pyproject.toml`), written again as the
  pipeline's `manifest.version`; bump-my-version changes both.
- **pixi's script runs the release** (`docs/run.py`): its plexembler is the
  version's, which bump-my-version writes, and its Nextflow the pipeline's
  (`tests/test_run.py` checks); a conda package added beside Nextflow brings
  back the conda-to-PyPI mapping (docs/contributing.md).
- **`plexembler run` passes everything else to Nextflow**: no pipeline
  parameter is named after one of its options (`pipeline`, `dry-run`), which
  would hide it (`tests/test_run.py` checks).
- **Every process has a `stub:` block**, and every local module and
  subworkflow an nf-test with a `-stub` test.
- **GPU runs cost money**: ask before starting one.

## Procedures

- **Adding a step:** its modules (nf-core's, or a local one), its subworkflow
  under `subworkflows/local/<name>/` with its nf-test, its options and
  publishing in `conf/modules.config` (a species' results under
  `species/<species>/`, a pool's under `pools/<pool>/`), its parameters in
  `nextflow.config` and `nextflow_schema.json` (then `docs/parameters.md`), its
  outputs in `docs/output.md`, its tools in `CITATIONS.md` and in
  `toolCitationText()`, its changelog entry, and its box ticked in
  `docs/plan.md`.
- **Adding a command:** its code, its tests (in-process), its section in
  `docs/cli.md` (`tests/test_docs.py` fails without it), and its changelog
  entry. A command that writes outputs takes `-o`, `--force`, `-p/--prefix`,
  `-v` and `-q`, and writes its log with `plexembler._log.run_log`.
- **Changing the test data:** edit `tests/data/make_test_data.py`, raise its
  `VERSION` when the files change, run it, and commit the files with it. The
  benchmark's data is never committed: raise `GENERATOR_VERSION` in
  `benchmarks/make_benchmark_data.py`, and `benchmarks.yml` publishes the new
  release `benchmark-data-v<N>`.
- **A milestone is done** when its boxes are ticked and its "Done when" passes
  on the benchmark: `python benchmarks/evaluate.py --milestone <N>`.
- **Reviewing a pull request:** as a whole, against these rules; run the tests
  at its head.
