# File formats

What plexembler reads and writes: its samplesheets, the files of its own
commands, its report and its logs. The files the tools it runs write (FASTA,
GTF, `.mcool`, …) keep their tools' formats; [outputs](output.md) lists them.

## Samplesheets

The species (`--input`) and the pooled libraries (`--libraries`): CSV files
with a header, checked against `assets/schema_input.json` and
`assets/schema_libraries.json`. [Usage](usage.md) describes every column.

## The names table { #names-table }

[`plexembler rename`](cli.md#rename) writes, beside the renamed contigs, a
table of their names: tab-separated, with a header.

| Column | Description |
|---|---|
| `name` | The contig's name, starting with its species: `<species>_<contig>` |
| `original_name` | Its name before (the first word of its FASTA header) |
| `length_bp` | Its length, in bases |

## The candidates table { #candidates-table }

What [`plexembler select`](cli.md#select) reads: tab-separated, with a header,
a row per candidate. A relative path is read from the table's folder; an empty
cell, or `NA`, is a check the candidate does not have, which leaves its rule
out.

| Column | Required | Description |
|---|---|---|
| `candidate` | yes | The candidate's name: the tool that made it (`flye`, `hifiasm`) |
| `fasta` | yes | Its sequence, FASTA, gzipped or not |
| `busco` | no | BUSCO's short summary (`.txt` or `.json`), or compleasm's `summary.txt` |
| `qv` | no | Merqury's `.qv` file |
| `stats` | no | gfastats' summary, tabular (`--tabular`) or not |

## The selection table { #selection-table }

What [`plexembler select`](cli.md#select) writes beside the kept sequence,
`<species>.<step>.selection.tsv`: tab-separated, with a header, a row per
candidate in the order of the candidates table; `NA` for a number the
candidate lacks.

| Column | Description |
|---|---|
| `candidate` | The candidate's name |
| `kept` | `yes` for the candidate kept, `no` for the others |
| `busco_complete`, `busco_single`, `busco_duplicated`, `busco_fragmented`, `busco_missing` | BUSCO's classes, in percent of `busco_n` |
| `busco_n` | The BUSCO lineage's number of genes |
| `qv` | Merqury's QV (`inf` when it found no error) |
| `contigs`, `contig_n50`, `total_length` | gfastats' number of contigs, contig N50 and total length, in bases |
| `reason` | For the kept candidate, the rule that kept it: `the most single-copy complete BUSCOs`, `the highest QV`, `the longest contig N50`, `named by the samplesheet`, `the only candidate`, or a tie on every rule |

## The identity pairs { #identity-pairs }

What [`plexembler ani`](cli.md#ani) writes for a pool, `<pool>.ani_pairs.tsv`:
tab-separated, with a header, a row per pair of species, each pair once (the
first species first in alphabetical order). Its matrix, `<pool>.ani.tsv`, has
the same identities, a row and a column per species, `100.00` on its diagonal.

| Column | Description |
|---|---|
| `species_a`, `species_b` | The pair's species |
| `ani` | skani's average nucleotide identity, in percent; `NA` below skani's floor, about 80 % |
| `aligned_fraction_a`, `aligned_fraction_b` | The share of each genome that skani aligned, in percent |
| `verdict` | `too close` (above `--warn-above`, 95 %), `same genus` (from `--genus-above`, 85 %), or `distinct` |

## The genome's names { #genome-names }

What [`plexembler finalize`](cli.md#finalize) writes beside a species' genome,
`<species>.names.tsv`: tab-separated, with a header, a row per sequence of the
genome, in its order: the chromosomes, longest first, then the organelle, then
the unplaced scaffolds, longest first.

| Column | Description |
|---|---|
| `name` | Its name in the genome: `<species>_chr01`, …, `<species>_chrM`, `<species>_unplaced_01`, … |
| `original_name` | Its name in the scaffolder's output (`scaffold_1`, …), or the contig's, for a species finished from its contigs |
| `length_bp` | Its length, in bases |
| `role` | `chromosome`, `organelle` or `unplaced` |

## The report of a species

`plexembler report` (milestone 6 of the [plan](plan.md)) writes the report of
one species, its sample, from the outputs of its run: `<species>.report.txt`,
`<species>.report_mqc.yaml` and `<species>.report.html`, with chromokit's
report kit (CONVENTIONS.md § 12). The pipeline writes it for every species at
its end. Its sections and keys are documented here as they are added, and
change only with its `format_version`.

### `<sample>.report.txt` { #reporttxt }

Every number of the report, for scripts: `[section]`, then a key and its value
per line; `[table: name]`, then a header and a row per line; `NA` when
missing. Its `format_version` changes when a key, a column, a unit or a
MultiQC column changes or goes.

| Section | Key | Unit | MultiQC column |
|---|---|---|---|

### `<sample>.report_mqc.yaml` { #report_mqcyaml }

[MultiQC](https://multiqc.info)'s custom content (`id: plexembler`):
the sample's row of its General Statistics table.

### `<sample>.report.html` { #reporthtml }

A page to read the run with. Its plots are drawn by plotly.js, which it loads
from its CDN: offline, it keeps its numbers and tables.

## Logs

Each command that writes outputs writes its log beside them,
`<sample>.plexembler_<YYYYMMDDHHMMSS>.log`, where the sample is the run's
`-p/--prefix` (the species, for `rename`). It starts with hicstuff's header,

```text
## plexembler: v0.0.0 log file
## date: 2026-11-26 10:15:30
## cmd: plexembler rename contigs.fa -p Scer -o Scer.contigs.fa.gz
## fasta: contigs.fa
## prefix: Scer
## output: Scer.contigs.fa.gz
## names: Scer.contigs.names.tsv
## force: False
---
```

then holds one line per record, `<date> :: <level> :: <message>`, the date as
`%Y-%m-%d,%H:%M:%S`, and ends with how long the run took (`Done after 0h 0m 2s`).
The file gets INFO and above, and DEBUG with `-v`; `-q` quiets stderr only.
