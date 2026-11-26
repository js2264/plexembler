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
