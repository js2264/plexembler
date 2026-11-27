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
| `length_bp` | Its length, in bases, as written: one copy's, for a contig cut to one |
| `original_length_bp` | Its length as assembled |
| `copies` | The copies of its unit it was assembled as, cut to one; `1` for a contig left whole |
| `circular` | `yes` when the assembler's table says it is circular, or when it was cut (a circle's copies, the last closing on the first); `no` when the table says it is linear; `NA` when nothing says |

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
| `kmers` | no | A polishing's changed k-mers: a header, `fixed` and `introduced`, and a row of their two counts (the pipeline's `MERYL_CHANGES`) |

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
| `kmers_fixed`, `kmers_introduced` | A polishing's changed k-mers that the held-out reads lack: the draft's that it removed, and its own |
| `reason` | For the kept candidate, the rule that kept it: `the most single-copy complete BUSCOs`, `the highest QV`, `the longest contig N50`, `named by the samplesheet`, `the only candidate`, or a tie on every rule; for a polishing, its changed k-mers (`the held-out reads lack …`), or that none was better than the draft |

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
| `original_name` | Its name in the scaffolder's output (`scaffold_1`, …; `scaffold_1.1`, `scaffold_1.2`, … for the sides of a scaffold whose unsupported joins were cut), or the contig's, for a species finished from its contigs |
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
| `report` | `format_version`, `sample`, `pool`, `files`, `plexembler_report_version`, `written` |  |  |
| `assembly` | `candidates`, `kept`, `reason`, `contigs`, `contig_n50` | bp (N50) |  |
| `polishing` | `kept`, `qv`, `reason` |  |  |
| `screen` | `keep`, `move`, `drop`, `organelle` | contigs |  |
| `genome` | `length_bp` | bp |  |
| `genome` | `length_mb` | Mb | Genome (Mb) |
| `genome` | `chromosomes` |  | Chromosomes |
| `genome` | `chromosomes_bp`, `organelles`, `unplaced` | bp; sequences |  |
| `quality` | `busco_complete` | % | BUSCO (%) |
| `quality` | `busco_single`, `busco_duplicated`, `busco_fragmented`, `busco_missing`, `busco_n` | %; genes |  |
| `quality` | `qv` | Phred | QV |
| `scaffolding` | `joins`, `supported`, `unsupported`, `thin` | joins |  |
| `annotation` | `genes` |  | Genes |
| `annotation` | `transcripts`, `genes_with_introns` |  |  |
| `annotation` | `proteins_busco_complete` | % | BUSCO, proteins (%) |
| `expression` | `libraries`, `pairs_on_genes` | pairs |  |

And its tables: `candidates` (each assembly candidate's `candidate`, `kept`,
`busco_complete`, `qv`, `contigs`, `contig_n50`, `total_length`), `libraries`
(the species' `pairs` and `share` of each pooled library), `chromosomes` (each
sequence's `name`, `length_bp`, `role`, `original_name`) and `expression`
(each library's `pairs_on_genes`).

### `<sample>.report_mqc.yaml` { #report_mqcyaml }

[MultiQC](https://multiqc.info)'s custom content (`id: plexembler_report`):
the species' row of its General Statistics table, the columns above.

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
