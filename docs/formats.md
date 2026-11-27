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

## The closest genomes { #closest-genomes }

What [`plexembler closest`](cli.md#closest) writes for a species,
`<species>.closest.tsv`: tab-separated, with a header, a row per genome of the
database that sourmash gather matched to the species' genome, best first
(`--top` of them, 5): the first is its closest known genome. A species with no
match has the header alone.

| Column | Description |
|---|---|
| `rank` | The match's rank, from 1 |
| `accession` | Its NCBI assembly accession (`GCA_000146045.2`), from its name; `NA` when the name has none |
| `name` | Its name in the database (`GCA_000146045.2 Saccharomyces cerevisiae S288C`) |
| `ani` | Its average nucleotide identity to the species' genome, in percent, estimated from the k-mers they share (sourmash's `max_containment_ani`) |
| `species_containment` | The share of the species' k-mers found in it, in percent |
| `match_containment` | The share of its k-mers found in the species' genome, in percent |
| `shared_bp` | The bases they share, estimated from their k-mers |
| `lineage` | Its lineage, its ranks joined by `;`, from the database's taxonomy (`--lineages`); `NA` without it |
| `sequence` | What became of the first match's sequence: `downloaded`; `offline` (NCBI could not be reached); `not_found` (NCBI does not know the accession); `failed` (any other answer); `not_ncbi` (its name has no accession); `not_downloaded` (`--no-download`). `NA` for the other matches |

## The dotplot's table { #dotplot-table }

What [`plexembler dotplot`](cli.md#dotplot) writes beside the dotplot,
`<species>.dotplot.tsv`: tab-separated, with a header, a row per sequence of
the species' genome, in its order, then a row per sequence of the other genome,
in its order. The numbers are those of minimap2's primary alignments.

| Column | Description |
|---|---|
| `genome` | `species` or `reference`, the other genome |
| `sequence` | The sequence's name in its genome |
| `label` | Its name on the dotplot: the species' without `<species>_`; the reference's from its FASTA header (`chrIV` for `… chromosome IV, complete sequence`, `chrM` for a mitochondrion) |
| `length_bp` | Its length, in bases |
| `aligned_bp` | Its bases aligned to the other genome, each counted once |
| `aligned_fraction` | Those bases over its length, in percent |
| `identity` | The alignments' matching bases over their columns, in percent; `NA` without one |
| `best_match` | The sequence of the other genome it aligns to most |
| `best_match_share` | The share of its aligned bases that align to that one, in percent: under 100, it aligns to several (a translocation, a fusion) |
| `strand` | `+` when most of what aligns to `best_match` is on the same strand, `-` on the other |

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
| `closest_genome` | `matches`, `accession`, `name`, `lineage`, `sequence` |  |  |
| `closest_genome` | `ani`, `species_containment`, `match_containment` | % |  |
| `closest_genome` | `species_aligned`, `reference_aligned`, `identity` | % |  |

And its tables: `candidates` (each assembly candidate's `candidate`, `kept`,
`busco_complete`, `qv`, `contigs`, `contig_n50`, `total_length`), `libraries`
(the species' `pairs` and `share` of each pooled library), `chromosomes` (each
sequence's `name`, `length_bp`, `role`, `original_name`), `expression` (each
library's `pairs_on_genes`), `closest_genomes` (each match's `rank`,
`accession`, `name`, `ani`, `species_containment`, `match_containment`,
`lineage`) and `alignment_to_closest` (the dotplot's table: each sequence's
`genome`, `sequence`, `label`, `length_bp`, `aligned_pct`, `identity`,
`best_match`, `best_match_share`, `strand`). The `closest_genome` section's
`species_aligned` and `reference_aligned` are the shares of each genome aligned
to the other, and its `identity` the species' sequences' identity, weighed by
their aligned bases.

### `<sample>.report_mqc.yaml` { #report_mqcyaml }

[MultiQC](https://multiqc.info)'s custom content (`id: plexembler_report`):
the species' row of its General Statistics table, the columns above.

### `<sample>.report.html` { #reporthtml }

A page to read the run with. Its plots are drawn by plotly.js, which it loads
from its CDN: offline, it keeps its numbers and tables. With the species'
closest known genome, its dotplot is one of them.

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
