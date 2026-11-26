# Command line

Every command of plexembler, and every option. `tests/test_docs.py`
checks this page against the command line: a command or an option that is
missing here, or that the command line does not have, fails the tests.

Each command has a section headed with its name, `` ### `name` { #name } ``
(`` #### `group name` `` for a command of a group),
with a table whose first column lists its options, short form first:
`` | `-o, --outdir` | … | ``. `-h, --help` and `--version` are left out.

## Commands

### `rename` { #rename }

Names every contig after its species, `<species>_<contig>`: the name of every
contig the pipeline makes starts with its species (see
[concepts](concepts.md#species-prefixed-contigs)). A contig whose name already
starts with `<species>_` keeps it. Sequences are written 80 bases a line.

```bash
plexembler rename contigs.fa -p Scer -o Scer.contigs.fa.gz
```

It writes the renamed FASTA, the table of names beside it (`Scer.contigs.names.tsv`:
[format](formats.md#names-table)), and its log, `Scer.plexembler_<YYYYMMDDHHMMSS>.log`.

| Option | Default | Description |
|---|---|---|
| `FASTA` | (required) | The contigs, as FASTA; gzipped or not. |
| `-p, --prefix` | (required) | The species' ID, which every new name starts with: a letter, then letters, digits and hyphens. It names the log. |
| `-o, --output` | (required) | The renamed contigs, as FASTA; gzipped when the name ends in `.gz`. |
| `-n, --names` | `<output>.names.tsv` | The table of new and original names; by default beside the output, its `.fa`/`.fasta`/`.fna` and `.gz` extensions replaced. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `select` { #select }

Keeps one candidate of a step for a species, by the candidates' checks, and
says why. The assembly step's rules, in order, each breaking the ties of the
one before: the most single-copy complete BUSCOs (complete minus duplicated;
within `--busco-tolerance` points is a tie), the highest QV (Merqury; within
`--qv-tolerance`), the longest contig N50. A number one candidate lacks leaves
its rule out for all of them; `--keep` keeps a candidate whatever its numbers,
which is what the samplesheet's `assembler` column does.

```bash
plexembler select candidates.tsv -p Scer -o Scer/
```

`candidates.tsv` has a row per candidate: its name (`candidate`), its sequence
(`fasta`), and its BUSCO summary (`busco`: BUSCO's short summary, `.txt` or
`.json`, or compleasm's `summary.txt`), its Merqury QV (`qv`) and its gfastats
summary (`stats`); a check may be empty or `NA`, and a relative path is read
from the table's folder ([format](formats.md#candidates-table)).

It writes the kept sequence (`Scer.assembly.fa.gz`), the table of candidates
(`Scer.assembly.selection.tsv`: [format](formats.md#selection-table)), its
MultiQC table (`Scer.assembly_selection_mqc.tsv`) and its log.

| Option | Default | Description |
|---|---|---|
| `CANDIDATES` | (required) | The table of candidates. |
| `-s, --step` | `assembly` | The step whose candidates these are, which sets the rules. |
| `-p, --prefix` | (required) | The species' ID, which names the outputs and the log. |
| `-o, --outdir` | (required) | The folder of the kept sequence, the table of candidates and the log. |
| `-k, --keep` | | A candidate to keep whatever its numbers: the samplesheet's choice. |
| `--busco-tolerance` | `0.5` | Single-copy BUSCOs within this many points of the best are a tie. |
| `--qv-tolerance` | `1.0` | QVs within this many units of the best are a tie. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `ani` { #ani }

The identity between the species of a pool, from skani's pairs
(`skani triangle --sparse`): a genome is named by its contigs' species. A pair
above `--warn-above` is too close for the pool's short reads to be split well
between its species, and is logged as a warning; a pair from `--genus-above`
is of the same genus, whose reads are split by identity as well as by
uniqueness (see [the plan](plan.md#m2-the-pool-split-v020)). skani reports no
pair below about 80 %: such a pair is `NA`, and distinct.

```bash
plexembler ani poolA.skani.tsv -p poolA -o poolA/ -s Scer -s Spar -s Kmar
```

It writes the matrix (`poolA.ani.tsv`), each pair with its aligned fractions
and verdict (`poolA.ani_pairs.tsv`: [format](formats.md#identity-pairs)), the
pairs for MultiQC (`poolA.ani_mqc.tsv`) and its log.

| Option | Default | Description |
|---|---|---|
| `SKANI` | (required) | skani's pairs, `skani triangle --sparse`. |
| `-p, --prefix` | (required) | The pool's ID, which names the outputs and the log. |
| `-o, --outdir` | (required) | The folder of the matrix, the pairs and the log. |
| `-s, --species` | | A species of the pool, listed even when skani paired it with none; once per species. |
| `--warn-above` | `95.0` | The identity, in percent, above which a pair is too close: a warning. |
| `--genus-above` | `85.0` | The identity, in percent, from which a pair is of the same genus. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |
