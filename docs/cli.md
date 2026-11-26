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
