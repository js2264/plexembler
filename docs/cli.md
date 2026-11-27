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

### `split` { #split }

Splits a pooled library between the species of its pool, from its reads'
alignments to every species' contigs at once (bwa-mem2; `-5SP` for Hi-C), as
bwa-mem2 writes them: a pair's records together. A pair goes to a species when
both mates' primary alignments are in it, unique (`--min-mapq`) and
near-identical (`--min-identity`, one minus the edit distance over the
alignment's columns): a read of a region missing from its own species' contigs
aligns, unique, to a sister species' copy, and only its identity tells. A pair
no species gets is counted: ambiguous (between which species, from bwa-mem2's
alternative hits), divergent, between species, or unmapped.

Each library keeps `--holdout` of its pairs apart, by a hash of the read's
name: the k-mers the polishing's QV guard counts against. A Hi-C library loses
its duplicate pairs; its unique pairs between two contigs are counted, the
contig screen's Hi-C signal; and its pairs between species estimate the
random-ligation rate, epsilon: their share over 1 - the sum of the species'
squared shares.

```bash
samtools view -h poolA.wgs.bam | plexembler split - -p poolA -o split/ -l wgs -t wgs
```

It writes each species' pairs, `<species>.<library>_R1.fastq.gz` and `_R2`,
and those held out, `<species>.<library>.holdout_R1.fastq.gz` and `_R2`, the
reads as they were sequenced; for Hi-C, the contacts between contigs
(`poolA.hic.contacts.tsv`: `contig_a`, `contig_b`, `pairs`); the split's table (`poolA.wgs.split.tsv`: each species' pairs,
share, held out and duplicates), its numbers (`poolA.wgs.split.json`: with the
ambiguous pairs of each pair of species, the pairs between species, and
epsilon), its MultiQC table and its log.

An RNA-seq library (`-t rnaseq`) is aligned by STAR to the pool's genomes,
whose sequences are each species' own: each species' pairs are written as the
alignments they are, `<species>.<library>.sam.gz`, its header holding that
species' sequences alone, ready to sort and count. STAR writes no alternative
hits: a multimapper's secondary records say between which species it is
ambiguous. RNA-seq keeps its duplicates and holds no pair out.

```bash
samtools view -h poolA.rna-ypd-1.bam | plexembler split - -p poolA -o split/ -l rna-ypd-1 -t rnaseq
```

| Option | Default | Description |
|---|---|---|
| `SAM` | (required) | The alignments, SAM, gzipped or not; `-` for stdin. |
| `-p, --prefix` | (required) | The pool's ID, which names the tables and the log. |
| `-o, --outdir` | (required) | The folder of each species' FASTQ (or SAM) files, the tables and the log. |
| `-l, --library` | (required) | The library's ID, which names each species' FASTQ (or SAM) files. |
| `-t, --type` | (required) | The library's type: wgs; hic, which drops duplicates, counts the contacts between contigs and estimates epsilon; or rnaseq, STAR's alignments, written as each species' SAM. |
| `--min-mapq` | `20` | A mate's lowest mapping quality: below, it is ambiguous. |
| `--min-identity` | `0.98` | A mate's lowest identity to its contig: below, it is divergent, a read of another species. |
| `--holdout` | `0.1` | The share of the library's pairs kept apart, for the QV guard. |
| `--dedup/--no-dedup` | on | Drop a Hi-C library's duplicate pairs. |
| `--compress-level` | `1` | The gzip level of the FASTQ files. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `digest` { #digest }

Cuts a species' Hi-C reads at their ligation junctions, into single reads to
polish its contigs with when its pool has no WGS library (see
[the plan](plan.md#m2-the-pool-split-v020)). A junction is the filled-in end of
one restriction fragment ligated to another's: for two of the library's
enzymes (an enzyme with itself included), the first's site up to the far end
of its overhang, then the second's from the near end of it (DpnII's GATCGATC,
HinfI's GANTANTC; four for Arima's kit). Each read is cut where a junction's
first end stops, as hicstuff does, so that each piece keeps its half of the
site; the pieces shorter than `--min-piece` are dropped; and every piece is
written once, as a single read named `<read>_<mate>.<piece>`.

```bash
plexembler digest Scer.hic_R1.fastq.gz Scer.hic_R2.fastq.gz -e DpnII,HinfI -p Scer -o Scer/
```

It writes the pieces (`Scer.digested.fastq.gz`), the numbers (`Scer.digest.json`:
the reads with no junction, one, and two or more, the pieces kept and dropped),
their MultiQC table and its log.

| Option | Default | Description |
|---|---|---|
| `FASTQS` | (required) | The Hi-C reads, FASTQ, gzipped or not: the first and the second reads, or one of them. |
| `-e, --enzymes` | (required) | The library's restriction enzymes, comma-separated (`DpnII,HinfI`), or `arima`. |
| `-p, --prefix` | (required) | The species' ID, which names the outputs and the log. |
| `-o, --outdir` | (required) | The folder of the pieces, their numbers and the log. |
| `--min-piece` | `40` | The shortest piece kept, in bases. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `screen` { #screen }

Keeps, moves or drops each contig of a pool, and names its organelles. A
species' nanopore run can hold another species' reads (bleed-through), which
assemble into contigs of the wrong species; the screen weighs three signals
for each contig: each species' nanopore depth on it (against that species'
depth on its own contigs), the share of it found in another species' contigs
at 99 % identity or more, and its Hi-C contacts with each species' contigs.

- `organelle`: `--organelle-depth` times its species' depth and under
  `--organelle-gc` GC: kept as `<species>_chrM` (`chrM2`, … for more);
- `drop`: most of it (`--duplicate-share`) found in another species' contigs, its own
  species' reads thin on it (under `--own-depth` of their depth), or that species' reads
  `--depth-ratio` times as deep on it, each against its own species' depth: that
  species' DNA, which it holds already;
- `move`: its own species' reads under `--own-depth` of their depth, another
  species' at `--other-depth` of theirs; or its own species' reads under
  `--hic-own-depth` and its Hi-C contacts with another species over
  `--hic-ratio` times its own species' (and `--min-contacts` of them): renamed
  after the species it goes to, its name kept (`Cgla_Scer_contig_5`). Hi-C alone
  moves no contig its own reads cover well: a species of one contig has no
  contacts of its own, and random ligation gives every contig some with the
  others;
- `keep`: the others.

```bash
plexembler screen poolA.fa.gz -r Scer=Scer.paf -r Kmar=Kmar.paf -i poolA.self.paf \
    -c poolA.hic.contacts.tsv -p poolA -o screen/
```

It writes each species' screened contigs (`<species>.screened.fa.gz`), the
verdicts (`poolA.screen.tsv`: each contig's length, GC, depth of each
species' reads, share found in each species' contigs, contacts with each
species, verdict, new name and reason), their MultiQC table and its log.

| Option | Default | Description |
|---|---|---|
| `CONTIGS` | (required) | Every species' contigs of the pool, named `<species>_<contig>`, FASTA. |
| `-r, --reads` | (required) | A species' nanopore reads aligned to the pool's contigs, `SPECIES=PAF`; once per species. |
| `-i, --identity` | | The pool's contigs aligned to each other (`minimap2 -x asm10 -X`), PAF. |
| `-c, --contacts` | | A Hi-C library's contacts between contigs ([`plexembler split`](#split)); once per library. |
| `-p, --prefix` | (required) | The pool's ID, which names the table and the log. |
| `-o, --outdir` | (required) | The folder of each species' screened contigs, the table and the log. |
| `--organelle-depth` | `2.0` | An organelle's depth, in times its species' contigs'. |
| `--organelle-gc` | `0.3` | An organelle's highest GC share. |
| `--duplicate-share` | `0.9` | The share of a contig found in another species' contigs that makes it theirs. |
| `--own-depth` | `0.2` | Its own species' reads' depth, in times their usual depth, under which they are thin on it. |
| `--other-depth` | `0.5` | Another species' reads' depth, in times their usual depth, at which they cover it. |
| `--hic-ratio` | `2.0` | How many times its Hi-C contacts with another species outnumber its own species'. |
| `--min-contacts` | `20` | The fewest Hi-C contacts with another species that move a contig. |
| `--hic-own-depth` | `0.5` | Its own species' depth on a contig, as a share of their depth, under which Hi-C contacts can move it. |
| `--depth-ratio` | `2.0` | Another species' depth on a duplicate, over its own species', each against their depth, from which the duplicate goes. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `finalize` { #finalize }

Names a species' scaffolds, as its genome is read: its chromosomes, longest
first (`<species>_chr01`, …), the `--n-chromosomes` longest scaffolds, or every
scaffold of `--min-chromosome` bases or more; its mitochondrion,
`<species>_chrM`, a scaffold made of the contig screen's `<species>_chrM`
alone (found through the AGP when the scaffolder renamed it); and the others,
longest first, `<species>_unplaced_01`, ….

```bash
plexembler finalize yahs_scaffolds_final.fa -a yahs_scaffolds_final.agp -p Scer -n 16 -o Scer/
```

It writes the genome (`Scer.fa.gz`), its sizes (`Scer.chrom.sizes`), the table
of names (`Scer.names.tsv`: `name`, `original_name`, `length_bp` and `role`,
chromosome, organelle or unplaced), the AGP with the new names (`Scer.agp`)
and its log.

| Option | Default | Description |
|---|---|---|
| `FASTA` | (required) | The scaffolds, FASTA, gzipped or not. |
| `-p, --prefix` | (required) | The species' ID, which starts every name and names the outputs. |
| `-o, --outdir` | (required) | The folder of the genome, its sizes, names, AGP and the log. |
| `-a, --agp` | | The scaffolder's AGP, written again with the new names. |
| `-n, --n-chromosomes` | | The species' number of nuclear chromosomes: its longest scaffolds, named chromosomes. |
| `--min-chromosome` | `100000` | Without `--n-chromosomes`, the shortest scaffold named a chromosome, in bases. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `junctions` { #junctions }

Scores each join of a genome's scaffolds by its Hi-C pairs: the species'
contact decay P(s), from the pairs whose mates lie in one contig, where no join
can mislead it; for each join, the pairs across it, a mate within `--window`
bases on each side, against those the decay expects there. A join is
`supported` at `--min-ratio` of them or more, `unsupported` below (a misjoin,
likely: the log warns of each), `thin` when the decay expects too few pairs to
judge it (`--min-expected`).

```bash
plexembler junctions Scer.pairs.gz -a Scer.agp -p Scer -o Scer/
```

The pairs are the species' deduplicated pairs on its genome (4DN pairs,
pairtools'), and the AGP the genome's, in the same names
([`finalize`](#finalize)'s). It writes the table of joins
(`Scer.junctions.tsv`: `scaffold`, `position`, the last base of the join's left
contig, `left_contig`, `right_contig`, `observed`, `expected`, `ratio` and
`verdict`), its MultiQC table and its log.

| Option | Default | Description |
|---|---|---|
| `PAIRS` | (required) | The species' pairs on its genome, 4DN pairs, gzipped or not. |
| `-a, --agp` | (required) | The genome's AGP, in the names of the pairs' chromosomes. |
| `-p, --prefix` | (required) | The species' ID, which names the outputs. |
| `-o, --outdir` | (required) | The folder of the table of joins, its MultiQC table and the log. |
| `-w, --window` | `20000` | How far from a join, on each side, a pair's mates may lie, in bases. |
| `--min-ratio` | `0.2` | The pairs across a join over those expected, below which it is unsupported. |
| `--min-expected` | `10.0` | The pairs a join is expected to have, below which it is too thin to judge. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `genes` { #genes }

Names a species' gene models and writes its annotation. It reads a predictor's
models, GTF (BRAKER's) or GFF3 (AUGUSTUS's), and names the genes in the
genome's order, `<species>_000010`, `<species>_000020`, … (a step of 10, room
for genes found later), their transcripts `<species>_000010-T1`, …. The stop
codon goes in the CDS.

```bash
plexembler genes Scer.softmasked.fa.gz braker.gtf -p Scer -o Scer/ -g 1 -s BRAKER3
```

It writes the annotation as GFF3 and GTF (`Scer.gff3.gz`, `Scer.gtf.gz`: each
gene, transcript, exon and CDS), the proteins (`Scer.proteins.faa.gz`),
translated with the species' genetic code, and the CDS (`Scer.cds.fna.gz`); a
table of the genes (`Scer.genes.tsv`: `gene`, `transcripts`, `chromosome`,
`start`, `end`, `strand`, `exons`, `protein_length` and `original_name`), its
MultiQC table and its log. A protein with a stop codon inside is warned about:
the genetic code may be wrong.

| Option | Default | Description |
|---|---|---|
| `GENOME` | (required) | The genome the models lie on, FASTA, gzipped or not. |
| `MODELS` | (required) | The gene models, GTF or GFF3. |
| `-p, --prefix` | (required) | The species' ID, which starts every gene's name and names the outputs. |
| `-o, --outdir` | (required) | The folder of the annotation and the log. |
| `-g, --genetic-code` | `1` | NCBI's translation table of the species' genes: 12 for the CUG-Ser clade, 26 for CUG-Ala. |
| `-s, --source` | `plexembler` | The GFF3 and GTF files' second column: what predicted the genes. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

