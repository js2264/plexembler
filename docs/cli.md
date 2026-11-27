# Command line

Every command of plexembler, and every option. `tests/test_docs.py`
checks this page against the command line: a command or an option that is
missing here, or that the command line does not have, fails the tests.

Each command has a section headed with its name, `` ### `name` { #name } ``
(`` #### `group name` `` for a command of a group),
with a table whose first column lists its options, short form first:
`` | `-o, --outdir` | … | ``. `-h, --help` and `--version` are left out.

## Commands

### `run` { #run }

Runs the pipeline with Nextflow, at plexembler's own version: `plexembler run
ARGS…` runs `nextflow run js2264/plexembler -r v<version> ARGS…`, and every
argument but its own options goes to Nextflow as it is: Nextflow's options
(`-profile`, `-resume`, `-c`, `-r`) and the pipeline's
[parameters](parameters.md) (`--input`, `--libraries`, `--outdir`).
plexembler brings Nextflow's launcher with it, PyPI's `nextflow`, which
fetches Nextflow at its first run ([installation](installation.md#with-uv)):
with Java, plexembler is all a machine needs to run the pipeline, but a
container engine.

```bash
plexembler run -profile docker \
    --input species.csv --libraries libraries.csv --outdir results
```

- **The pipeline, and its revision.** A release of plexembler (PyPI's, as
  `uvx plexembler` takes it) runs `js2264/plexembler` at its tag,
  `v<version>`: the package's version is the pipeline's. Installed from a folder that holds the
  pipeline (a clone, by `uv sync` or `pip install -e .`), it runs that folder
  as it is; installed from the git repository by pip or uv, that repository at
  the commit installed. `-r` runs another revision of the repository,
  `--pipeline` another pipeline.
- **Nextflow** is the one beside plexembler's Python, the launcher that uv
  and pip install with it (conda, its Nextflow), else the one on `PATH`:
  26.04.0 or newer, which needs Java 17 or newer, the machine's, or jdk4py's
  with the `java` extra (`plexembler[java]`), which plexembler gives it
  (`NXF_JAVA_HOME`, unless set). plexembler hands
  its process over to it, so that Nextflow's exit status and logs are the
  run's (`.nextflow.log` where it runs, `<outdir>/pipeline_info/`), and Ctrl-C
  stops Nextflow as it would without plexembler.
- **Its own options** are read word for word, wherever they are; every
  argument after `--` goes to Nextflow, so that `plexembler run -- --help`
  lists the pipeline's parameters, which `--help` alone does not.
- **No log of its own,** and none of the options of the commands that write
  outputs (`-o`, `--force`, `-p`, `-v`, `-q`): the run's outputs and log are
  Nextflow's, and a single-dash option is Nextflow's, as `-profile` is.

`--dry-run` prints the command:

```console
$ plexembler run --dry-run -profile docker --input species.csv --outdir results
/usr/local/bin/nextflow run js2264/plexembler -r v0.1.0 -profile docker --input species.csv --outdir results
```

| Option | Default | Description |
|---|---|---|
| `NEXTFLOW_ARGS` | | Nextflow's options and the pipeline's parameters, as `nextflow run` takes them. |
| `--pipeline` | the folder or the commit plexembler was installed from, else `js2264/plexembler` | The pipeline: a folder (a clone), run as it is, or a repository, which Nextflow pulls, run at `v<version>` unless `-r` names a revision. |
| `--dry-run` | off | Print the Nextflow command, and run nothing. |

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

**A circle written twice.** An assembler can write a circular genome, a
mitochondrion or a plasmid, as two or more copies of itself joined end to end:
on some runs, Flye wrote the benchmark's 20,063-bp *N. glabratus*
mitochondrion as one 40,128-bp contig, which QUAST reads as a misassembly. A
contig that is such copies is cut to its first copy, before anything else sees
it, and the log says so: its copies, the unit's length and how alike they are.
A contig is `k` copies of a unit when its 21-mers, but the last copy's, are
found again one unit on (within as many bases as the copies may differ by),
at `--min-identity` or more (the share found is about the identity to the
21st power), and when it is a whole number of units long: a circle's copies
are, the last one closing on the first, and a tandem repeat of a linear
chromosome seldom is. The unit is `--min-unit` bases or more, and a contig that
is a tandem repeat of a shorter one throughout, a satellite, stays whole. With
`--assembly-info`, Flye's table of its contigs, a contig Flye says is linear
(`circ.` N), a tandem repeat such as the rDNA's, is never cut; its circular
ones (`circ.` Y, a doubled one with a `mult.` of 2 or more) and every contig of
an assembler without such a table (hifiasm) are judged by their sequence.

```bash
plexembler rename assembly.fasta -p Cgla -o Cgla.contigs.fa -a assembly_info.txt
```

| Option | Default | Description |
|---|---|---|
| `FASTA` | (required) | The contigs, as FASTA; gzipped or not. |
| `-p, --prefix` | (required) | The species' ID, which every new name starts with: a letter, then letters, digits and hyphens. It names the log. |
| `-o, --output` | (required) | The renamed contigs, as FASTA; gzipped when the name ends in `.gz`. |
| `-n, --names` | `<output>.names.tsv` | The table of new and original names; by default beside the output, its `.fa`/`.fasta`/`.fna` and `.gz` extensions replaced. |
| `-a, --assembly-info` | | The assembler's table of its contigs, Flye's `assembly_info.txt`: a contig it says is linear (`circ.` N) is never cut. |
| `--collapse/--no-collapse` | on | Cut a contig that is two or more copies of a unit, end to end (a circle its assembler went round twice), to its first copy. |
| `--min-unit` | `1000` | The shortest unit, in bases, that a contig is looked at as copies of. |
| `--min-identity` | `0.99` | The lowest identity between a contig's copies, from their 21-mers, at which it is cut. |
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

The polishing step's candidates are the draft, first, then each polishing.
Where a polishing changed the draft, its `kmers` count the k-mers that reads
held out of the polishing lack: the draft's that it removed (fixed) and its own
(introduced). The polishing kept is the one that fixed the most more than it
introduced, if one did; otherwise the draft. Only the changed k-mers count, so
that a k-mer the held-out reads miss by chance, where no polishing touched the
draft, cannot decide. Without `kmers` for every polishing, one is kept when its
QV beats the draft's by more than `--qv-tolerance`.

```bash
plexembler select candidates.tsv -p Scer -o Scer/
```

`candidates.tsv` has a row per candidate: its name (`candidate`), its sequence
(`fasta`), and its BUSCO summary (`busco`: BUSCO's short summary, `.txt` or
`.json`, or compleasm's `summary.txt`), its Merqury QV (`qv`) and its gfastats
summary (`stats`), and, for a polishing, the k-mers it changed (`kmers`); a
check may be empty or `NA`, and a relative path is read from the table's folder
([format](formats.md#candidates-table)).

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
| `--qv-tolerance` | `1.0` | QVs within this many units of the best are a tie; a polishing without k-mer counts must beat the draft's QV by more. |
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
depth on its own contigs, those its reads cover the most: a whole genome of
bleed-through, assembled, can outweigh a species' own), the share of it found
in another species' contigs at 99 % identity or more (and in a longer contig of
its own species at 90 %), and its Hi-C contacts with each species' contigs.

- `organelle`: `--organelle-depth` times its species' depth and under
  `--organelle-gc` GC: kept as `<species>_chrM` (`chrM2`, … for more);
- `drop`: most of it (`--duplicate-share`) found in another species' contigs, its own
  species' reads thin on it (under `--own-depth` of their depth), or that species' reads
  `--depth-ratio` times as deep on it, each against its own species' depth: that
  species' DNA, which it holds already. Or, for a species given `--haploid`,
  `--own-duplicate-share` of it found in a longer contig of its own, 90 %
  identical or more, and its own reads under `--own-duplicate-depth` of their
  depth: a stretch of the genome its contigs hold already, which the assembler
  wrote twice (a repeat of its graph, in a contig and alone), and to which the
  reads go less, as they go to the better copy. A diploid's such contigs are
  its haplotigs, and a hybrid's its other subgenome: the screen leaves them;
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
species' reads, share found in each species' contigs (in its own species',
the share in its longer contigs), contacts with each species, verdict, new
name and reason), their MultiQC table and its log.

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
| `--haploid` | | A haploid species of the pool, whose contigs another of its contigs holds already go (`--own-duplicate-share`); once per species. |
| `--own-duplicate-share` | `0.5` | The share of a contig found in a longer contig of its own species, 90 % identical or more, from which it is a duplicate of it. |
| `--own-duplicate-depth` | `0.75` | Its own species' depth on a duplicate of one of its contigs, as a share of their depth, under which it goes. |
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

With `--cut`, the table of the scaffolder's joins that
[`junctions --contigs`](#junctions) scored, it first cuts each join that table
finds `unsupported`: the scaffold's two sides become two scaffolds,
`scaffold_1.1` and `scaffold_1.2` (and so on, for more cuts in one scaffold),
the gap between them dropped, and are named as any other. The log gives each
cut, its pairs across and those expected. A `thin` join is kept.

```bash
plexembler finalize yahs_scaffolds_final.fa -a yahs_scaffolds_final.agp --cut Scer.junctions.tsv -p Scer -n 16 -o Scer/
```

| Option | Default | Description |
|---|---|---|
| `FASTA` | (required) | The scaffolds, FASTA, gzipped or not. |
| `-p, --prefix` | (required) | The species' ID, which starts every name and names the outputs. |
| `-o, --outdir` | (required) | The folder of the genome, its sizes, names, AGP and the log. |
| `-a, --agp` | | The scaffolder's AGP, written again with the new names. |
| `-n, --n-chromosomes` | | The species' number of nuclear chromosomes: its longest scaffolds, named chromosomes. |
| `--min-chromosome` | `100000` | Without `--n-chromosomes`, the shortest scaffold named a chromosome, in bases. |
| `--cut` | | The table of the scaffolder's joins from `plexembler junctions --contigs`: each it finds unsupported is cut first, its two sides two scaffolds (needs `--agp`). |
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

Both counts leave out the positions the pairs cannot be trusted on: each side
of a join is cut into bins of 1 kb, and a bin counts when its mates' coverage
is half to twice the genome's mean. A contig often ends in a repeat: where its
copies all assembled, no mate maps there uniquely, and a true join would look
unsupported; where the assembly collapsed them, or misses a copy, another
copy's reads map there, and their contacts can make a misjoin look supported,
which is how a scaffolder makes one. A join with no such bin on a side is
`thin`.

```bash
plexembler junctions Scer.pairs.gz -a Scer.agp -p Scer -o Scer/
```

The pairs are the species' deduplicated pairs on its genome (4DN pairs,
pairtools'), and the AGP the genome's, in the same names
([`finalize`](#finalize)'s). Or, with `--contigs`, its pairs on the contigs the
scaffolder laid out, and the scaffolder's AGP: each mate is placed on its
scaffold through the AGP (a contig broken into pieces, or reversed, included),
which scores the scaffolder's joins before the genome is written, so that
[`finalize --cut`](#finalize) can cut those the pairs do not support. A mate
on a part of a contig the AGP leaves out is dropped. It writes the table of joins
(`Scer.junctions.tsv`: `scaffold`, `position`, the last base of the join's left
contig, `left_contig`, `right_contig`, `observed`, `expected`, `ratio`,
`verdict`, and `left_seen` and `right_seen`, the share of each side's window
whose pairs count), its MultiQC table and its log.

| Option | Default | Description |
|---|---|---|
| `PAIRS` | (required) | The species' pairs on its genome, 4DN pairs, gzipped or not. |
| `-a, --agp` | (required) | The genome's AGP, in the names of the pairs' chromosomes. |
| `-p, --prefix` | (required) | The species' ID, which names the outputs. |
| `-o, --outdir` | (required) | The folder of the table of joins, its MultiQC table and the log. |
| `-w, --window` | `20000` | How far from a join, on each side, a pair's mates may lie, in bases. |
| `--min-ratio` | `0.2` | The pairs across a join over those expected, below which it is unsupported. |
| `--min-expected` | `10.0` | The pairs a join is expected to have, below which it is too thin to judge. |
| `--contigs` | off | The pairs are on the AGP's contigs: each mate is placed on its scaffold through the AGP (a scaffolder's joins, before the genome is written). |
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

### `expression` { #expression }

Writes a species' expression tables from featureCounts' table of its read
pairs of each pooled RNA-seq library ([`split --type rnaseq`](#split) gave it
its own pairs alone): the counts, genes × libraries, after each gene's length
(its exons', without overlaps, as featureCounts measures it); and the TPM,
within the species: each gene's pairs over its length, as a share of the
library's sum, times a million.

```bash
plexembler expression Scer.featureCounts.tsv -p Scer -o Scer/ -l rna-ypd-1=YPD_1 -l rna-gal-1=GAL_1
```

It writes `Scer.counts.tsv` (`gene`, `length_bp`, a column per library) and
`Scer.tpm.tsv` (`gene`, a column per library), each library's column named by
`--label`, `<condition>_<replicate>` in the pipeline; its MultiQC table (each
library's pairs on genes and genes detected) and its log.

| Option | Default | Description |
|---|---|---|
| `COUNTS` | (required) | featureCounts' table, a column per library's BAM (`<species>.<library>.bam`). |
| `-p, --prefix` | (required) | The species' ID, which starts its alignments' names and names the outputs. |
| `-o, --outdir` | (required) | The folder of the tables and the log. |
| `-l, --label` | | A library's column name, LIBRARY=LABEL (rna-ypd-1=YPD_1); the library's ID otherwise. |
| `--min-count` | `10` | The fewest pairs that make a gene detected, in MultiQC's table. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `closest` { #closest }

Names a species' closest known genome from `sourmash gather`'s matches of its
genome in a database of genomes (sourmash's prepared databases of NCBI's
genomes: [usage](usage.md#the-closest-known-genome)), best first; with the
database's taxonomy (`--lineages`), each match's lineage. Then downloads the first
match's sequence from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/),
by the assembly accession its name starts with (`GCA_000146045.2 Saccharomyces
cerevisiae S288C`): the reference of the species' [`dotplot`](#dotplot).

The download is the command's only use of the network, and it never fails the
command: without network, for an accession NCBI does not know, or a match
whose name has no accession, it warns, writes the table without the sequence
and ends without an error. The species then has no dotplot.

```bash
plexembler closest Scer.gather.csv.gz -p Scer -o Scer/ -l ncbi-eukaryotes.2025.01.lineages.csv
```

It writes the matches (`Scer.closest.tsv`: [format](formats.md#closest-genomes)),
the first one's sequence when it was downloaded (`Scer.closest.fa.gz`, NCBI's
FASTA, gzipped), MultiQC's table (`Scer.closest_mqc.tsv`) and its log. An empty
gather file (`sourmash gather --create-empty-results`) is no match: the table
has its header alone.

| Option | Default | Description |
|---|---|---|
| `GATHER` | (required) | sourmash gather's matches (CSV, gzipped or not). |
| `-p, --prefix` | (required) | The species' ID, which names the outputs. |
| `-o, --outdir` | (required) | The folder of the table, the closest genome's sequence, MultiQC's table and the log. |
| `-l, --lineages` | | The database's taxonomy, sourmash's lineages (CSV, or gzipped CSV): each match's lineage. |
| `-n, --top` | `5` | How many matches the table keeps, best first. |
| `--download/--no-download` | `--download` | Download the closest genome's sequence from NCBI, for the dotplot. |
| `--timeout` | `60.0` | How long to wait for NCBI at each step of the download, in seconds. |
| `--attempts` | `3` | How many times to try NCBI when it cannot be reached or is busy. |
| `--api` | `https://api.ncbi.nlm.nih.gov/datasets/v2` | NCBI Datasets' API, which the sequence is downloaded from. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `dotplot` { #dotplot }

Draws a species' genome against another, its closest known genome
([`closest`](#closest)), from minimap2's alignments of the first (the query) to
the second (the target): the other genome's sequences along x, in their order,
named from their FASTA headers (`chrIV` for NCBI's `… chromosome IV, complete
sequence`, `chrM` for a mitochondrion); the species' along y, each beside the
sequence it aligns to most, turned round, and marked (−), when most of it
aligns to that one's other strand. The chromosomes the two genomes share make a
diagonal; an inversion, a translocation or a fusion stands out of it.
Alignments on the same strand are blue, on the other vermilion.

```bash
minimap2 -x asm20 -c Scer.closest.fa.gz Scer.fa.gz > Scer.closest.paf
plexembler dotplot Scer.closest.paf -p Scer -o Scer/ -g Scer.fa.gz -r Scer.closest.fa.gz -c Scer.closest.tsv
```

It writes the dotplot (`Scer.dotplot.svg`), the same as a plotly figure, which
the species' [report](#report) draws (`Scer.dotplot.json`), the table of every
sequence of both genomes and how it aligns to the other
(`Scer.dotplot.tsv`: [format](formats.md#dotplot-table)), and its log. The
table counts every primary alignment; the plot draws those of `--min-length`
bases or more, secondary ones included.

| Option | Default | Description |
|---|---|---|
| `PAF` | (required) | minimap2's alignments of the species' genome (query) to the other (target). |
| `-p, --prefix` | (required) | The species' ID, which names the outputs and the y axis. |
| `-o, --outdir` | (required) | The folder of the dotplot, its figure, its table and the log. |
| `-g, --genome` | | The species' genome (FASTA), the alignments' query: its sequences without an alignment are drawn too. |
| `-r, --reference` | | The other genome (FASTA), the alignments' target: its sequences' order, and names from their headers. |
| `-c, --closest` | | plexembler closest's table, whose first match names the other genome on the x axis. |
| `--min-length` | `1000` | The shortest alignment drawn, in bases; the table counts every one. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `report` { #report }

Writes a species' report from the files the pipeline wrote for it, each known
by its name: its assembly's candidates and choice
(`<species>.assembly.selection.tsv`), its polishing's
(`<species>.polish.selection.tsv`), its share of each pooled library
(`<pool>.<library>.split.tsv`), the screen's verdicts on its contigs
(`<pool>.screen.tsv`), its chromosomes (`<species>.names.tsv`), BUSCO and
Merqury on its genome (`short_summary.*.<species>.fa.txt`,
`<species>.genome.<k-mers>.qv`), its joins (`<species>.junctions.tsv`), its
genes (`<species>.genes.tsv`), BUSCO on its proteins
(`short_summary.*.<species>.proteins.faa.txt`), its pairs on genes
(`<species>.counts.tsv`), its closest known genome (`<species>.closest.tsv`)
and how it aligns to it (`<species>.dotplot.tsv`, and the dotplot,
`<species>.dotplot.json`, which the page draws). A file it is not given leaves
its numbers `NA`.

```bash
plexembler report Scer.* poolA.* short_summary.* -p Scer --pool poolA -o Scer/
```

It writes chromokit's three files ([format](formats.md#the-report-of-a-species)):
`Scer.report.txt`, every number; `Scer.report_mqc.yaml`, the species' row of
MultiQC's General Statistics table (its genome's length, chromosomes, BUSCO
and QV, its genes, BUSCO on its proteins); `Scer.report.html`, its page; and
its log.

| Option | Default | Description |
|---|---|---|
| `FILES` | (required) | The species' files, and its pool's. |
| `-p, --prefix` | (required) | The species' ID, which its files' names start with and which names the report. |
| `--pool` | | The species' pool, which starts its pool's files' names. |
| `-o, --outdir` | (required) | The folder of the report and the log. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |

### `browser` { #browser }

Writes `chromosomix.store.json`, [chromosomix](https://chromosomix.app)'s
store of the results: a row per species' genome and Hi-C map, each with its
tracks (genes, RNA-seq coverage, nanopore depth) placed, laid out as cards,
one per species, every file named by its path in the results. They then open
in chromosomix as they are: its app with the results folder connected, or
chromosomix.app with `?store=<url>` once they are on a bucket (`--base`).
Each file is known by its name: a `.mcool` is the species' map, a FASTA its
genome, a GFF3, GTF or bigWig one of its tracks.

```bash
plexembler browser -s Scer=poolA -f Scer=species/Scer/hic/Scer.mcool -f Scer=species/Scer/genome/Scer.fa.gz -o results/
```

| Option | Default | Description |
|---|---|---|
| `-s, --species` | (required) | A species and its pool, SPECIES=POOL; once per species. |
| `-f, --file` | | A species' map, genome or track, by its path in the results; once per file. |
| `-o, --outdir` | (required) | The top of the results, where the store goes, and the log's folder. |
| `-n, --name` | `plexembler` | The store's name, as chromosomix's picker shows it. |
| `--base` | | Where the results are served (https://…/results/), when they are on a bucket. |
| `--force` | off | Overwrite outputs that exist already. |
| `-v, --verbose` | off | Debug messages, in the log and on stderr. |
| `-q, --quiet` | off | Only warnings and errors on stderr; the log keeps everything. |
