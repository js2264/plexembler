# Implementation plan

How plexembler gets from its scaffold to the pipeline that [concepts](concepts.md)
describes. The work is cut into **milestones**, each a release that runs on real
data, and **tracks** that run beside them. A box is ticked in the pull request
that delivers it; a decision taken along the way is written where it applies.

```mermaid
flowchart LR
  M0[M0 Foundations] --> M1[M1 Contigs<br>v0.1.0]
  M1 --> M2[M2 The pool, split<br>v0.2.0]
  M2 --> M3[M3 Chromosomes<br>v0.3.0]
  M3 --> M4[M4 Annotation<br>v0.4.0]
  M4 --> M5[M5 Expression<br>v0.5.0]
  M3 --> M6[M6 Reports and extras<br>v0.6.0]
  M5 --> M6
  M6 --> C[The complete workflow<br>a 0.x release]
  T2[T2 Test data, benchmark pool] -.-> M1
  T1[T1 The package's container] -.-> M1
  T7[T7 chromosomix, metacooler] -.-> M2
```

## Versions: 0.x, for good

plexembler stays at 0.x, the complete workflow included: there is no 1.0.0 in
this plan. CONVENTIONS.md § 18 allows a minor release before 1.0 to break the
interface, and that is the freedom kept.

- **A milestone is a minor release** (v0.1.0 for M1, … v0.6.0 for M6), and so
  is every feature after them: the minor number has no ceiling (v0.10.0 follows
  v0.9.0). `uv run bump-my-version bump minor` or `bump patch`, never `bump
  major`.
- **The complete workflow** is the first release in which every box of M1 to M6
  is ticked and `test_full` is green: a 0.x release like the others, named so
  in its changelog section and in the README.
- **A patch release only fixes**, and never breaks a run.
- **A minor release may break**, said where it does: a *Breaking changes*
  section in its changelog, and `docs/upgrading.md`, a page on moving a run to
  it (CONVENTIONS.md § 9). What can be renamed is first deprecated, with a
  warning, for one minor release (§ 18): a parameter, a samplesheet column, an
  output path, a key of the report.
- **Runs pin their release** (`nextflow run js2264/plexembler -r v0.7.2`):
  [usage](usage.md) says so, and every report names the release it was made
  with.
- **nf-core releases a pipeline from 1.0.0**, so plexembler stays
  nf-core-compatible rather than becoming an nf-core pipeline (T6).

## Decisions

Taken on 2026-11-26, from the maintainer's answers; each is written again where
it applies.

| Topic | Decision |
|---|---|
| Nanopore data | Mixed or unknown chemistries: a samplesheet column `chemistry` (`r10`, `r9`, or `NA`, told from the reads' median quality), R10.4.1 assumed until then |
| Strains | Haploid or homozygous: one collapsed assembly per species. hifiasm, purge_dups and HapHiC are opt-in |
| Hi-C | Arima's kit (DpnII and HinfI): its four junctions are the default; YaHS and instaGRAAL, in its `cpu` mode, are the scaffolders |
| Compute | Workstations with and without a GPU, SLURM clusters with and without GPU nodes: every GPU step is opt-in (`-profile gpu`) and has a CPU path |
| GPU steps | Not tested yet: `dorado polish`, instaGRAAL's MCMC and Helixer are marked experimental, in the docs and in their parameters' help, until their tests (nf-test's tag `gpu`) have been run on a GPU workstation. CI tests the CPU paths |
| The species of a pool | Of the same genus at most, some pairs at 85 to 95 % ANI: a read goes to a species only when its alignment is unique (`--deconvolve_mapq`) and near-identical (`--deconvolve_min_identity`), and the ambiguous share of each pair of species is reported (M2). Above 95 %, the pair is warned about (M1) |
| RNA-seq | Counts and TPM per species and library; no differential expression (M5) |
| Delivery | A pull request per milestone into `main`, squash-merged by the maintainer: each milestone on its own branch, its pull request stacked on the previous milestone's until that one is merged, then rebased onto `main` |
| Benchmark | Public data of yeasts with T2T references, mixed in silico (T2) |
| Benchmark sets | Three, from the easiest to the hardest (2026-11-27): each plants more of the repeats yeast genomes hold and gives worse long reads and Hi-C, R10 in all three. Milestone 4 counts the characterized genes, not the uncharacterized ORFs (T2) |
| Upstream | hicstuff, instaGRAAL, chromosomix and metacooler are attached: their changes are prepared there (T4, T7), and worked around here until they are released |
| The package's container | An image on GHCR per release, `ghcr.io/js2264/plexembler:X.Y.Z`, built from the released wheel; PyPI later (T1) |

## How a step is built

Every box below is done when all of this is:

1. **Its tools.** An nf-core module (`nf-core modules install`), or a local
   module in nf-core's form: `modules/local/<tool>/`, a conda environment and
   its Seqera container, a `stub:` block, versions in the `versions` topic, a
   `meta.yml`. Logic of plexembler's own is a `plexembler <command>` (Python,
   tested in-process, 80 % covered) that a local module runs: never a `bin/`
   script.
2. **Its subworkflow**, under `subworkflows/local/<name>/`, with nf-test: a real
   test and a `-stub` test, on `tests/data/` (the plumbing) and, from M1, on the
   benchmark pool's small set (T2).
3. **Its configuration.** Arguments and publishing in `conf/modules.config` (a
   species' results under `species/<species>/`, a pool's under `pools/<pool>/`),
   resources by label in `conf/base.config`, parameters in `nextflow.config` and
   `nextflow_schema.json` (then `docs/parameters.md`).
4. **Its documentation.** Its outputs in `docs/output.md`, its parameters and
   columns in `docs/usage.md`, its tools in `CITATIONS.md` and
   `toolCitationText()`, its entry in `docs/changelog.md`.
5. **Its numbers in MultiQC**, through MultiQC's own module for the tool or the
   step's `_mqc` file, and, from M6, in the species' report.

## M0 · Foundations

The repository, and what every later step relies on. Not released.

- [x] The pipeline on nf-core's template (nf-core/tools 4.2.2, strict syntax,
      topic channels, nf-schema 3, Nextflow ≥ 26.04.0), the package on
      chromokit's copier template, the docs site in chromokit's look.
- [x] The species and pooled libraries samplesheets: relative paths read from
      the samplesheet's folder, rows of a species or a library merged, pools
      checked against each other.
- [x] The reads' QC: NanoPlot per species, FastQC per pooled library, MultiQC.
- [x] `plexembler rename` (`<species>_<contig>`), and the command line's log
      (CONVENTIONS.md § 11).
- [x] CI: the package (`ci.yml`), the pipeline (`pipeline.yml`: lint, nf-core
      lint, nf-test on Nextflow 26.04.0 and the newest), the docs, the
      changelog check, the release, Dependabot.
- [x] Unknown values (`NA`) in the species samplesheet, and its `chemistry`
      column.
- [ ] GitHub's settings (T5).
- [x] The benchmark pool (T2): its three sets, `benchmark-data-v4`, released by `benchmarks.yml` (v1 never was; v2 was, by the workflow's first run).
- [x] The package's container, which M1's local modules need (T1).

## M1 · Contigs, per species · v0.1.0

**Goal:** for each species, its nanopore contigs, polished, with their checks;
usable as a draft genome on its own.

- [x] **Nanopore input.** FASTQ or Dorado's unaligned BAM (`samtools fastq -T '*'`,
      keeping its tags), a species' runs as one FASTQ; filtered by length and
      quality (`chopper`, the quality by chemistry); capped at
      `--ont_max_depth` (100×) of the longest and best reads (`filtlong
      --target_bases`). The samplesheet's `ont` takes `.bam`. The species'
      `chemistry` sets Flye's mode (`--nano-hq` for `r10`, `--nano-raw` for
      `r9`) and Medaka's model; `NA` is told from the reads' median quality
      (NanoPlot; Q15 and above: `r10`).
- [x] **What the samplesheet leaves unknown** (`NA`): the genome size, from the
      nanopore reads (LRGE), before the depth cap needs it; the k-mer spectrum
      of the species' nanopore reads fitted by GenomeScope2 as a diploid's,
      reported (nanopore errors inflate its heterozygosity: a hint, not an
      answer), the species treated as haploid. Its own WGS reads, once M2 has
      split them, are fitted again, with Smudgeplot. Each estimate is reported
      as one, beside what the samplesheet gave when it gave something.
- [x] **Assemblers, as candidates** (`--assemblers flye`, by default): Flye,
      `--no-alt-contigs` for a haploid species; hifiasm `--ont` opt-in
      (`--assemblers flye,hifiasm`), for R10 reads, `-l0` for a haploid species
      (nf-core: `flye`, `hifiasm`, `gfatools/gfa2fa`).
- [x] **Contigs named after their species:** `plexembler rename`, in a local
      module `PLEXEMBLER_RENAME` (plain FASTA: Medaka indexes it).
- [x] **Nanopore polishing** (`--ont_polisher medaka|dorado|none`): Medaka
      (nf-core `medaka`); `dorado polish` (local, `process_gpu`, experimental
      until tested on a GPU) for Dorado's BAMs.
- [x] **The genome checks**, a subworkflow `GENOME_QC` that every later stage
      reuses: `gfastats` (contiguity); BUSCO in genome mode, with the
      samplesheet's lineage or `--busco_lineage` (BUSCO's choice by default),
      offline with `--busco_db`; Merqury's QV against the pool's WGS k-mers
      (valid per species: an assembly k-mer missing from the pool is an error)
      or, for a pool without WGS, the species' own nanopore k-mers; its
      completeness against the species' own nanopore k-mers (the pool's short
      reads mix the species) (subworkflow `KMERS`); the assembly graph (GFA).
      compleasm's summary is read by `plexembler select`, for a later option.
- [x] **One candidate kept:** `plexembler select --step assembly`: BUSCO
      complete minus duplicated, then QV, then contig N50, each with a
      tolerance; a samplesheet column `assembler` names one instead. Every
      candidate stays under `species/<species>/candidates/assembly/`.
- [x] **Identity between the species of a pool:** `skani triangle`, then
      `plexembler ani` → a matrix and each pair's verdict; a warning above 95 %
      ANI, where the pool's short reads stop telling two species apart, and a
      pair of the same genus named as one from 85 %.

**Done when**, on the benchmark pool (T2): each species' BUSCO is within 1 point
of its reference's; QV ≥ 40 after polishing; QUAST against the T2T references
finds no misassembly over 10 kb in the kept candidate.

## M2 · The pool, split · v0.2.0

**Goal:** the pooled libraries split between the species; the contigs polished
with their own species' short reads (WGS, or Hi-C cut at its junctions), and
screened.

- [x] **`POOL_DECONVOLVE`**, the pool-level subworkflow
      ([concepts](concepts.md#deconvolution)): the contigs of every species of a
      pool in one reference; `bwa-mem2` (`-5SP` for Hi-C); pairs sorted by name;
      `plexembler split` (pysam, in a local module) gives each pair to the species
      of its best unique alignments: WGS proper pairs, Hi-C pairs whose mates
      both land in one species (pairs between species counted, not kept), at
      `--deconvolve_mapq` (20) and an identity of at least
      `--deconvolve_min_identity` (98 %): a read of a region missing from its
      own species' contigs would otherwise go, unique, to a sister species'
      copy, 85 to 95 % identical (Decisions). It writes each species' BAM (and FASTQ where a
      tool wants reads) and the pool's numbers: each species' share of each
      library, the ambiguous fraction of each pair of species, the pairs between
      species and the random-ligation rate ε (Hi-C: the pairs between species
      over 1 − Σpᵢ²), and, for Hi-C, the unique pairs between contigs (the
      screen's Hi-C signal); 10 % of each library's pairs held out, by a hash of
      their names, for the QV guard. Duplicates are dropped by `plexembler split`
      itself (the mates' 5' ends), not by `pairtools dedup`. To come: a pool's
      channels grouped with `groupKey(pool, n_species)`, so that a pool starts as
      soon as its own species are ready (today, every pool waits for all).
- [x] **The polishing reads** (`--polish_reads auto|wgs|hic`, `auto` by
      default): a species' WGS pairs when its pool has a WGS library; when it has
      none, or for a species whose WGS depth is under `--min_polish_depth`
      (30×), its **Hi-C reads, cut at their ligation junctions**. A Hi-C read is
      chimeric at a junction, but the pieces between junctions are genome, read
      with the same chemistry as WGS. `plexembler digest` cuts each read with
      hicstuff's rule (`cutsite_read`: 4 bases into the junction, so that each
      piece keeps its half of the site), from the junctions of the library's
      `enzymes` (DpnII's `GATCGATC`; Arima's four); drops the pieces under
      `--min_piece` (40 bp); and writes every piece once, as single-end reads:
      a Hi-C pair is two loci, not a fragment, and its pairing would mislead a
      polisher. Before cutting, duplicate pairs go (Hi-C libraries carry more
      PCR duplicates than WGS: `pairtools dedup` on the pool's Hi-C alignments),
      and a piece takes the species its pair was given. The reads with 0, 1 and
      2 or more junctions are counted: a ligation-efficiency figure for the
      report. (hicstuff's own `cutsite` pairs the pieces again; a mode of it
      that writes them once, unpaired, would replace `plexembler digest`, T4.)
      To come: `both`, and `--min_polish_depth`, Hi-C for a species whose WGS
      depth is too low.
- [x] **Short-read polishing** (`--sr_polisher polypolish,pypolca|none`):
      Polypolish on `bwa mem -a` alignments of the species' reads, single-end for
      Hi-C pieces (its insert-size filter needs true pairs), alignments with
      more than 5 clipped bases dropped (a junction the cutting missed); then
      pypolca. Hi-C's coverage is uneven, highest near restriction sites: a base
      it does not cover is left as it was, not changed. The guard: Merqury's QV
      before and after on k-mers of 10 % of the pairs, held out of the polishing
      (the polishing reads' own k-mers would flatter it, and Hi-C's uneven
      coverage makes a QV from its k-mers pessimistic); the polished version is
      kept only if it is better (`plexembler select --step polish`: the draft,
      Polypolish's and pypolca's QVs). pypolca runs `--careful`.
- [ ] **A polishing guard that sees past the held-out pairs' sampling.** On the
      benchmark pool the held-out pairs are 3 to 4× of each genome, at which 3
      to 7 % of a genome's 21-mers are missing from them by chance: the QV reads
      about 24 to 27 whatever the genome's true quality, and a polishing's gain or
      harm moves it by hundredths. There, measured against the references, the
      WGS pairs' polishing fixed 2 of Scer's 3 indels and 2 of Spar's 9 errors,
      and the Hi-C reads' added errors (Spar: 36 mismatches and 10 indels;
      Scer: 20 mismatches and 2 indels); the QV saw neither, and the guard's
      tolerance (1 QV) kept every draft: the small gains lost, the harm
      avoided. The second set (WGS at 100×, 10× held out) showed the gain lost
      too: Scer's polishing took its errors against the reference from 453 to
      267, and its held-out QV from 29.66 to 30.63 alone, under the tolerance.
      The polishing's tolerance is now 0.1, under which the first set's Hi-C
      polishing (0.00) stays out. To try still: the assembly-only k-mers the
      candidates differ in, rather than their QVs; a larger hold-out when the
      coverage allows it.
- [ ] **Micro-C** (`enzymes` is `none`) has no junction sequence to cut at:
      such a pool's species stay polished with nanopore reads alone, by default.
      Polishing with the clipped local alignments of its reads (primary and
      supplementary segments, the clipped bases left out) is tried on the
      benchmark first, and offered only if it helps.
- [ ] **The pool's Hi-C map, species as genomes:** every Hi-C pair of the
      pool, those between species included, binned on the pool's contigs as a
      `.metacool` whose genomes are the species (`metacooler create pairs`, the
      contig-to-species table from the contigs' names). `metacooler qc`, `ps
      --by genome` and `network` give each species' contacts, P(s) and the
      contacts between species, the noise floor's raw material; chromosomix
      opens it as it opens a metagenome's map.
- [x] **The contig screen,** `plexembler screen`, which weighs three signals for
      each contig: the nanopore reads of every species of the pool mapped on the
      pool's contigs (`minimap2 -x map-ont`, the pool's cross-talk matrix); its
      identity to the other species' contigs (`minimap2 -x asm10`); and its Hi-C
      contacts outside its species, from `plexembler split`'s unique pairs
      between contigs (`chromosomix magpurity` on the pool's `.metacool`, once
      the pool has one). A verdict per contig, keep, move, drop or organelle,
      with its evidence and thresholds; a move renames the contig
      `<species>_<contig>` (`metacooler reassign` for the pool's map, then).
- [x] **Organelles and plasmids:** the mitochondrial contig, by circularity, GC,
      depth and `miniprot` hits of a small set of yeast mitochondrial proteins
      (`assets/`), kept apart as `<species>_chrM`. Today by depth and GC (the
      screen's `organelle` verdict); circularity (Flye's graph) and miniprot's
      hits to come.
- [ ] **Foreign DNA:** FCS-adaptor always; FCS-GX with `--fcsgx_db` (nf-core's
      subworkflow `fasta_clean_fcs`); Tiara's labels (eukaryote, prokaryote,
      organelle) for the report.
- [ ] **Haplotigs:** purge_dups for a species whose new column `purge` says so;
      never by default, never for a hybrid (it would remove a subgenome).
- [ ] **The pool's Hi-C on the screened contigs** (a second `POOL_DECONVOLVE`),
      which M3 scaffolds with.

**Done when**, on the benchmark pool: ≥ 99.5 % of the assigned pairs go to their
true species (the reads' names carry it) for pairs of species below 90 % ANI;
ε within 10 % of the simulated rate; every planted bleed-through contig moved
or dropped; no species' QV lower after polishing; polishing with digested Hi-C
alone gains at least 80 % of what WGS gains, species by species.

## M3 · Chromosomes · v0.3.0

**Goal:** each species' genome, scaffolded with its Hi-C pairs and finished, and
its Hi-C maps.

- [x] **Pairs per species:** `pairtools parse/sort/dedup` on the species'
      Hi-C alignments to its genome (4DN `.pairs.gz`, with the `#chromsize`
      lines instaGRAAL reads; pairtools' stats in MultiQC). On the contigs, for
      instaGRAAL, to come with it.
- [ ] **Scaffolders, as candidates** (`--scaffolders yahs,instagraal,haphic`):
      YaHS done, the only one for now (its candidate in
      `candidates/scaffold/yahs/`); instaGRAAL waits for its CPU image (T4):
      - YaHS (nf-core `yahs`), its motifs from the library's `enzymes` (none for
        Micro-C), in two passes, the second (`-a`) on the first's scaffolds;
      - instaGRAAL, a local module `INSTAGRAAL_ENDTOEND`:
        `instagraal-endtoend contigs.fa pairs.gz -e <enzymes> --mode cpu`, or
        `--instagraal_mode default` (MCMC and refine) on a GPU (`process_gpu`,
        experimental until tested on a GPU);
        not for a library made without an enzyme;
      - HapHiC (local), opt-in, for species with a `ploidy` of 2 or more and
        an `n_chromosomes`.
- [ ] **Layouts as AGP:** `plexembler agp`: instaGRAAL's `final_info_frags.txt`
      to AGP (until instaGRAAL writes one, T4); AGP to FASTA, checked against the
      contigs.
- [x] **Each join scored:** `plexembler junctions`: the pairs across each
      join against the species' own contact decay P(s), from the pairs within
      its contigs, and the joins they do not support (`genome/qc/junctions/`).
      `instagraal-refine` as a second opinion, with instaGRAAL.
- [ ] **Telomeres:** `tidk explore` finds the repeat when the samplesheet has
      none, `tidk search` counts the chromosomes with telomeres at both ends
      (nf-core's subworkflow `fasta_explore_search_plot_tidk`).
- [ ] **One candidate kept:** `plexembler select --step scaffold`: no unsupported
      join (unless every candidate has some), then the most T2T chromosomes,
      then the number of scaffolds nearest `n_chromosomes` (a term left out when
      it is `NA`, and HapHiC, which needs it, skipped), then N50; a samplesheet
      column `scaffolder` names one instead. The karyotype found, T2T
      chromosomes first, is reported against the one given.
- [ ] **Curation:** a samplesheet column `curated_agp` replaces the kept layout
      with a curated one, checked against the contigs; each candidate's `.pretext`
      and `.hic` with its AGP are written for PretextView and Juicebox.
- [x] **Finishing:** `plexembler finalize`: chromosomes named by size
      (`<species>_chr01`, …, the samplesheet's `n_chromosomes` of them, or every
      scaffold of `--min_chromosome` or more), `<species>_chrM`,
      `<species>_unplaced_01`; bgzip, `.fai`, `.gzi`, `chrom.sizes`, the names'
      table, the AGP in the new names; a species without Hi-C finished from its
      contigs. To come: gap closing (TGS-GapCloser, local) and telomere
      extension (teloclip, local) behind `--close_gaps` and
      `--extend_telomeres`; a related reference (`--orient_reference`) to
      orient the chromosomes (`ragtag`), never reorder them.
- [ ] **The pool on the final genomes** (a third `POOL_DECONVOLVE`: WGS, Hi-C).
      The second, on the screened contigs, is done: the pairs each species is
      scaffolded and mapped with.
- [ ] **Maps:** `cooler cload pairs`, `zoomify --balance` (`.mcool`, done:
      `species/<species>/hic/`); to come, `.hic`
      (hictk, local); `.pretext` and PNG (`pretextmap`, `pretextsnapshot`); the
      pool's map on the final genomes, a `.metacool` again, where the blocks
      between species show the noise.
- [ ] **The maps' checks, by chromosomix's command line:** `chromosomix
      pairsqc` per species (pair classes against the library's enzymes, distance
      to the cut, the four orientations' P(s)), `chromosomix stats` (pairtools'
      counts, for MultiQC), `chromosomix mapreport` (the map's report, its
      MultiQC sections), `chromosomix expected` (the P(s) that
      `plexembler junctions` scores joins against).

**Done when**, on the benchmark pool: ≥ 90 % of the reference chromosomes are
T2T in the kept candidate; no misjoin against the T2T references; the
benchmark's planted misjoins are flagged by `plexembler junctions`.

## M4 · Annotation · v0.4.0

**Goal:** each species' genes (GTF, GFF3, proteins), repeats, non-coding RNAs and
mitochondrial genes.

- [x] **Repeats:** RepeatModeler2 (a library per species), RepeatMasker →
      soft-masked genome (`--repeats`). The repeat landscape, to come.
- [x] **RNA-seq, split:** STAR on the pool's final genomes (introns capped at
      5 kb, `--alignIntronMax 5000 --alignMatesGapMax 5000`: compact genomes get
      false splices otherwise), `plexembler split --type rnaseq` (unique hits,
      each species' pairs kept as alignments, already on its genome: no second
      alignment); each species' libraries merged for the gene models.
- [x] **Genes:** BRAKER3 `--fungus` (nf-core `braker3`) with the species' RNA
      alignments and OrthoDB's fungal proteins (`--orthodb_proteins`, given;
      its download to the cache, to come), as hints to AUGUSTUS's model of
      `--augustus_species`, without training: on the benchmark's second set,
      GeneMark-ET failed on S. cerevisiae's 14 introns, and the models BRAKER3
      trained on Spar and Cgla found 88.8 and 87.4 % of their genes, where
      S. cerevisiae's model found 94.0 and 95.1 % ab initio. `--braker_train`
      trains on the species (far from any model, rich in introns). Without
      RNA-seq, BRAKER on the proteins alone; with neither, AUGUSTUS ab initio
      with a trained species (`--augustus_species`, local module). To come: GALBA, if it does better
      than BRAKER on proteins on the benchmark; Helixer (local, GPU,
      experimental until tested on one) as a second candidate
      (`--gene_predictors braker,helixer`), scored by BUSCO's proteins and
      OMArk.
- [x] **Formats:** `plexembler genes`, in place of AGAT: GTF and GFF3, locus
      tags `<species>_000010`, proteins translated with the species'
      `genetic_code` (CUG is a sense codon in tables 1, 12 and 26 alike: gene
      structures do not change, proteins do; a protein with a stop inside is
      warned of), CDS, a table of the genes. AGAT's statistics, to come with
      the checks.
- [ ] **The genetic code, checked:** Codetta (Shulgina and Eddy 2021) infers each
      species' nuclear code from its genome: it is the code used when the
      samplesheet says `NA`, and a warning when it differs from the one given.
- [ ] **Non-coding RNAs:** tRNAscan-SE (`-E`), barrnap (`--kingdom euk`: the rDNA
      locus), Infernal with Rfam. tRNAscan-SE and barrnap done; Infernal, to
      come (`--rfam_cm`).
- [ ] **Functions:** eggNOG-mapper (orthologue names, GO, KEGG, EC, Pfam) →
      `functions.tsv`; InterProScan with `--interproscan_db`.
- [ ] **The mitochondrial genome's genes,** with the yeast mitochondrial code
      (table 3): MitoHiFi in its contig mode, or MFannot; the one that annotates
      the benchmark's mitochondria best.
- [ ] **Annotation checks:** BUSCO's proteins, OMArk (local), AGAT's statistics,
      the RNA-seq reads in genes. BUSCO's proteins done; the share of the
      RNA-seq reads in genes comes with M5's counts.

**Done when**, on the benchmark species with reference annotations: ≥ 90 % of
the characterized reference genes found, the uncharacterized ORFs left out
(Decisions; a predicted CDS over 80 % of each reference CDS, at 95 % identity:
the genome's coordinates are not the reference's, which gffcompare would need);
BUSCO's proteins within 2 points of BUSCO's genome.

## M5 · Expression and tracks · v0.5.0

**Goal:** each species' expression tables, counts and TPM, and its coverage
tracks. Differential expression is out of scope (Decisions).

- [x] **Counts:** featureCounts on each RNA-seq library's per-species
      alignments, its strandedness `--rnaseq_strandedness` (from Salmon's
      library type, to come).
- [ ] **TPM:** Salmon on the pool's transcripts, every species with its genome
      as decoy, so that a read of a gene shared by two species is shared by
      Salmon's EM; TPM within each species, from that species' reads alone.
      `plexembler expression` writes `<species>.counts.tsv` and
      `<species>.tpm.tsv`: genes × libraries (`<condition>_<replicate>`), and each
      gene's ambiguous fraction. Done: `plexembler expression`'s counts and
      TPM, from featureCounts' counts and the genes' lengths, the split's
      pairs alone; to come: Salmon's EM for the shared reads, and each gene's
      ambiguous fraction.
- [ ] **Tracks:** bigWig of the RNA-seq (by strand), the WGS and the nanopore
      reads (`deeptools bamCoverage`, `mosdepth`). The RNA-seq's and the
      nanopore reads' done; the WGS', to come.

**Done when**, on the benchmark pool's simulated RNA-seq (known TPMs):
Spearman's ρ ≥ 0.9 for every species.

## M6 · Reports and extras · v0.6.0

- [x] **`plexembler report`**, on chromokit's report kit (CONVENTIONS.md § 12):
      per species, `report.txt`, its MultiQC row and its page: the summary;
      its share of the pool; the candidates and the choice; contamination; a
      karyotype (chromosomes, telomeres, centromeres, rDNA, chrM); BUSCO, QV,
      OMArk; the annotation; the expression; versions and methods. A pool's
      report beside it (`pools/<pool>/`). Done on the kit's own page: the
      assembly and polishing, the share of each library, the screen, the
      chromosomes, BUSCO and QV, the joins, the genes and BUSCO on the
      proteins, the expression. To come: its own template and plots, the
      telomeres, centromeres and rDNA in the karyotype, the pool's report.
- [ ] **Contamination, further:** EukCC; a GC × depth plot of the contigs,
      labelled by Tiara, FCS-GX and the screen.
- [ ] **Ploidy, heterozygosity, aneuploidy:** Clair3 (nanopore) and bcftools (WGS)
      → allele balance; `mosdepth` by chromosome → copy number; Smudgeplot on the
      species' WGS reads.
- [ ] **Centromeres from Hi-C:** `plexembler centromeres`, from the
      inter-chromosomal contacts of the Rabl configuration, as centurion does
      (Varoquaux et al. 2015; Marie-Nelly et al. 2014) → `centromeres.bed`.
- [ ] **Each species' identity:** ITS (ITSx, local) and the D1/D2 domain of the
      28S rRNA, against a local UNITE or NCBI set (BLAST): sample swaps show.
- [ ] **Comparative**, per pool (`--comparative`): OrthoFinder (orthogroups,
      species tree), synteny (GENESPACE, local), the ANI matrix of M1.
- [x] **A browser for the delivery: chromosomix.** Done: the store, checked
      by chromosomix's own `check_store` in the package's tests; the
      species' workspaces, which name absolute URLs only, wait for
      `--browser_base`. `plexembler browser` writes,
      for each species, a workspace (`<species>.chromosomix.json`: its map with
      its genome, genes, repeats, telomeres, centromeres and coverage tracks
      placed, a colour scale, a first window), and, at the top of the results, a
      store (`chromosomix.store.json`, `layout: cards`, one card per species,
      grouped by pool) whose rows name every map, genome, track and workspace
      by a relative path. The results then open in chromosomix as they are: its
      web app, desktop app or single-file HTML with the results folder
      connected, or chromosomix.app with `?store=<url>` once the folder is on a
      bucket; the species' report embeds its view (`?embed=1`) where the results
      are served. Checked with chromosomix's own `scripts/make_store.py check`.
- [ ] **The complete workflow** (a 0.x release, see [Versions](#versions-0x-for-good)):
      `test_full` green on the benchmark, every page written, a metro map
      (`manifest.diagram`), the benchmarks page (CONVENTIONS.md § 20).

## Tracks

### T1 · The package's container

The local modules that run `plexembler <command>` need it in a container: an
image on GHCR per release (Decisions), a deviation from CONVENTIONS.md § 15,
which publishes images only for what pip cannot install. Its reason: the
pipeline's own steps need a container before the package is anywhere else.

- [x] `Dockerfile` and `.dockerignore` (the released wheel, on a pinned
      micromamba base; chromokit's template with `image: true`), and
      `release.yml`'s image job: `ghcr.io/js2264/plexembler:X.Y.Z` and
      `latest`, linux/amd64 and linux/arm64, `docker run … plexembler --version`
      before it is pushed.
- [x] PyPI's publishing, opt-in: `release.yml` publishes to PyPI only once the
      repository variable `PUBLISH_PYPI` is `true` (after PyPI's trusted
      publisher is set up), so that a release does not fail before.
- [x] `pipeline.yml` builds `ghcr.io/js2264/plexembler:<version>` from the
      checkout before nf-test, so that a pull request tests its own code.
- [x] The local modules that run `plexembler`: `container
      "ghcr.io/js2264/plexembler:${workflow.manifest.version}"`, and for conda
      an environment that installs the package at the release's tag (M1).

### T2 · Test data and the benchmark pool

- [x] `tests/data/`: two species of a few dozen kilobases and their pool, for
      the plumbing (`make_test_data.py`).
- [x] **`benchmark-data-v1`**, a release asset (CONVENTIONS.md § 13, § 20),
      and its generator beside it (`benchmarks/make_benchmark_data.py`): two or
      three chromosomes and the mitochondrion each of three yeasts with complete
      references, two of a genus (*S. cerevisiae*, *S. paradoxus*) and one far
      (*N. glabratus*); nanopore reads by Badread (R10 error model) with
      planted bleed-through; a pooled WGS by InSilicoSeq; a pooled Arima Hi-C
      simulated by the generator itself (ligation junctions in the reads,
      contacts decaying with distance, noise between species, PCR duplicates).
      `-profile benchmark`; `benchmarks/evaluate.py` checks a milestone's "Done
      when"; `.github/workflows/benchmarks.yml` fetches the release, or makes
      and publishes it, and runs both.
- [x] **`benchmark-data-v2`**: a 9 kb repeat planted every 200 kb and
      nanopore reads under 7 kb, so that the drafts break where the scaffolds
      must join them; errors the nanopore reads share, for the polishing;
      the references' annotations (M4) and RNA-seq at known levels (M5);
      misjoins laid out at evaluation (M3). What it found (js2264/plexembler#8
      holds the numbers): the polishing guard too strict for the held-out
      QV's floor (its tolerance now 0.1); YaHS's own resolutions joining
      nothing below a contig N50 of 100 kb, and a clear join left apart (its
      resolutions given, from 1 kb, and a second pass on its scaffolds);
      BRAKER3's training failing on a yeast's few introns (the clade's
      AUGUSTUS model, `--braker_train` to train). And what it could not
      measure: reads under 7 kb collapse the yeasts' own subtelomeric
      paralogs as well as the planted repeat, so chromosome ends went
      missing before any scaffolding.
- [x] **`benchmark-data-v3`**: three sets, from the easiest to the hardest,
      each planting more of the repeats yeast genomes hold, S288C's own (Ty1
      upstream of tRNA genes, Ty1 in tandem, Y' arrays inside the telomeres,
      an rDNA array), with R10 reads and Hi-C less able to resolve them
      (`benchmarks/README.md`); each reference gene characterized or not,
      milestone 4 counting the characterized ones; milestone 3 leaving out the
      tandem arrays no read crosses. `benchmarks.yml` runs the three, the easy
      set gating, the medium and hard ones reported.
- [x] **`benchmark-data-v4`**: the easy set's nanopore reads without
      Badread's artifacts (junk, random and chimeric reads, glitches), its
      recipe for good reads. Its junk reads, a short motif repeated over the
      whole read (up to 37 kb), stalled Flye for hours on `Scer`: two of them
      share millions of k-mer matches, which its overlapper chains in
      quadratic time. The
      medium and hard sets keep Badread's defaults, their reads v3's (the
      seed kept).
- [ ] **Long junk nanopore reads** (js2264/plexembler#8): a short motif
      repeated over tens of kb, which real runs have too, stalls Flye as
      above. Drop such reads before assembling (a read mostly one
      short-period repeat), then give the easy set Badread's junk reads back.
- [ ] **What the medium set found** (js2264/plexembler#8), the next pull
      requests; the medium set gates the benchmarks once they pass it:
      - YaHS joins two chromosomes' ends through a subtelomeric Y' the
        assembly collapsed (Scer, in its first pass at 5 kb; Spar on the hard
        set): `plexembler junctions` flags the join (6 % of the pairs
        expected across it), and nothing cuts it. Cut the joins it flags, or
        keep the scaffolder candidate without them (M3, "One candidate kept").
      - The Hi-C polishing is never kept: the held-out pairs of a digested
        Hi-C library are too thin for the guard's QV to see a gain (M2).
      - The first split's epsilon, on the contigs before the screen, is 14 %
        low (0.0259 of 0.03); the second split's is within 10 % (M2).
      - A circular contig that Flye writes as two copies of itself (Cgla's
        mitochondrion, once, on the first set) is kept as it is.
- [ ] **`test_full`**: an in-silico pool of public data (nanopore, Hi-C, Illumina,
      RNA-seq) of five to ten yeasts with T2T references, subsampled and mixed,
      on Zenodo; run by hand on a cluster before each release, its numbers on
      the benchmarks page.

### T3 · Nextflow's language

- [x] Strict syntax (the default since Nextflow 26.04), topic channels for the
      versions, nf-schema 3.
- [ ] **Workflow outputs** (`output {}` and `publish:`) in place of `publishDir`,
      once nf-core's template has them: their index becomes the delivery's
      manifest (`<outdir>/manifest.csv`).
- [ ] **Static types and records:** nf-core adds them to its template in Q4 2026
      and makes them mandatory in Q2 2027. Local modules take one meta map and
      give one back, so that their rewrite as records is mechanical;
      `nf-core pipelines sync` brings the template's changes, `uvx copier update`
      chromokit's.

### T4 · instaGRAAL and hicstuff, upstream

- [ ] An AGP beside `final_genome.fa` (then `plexembler agp` loses that job).
- [ ] A seed for the MCMC, for reproducible runs and tests.
- [ ] Versioned image tags, a CPU-only image for `--mode cpu`, and a bioconda
      recipe (then Seqera Containers and BioContainers for free).
- [ ] nf-core modules for `instagraal-endtoend` and `hicstuff pipeline`.
- [ ] hicstuff: a `cutsite` mode that writes each piece once, unpaired, for
      polishing (M2), so that `plexembler digest` can go.

### T7 · chromosomix and metacooler

The project's own tools for maps (CONVENTIONS.md § 22), which plexembler uses
from M2 on.

- [ ] **metacooler** in a container: its bioconda recipe (metacooler's item 1.1),
      or the Seqera Containers image of its PyPI release until then.
- [ ] **chromosomix's command line** in a container: its npm package needs Node
      22.12 or newer, which neither bioconda nor Seqera Containers' pip path
      gives. An image published by chromosomix's release
      (`ghcr.io/js2264/chromosomix:<version>`, CONVENTIONS.md § 15) is the clean
      way; a local module that installs a pinned version from npm at run time
      is not, since a cluster's nodes may have no network.
- [ ] **The workspace and store formats** pinned to a chromosomix version, and
      checked in plexembler's tests with that version's `make_store.py check`.

### T5 · GitHub and the infrastructure

- [ ] `main` as the default branch (the repository was empty: the first branch
      pushed may have become the default).
- [ ] The default branch's ruleset (CONVENTIONS.md § 16): pull requests,
      squash-merging, and the required checks: CI's `lint` and stable `test`
      legs, the docs build, the pipeline's lint and nf-test.
- [ ] Pages from GitHub Actions; a `CODECOV_TOKEN` secret (koszullab's
      organisation secret does not reach js2264's repositories); Zenodo's
      integration; CodeQL and the dependency review (CONVENTIONS.md § 19).

### T6 · nf-core

- [ ] Each local module that is useful beyond plexembler, contributed to
      nf-core/modules once stable (instaGRAAL, hicstuff, metacooler, chromosomix,
      HapHiC, hictk, `dorado polish`, OMArk, Helixer).
- [ ] plexembler itself stays here, nf-core-compatible: nf-core releases a
      pipeline from 1.0.0, which plexembler does not plan
      ([Versions](#versions-0x-for-good)), and nf-core's layout would replace some
      of the conventions' choices (`CHANGELOG.md`, the docs' site).

## Questions that change defaults

All answered ([Decisions](#decisions)):

1. ~~**The nanopore chemistry and basecalling**~~: mixed or unknown, a
   `chemistry` column.
2. ~~**The strains**~~: haploid or homozygous.
3. ~~**The Hi-C protocol of each pool**~~: Arima.
4. ~~**Where it runs**~~: all four of workstation and cluster, with and without
   GPUs.
5. ~~**How close the species of a pool are**~~: of the same genus at most, some
   pairs at 85 to 95 % ANI.
6. ~~**The RNA-seq design**~~: counts and TPM; no differential expression.

## Risks

| Risk | What it would do | What the plan does about it |
|---|---|---|
| Two species of a pool too close | Their short reads cannot be told apart in conserved regions; polishing could mix them | The ANI warning (M1); ambiguous reads go to neither species (M2); the QV guard keeps a polishing only when it helps |
| A species scarce in a pooled library | Too few Hi-C pairs or RNA reads for it | Each species' share and depth reported (M2); its steps run with a warning below thresholds |
| Bleed-through between nanopore runs | Contigs of one species in another's assembly | The three-signal screen (M2) |
| The Rabl configuration | Scaffolders join chromosomes at their centromeres | Each join scored against the species' P(s) (M3); and the same signal annotates the centromeres (M6) |
| An MCMC run differs from the next | instaGRAAL's results change run to run | `cpu` mode by default; a seed upstream (T4) |
| Large databases (FCS-GX, InterProScan) | Steps that cannot run on a workstation | Optional, behind their `--*_db` parameter; skipped with a note in the report |
| nf-core's template and language move | Rewrites of the template's files | `nf-core pipelines sync`; local modules ready for records (T3) |
