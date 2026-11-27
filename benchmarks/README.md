# The benchmark pool

What each milestone's "Done when" ([docs/plan.md](../docs/plan.md)) is checked
on: a pool of yeasts whose truth is known, simulated from their references, in
three sets, from the easiest to the hardest.

## The species

`make_benchmark_data.py` writes each set from three references at NCBI and
their annotations:

| Species | Strain | Reference | Kept | Genes (characterized) |
|---|---|---|---|---|
| `Scer` | *Saccharomyces cerevisiae* S288C | GCF_000146045.2 | chrI, chrIII, chrVI, chrM | 391 (339) |
| `Spar` | *S. paradoxus* CBS432 | GCA_002079055.1 | chrI, chrIII, chrVI, chrM | 349 (287) |
| `Cgla` | *Nakaseomyces glabratus* CBS138 | GCF_000002545.3 | chrA, chrB, chrM | 412 (371) |

`Scer` and `Spar` are a pair of a genus, about 88 % identical: their pooled
short reads are the hard ones to split (docs/plan.md, Decisions). `Cgla` is far
from both. Genes are the nuclear protein-coding ones of the references'
annotations; a gene is characterized unless its annotation calls it
"uncharacterized protein" (`Scer`) or "hypothetical protein" (`Spar`), or, for
`Cgla`, whose annotation calls every gene uncharacterized, unless its note names
an *S. cerevisiae* homolog that is itself characterized
(`truth/<species>.genes.tsv`, each gene's class and why).

## The three sets

Each set plants more of the repeats yeast genomes hold, and gives nanopore
reads and Hi-C less able to resolve them. The repeats are S288C's own: its Ty1
element YARCTy1-1 (chrI, 5,925 bp from LTR to LTR), which lands upstream of a
tRNA gene, as Ty1 does, and duplicates the 5 bases it lands on; the Y' element
at chrXIV's left end (6,617 bp), planted just inside a chromosome's telomeric
repeat, tandem copies separated by a short telomeric repeat; and a unit of its
rDNA (chrXII, 9,137 bp). Each species' families are S288C's diverged (`Spar`'s
by 8 %, `Cgla`'s by 25 %), each copy diverged a little more from its family's
(Ty1 0.2 %, Y' 0.5 %, rDNA 0.05 %). No repeat lands in a gene.

| | `easy/` | `medium/` | `hard/` |
|---|---|---|---|
| Ty1 elements | one per chromosome | one every 100 kb | one every 70 kb, and three in tandem once per chromosome |
| Y' elements | none | at half the chromosomes' ends, one or two in tandem | at three quarters of the ends, two to four in tandem |
| rDNA | none | none | twelve units in tandem on each species' longest chromosome |
| Nanopore reads (R10's models) | 99.5 % identity, N50 about 25 kb, 50× | 98 %, about 11 kb, 40× | 96 %, about 6 kb, 30× |
| Their shared errors | one every 50 kb | one every 20 kb | one every 10 kb |
| Badread's junk, random and chimeric reads, its glitches | none (its recipe for good reads) | 1 % each, a glitch every 10 kb (its defaults) | the same |
| Hi-C pairs, random ligations | 1,000,000, 1 % | 400,000, 3 % | 150,000, 6 % |
| Hi-C between chromosomes, unligated, duplicates | 10 %, 3 %, 3 % | 12 %, 5 %, 5 % | 15 %, 8 %, 10 % |

| Reads | How | Truth in the read's name |
|---|---|---|
| Nanopore, per species (`ont/<species>.fastq.gz`) | Badread 0.4.1, R10.4.1's error and quality models (`nanopore2023`), at the set's identity, length, depth and artifacts; from the genome with its shared errors (70 % of them a homopolymer a base too long or too short, the others substitutions, outside the planted repeats), as a nanopore consensus keeps them, for the short reads' polishing to fix; bleed-through: 15 % of `Spar`'s reads are `Cgla`'s (enough to assemble some), 1 % of the others' another species' | the sequence, strand and position |
| WGS, pooled (`pool/wgs_R{1,2}.fastq.gz`) | InSilicoSeq 2.0.1, NovaSeq's errors, 100× of each species: its held-out pairs (10 %) see a genome's QV past 40 | the sequence |
| Hi-C, pooled (`pool/hic_R{1,2}.fastq.gz`) | here: Arima's DpnII and HinfI, the set's pairs; each pair two restriction fragments' ends ligated, a read through the junction when the molecule is short; contacts decaying as 1/s; random ligations a partner from anywhere in the pool, by abundance (ε) | both ends, their sites and enzymes |
| RNA-seq, pooled (`pool/rna-<condition>-<replicate>_R{1,2}.fastq.gz`) | here: two conditions (YPD, heat) of two replicates, 500,000 pairs of 100 bp each, reverse-stranded (dUTP); 40, 35 and 25 % of the pairs `Scer`'s, `Spar`'s and `Cgla`'s; each gene's fragments (250 bp on average) drawn by its level and its length, its level log-normal, 15 % of the genes changed in heat (a log2 fold from N(0, 1.5)), each replicate's noise 15 % | the species and the gene |

Each species' mitochondrion is 5 copies per nuclear genome in every library
but the RNA-seq. `truth.json` records every number above, the set's level and
where each repeat sits (its family, copies, and whether it is a tandem array);
`species.csv` and `libraries.csv` are the run's samplesheets. The truth:
`references/<species>.fa` (the genome, its repeats in),
`references/<species>.gtf` (its genes, in the genome's coordinates, each
protein-coding gene's `gene_class`), `references/<species>.cds.fa` (its nuclear
protein-coding genes' CDS), `truth/<species>.genes.tsv` (those genes, each
characterized or not), `truth/<species>.tpm.tsv` (each gene's true TPM in each
RNA-seq library) and `truth/<species>.ont_errors.tsv` (the nanopore reads'
shared errors).

The sets are deterministic: the same `GENERATOR_VERSION` writes the same reads,
and a version keeps the seed of the one before (`SEED`), so that changing one
set leaves the others' reads as they were.
They are published, with that version, as the release `benchmark-data-v<N>`
(`<set>.tar`, and `<set>.sha256`, each file's checksum), which
[the benchmarks' workflow](../.github/workflows/benchmarks.yml) downloads, or
makes and publishes when it does not exist.

```bash
# A set, without the simulators
gh release download benchmark-data-v4 --repo js2264/plexembler --dir benchmarks/data --pattern 'medium.*'
tar -xf benchmarks/data/medium.tar -C benchmarks/data
(cd benchmarks/data && sha256sum -c medium.sha256)
```

Or make them (`--level` for one set):

```bash
# The simulators: Badread and InSilicoSeq
micromamba create -p .bench -c conda-forge -c bioconda python=3.12 badread=0.4.1 insilicoseq=2.0.1
PATH=$PWD/.bench/bin:$PATH python benchmarks/make_benchmark_data.py --outdir benchmarks/data --threads 4
```

## A run, and its checks

```bash
nextflow run . -profile benchmark_medium,docker --outdir results
python benchmarks/evaluate.py --results results --data benchmarks/data/medium --milestone 1
```

The profiles `benchmark_easy`, `benchmark_medium` and `benchmark_hard` run the
pipeline on a set (`benchmark` is the medium one), each task within GitHub's
runner (4 CPUs, 12 GB). `evaluate.py` runs BUSCO on the references, QUAST and
minimap2 against them, and plexembler's own `junctions` (`--busco`,
`--quast`, `--minimap2`, `--plexembler`: how to run them, on PATH by default)
and writes `results/benchmark/m<N>.md`; it exits 1 when a check fails (a check
with nothing to judge reads "not judged"). Each milestone's checks are its
"Done when":

- milestone 1, the kept contigs: BUSCO against the reference's,
  misassemblies; and the QV after polishing, against the reference (QUAST's
  mismatches and indels);
- milestone 2, the two splits (on the contigs, then on the screened contigs)
  and their epsilon, the screened contigs and the polishing: its QV, and its
  errors against the reference, QUAST's (`--hic-results`: a second run,
  polished with the Hi-C reads, `--polish_reads hic`);
- milestone 3, each genome's chromosomes against the reference's (whole, end
  to end; the planted tandem arrays, which no read crosses, left out, an array
  at a chromosome's end moving that end inside it), its misassemblies and its
  flagged joins; and, on its own Hi-C pairs, its chromosomes cut in two and
  laid back (each join supported), and planted misjoins, an inversion and a
  translocation (each flagged);
- milestone 4, its characterized reference genes found among its predicted
  ones (a CDS covering 80 % of the reference's at 95 % identity), BUSCO on its
  proteins against BUSCO on its genome;
- milestone 5, Spearman's rho between its TPM and the true TPM, in every
  RNA-seq library, on its genes paired with the reference's.

## In CI

[`benchmarks.yml`](../.github/workflows/benchmarks.yml) runs the pipeline on
each set, a job each, in parallel (Docker, on GitHub's runner, the package's
image built from the checkout), then the checks: every set's and every
milestone's for every release tag (`v*`), and by hand (Actions, Benchmarks,
Run workflow: a branch, one set or all, one milestone or all, and whether to
publish a set the run had to make). The easy set must pass; the medium and
hard sets measure how far the pipeline gets, their failures reported, not the
run's, until the pipeline passes the medium set too (docs/plan.md, T2). Each
milestone's table lands in the run's summary, its report in the run's
artifact. A run takes hours, so a pull request has the pipeline's
nf-tests instead (`pipeline.yml`).

## To come

`test_full` is real data: public reads of five to ten yeasts with complete
references, mixed in silico (docs/plan.md, T2). The sets' next versions:
misjoins the scaffolder itself makes (planted ones are laid out where they are
scored), a species of R9.4.1 reads, for the chemistries' mix, and Badread's
junk reads in the easy set too, once the pipeline drops them before assembling
(docs/plan.md, T2).
