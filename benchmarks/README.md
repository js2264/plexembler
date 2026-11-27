# The benchmark pool

What each milestone's "Done when" ([docs/plan.md](../docs/plan.md)) is checked
on: a pool of yeasts whose truth is known, simulated from their references.

## The small set

`make_benchmark_data.py` writes it to `data/small/`, from three references at
NCBI and their annotations:

| Species | Strain | Reference | Kept | Genes | Size |
|---|---|---|---|---|---|
| `Scer` | *Saccharomyces cerevisiae* S288C | GCF_000146045.2 | chrI, chrIII, chrVI, chrM | 391 | 0.93 Mb |
| `Spar` | *S. paradoxus* CBS432 | GCA_002079055.1 | chrI, chrIII, chrVI, chrM | 349 | 0.99 Mb |
| `Cgla` | *Nakaseomyces glabratus* CBS138 | GCF_000002545.3 | chrA, chrB, chrM | 412 | 1.05 Mb |

`Scer` and `Spar` are a pair of a genus, about 88 % identical: their pooled
short reads are the hard ones to split (docs/plan.md, Decisions). `Cgla` is far
from both. Genes are the nuclear protein-coding ones of the references'
annotations.

Each species' genome is its reference with a repeat planted in it, so that
its assembly breaks and its scaffolding has joins to make: a 9 kb sequence
of its own, identical in every copy, a copy for every 200 kb of each nuclear
chromosome (one at least), between genes, 40 kb or more from the ends. No
nanopore read is longer than 7 kb, so none crosses a copy. And the nanopore
reads share errors, as a nanopore consensus keeps them, for the short reads'
polishing to fix: one every 20 kb, outside the repeats, 70 % of them a
homopolymer a base too long or too short, the others substitutions.

| Reads | How | Truth in the read's name |
|---|---|---|
| Nanopore, per species (`ont/<species>.fastq.gz`) | Badread 0.4.1, R10.4.1's error and quality models (`nanopore2023`), reads of 98 % identity on average, drawn 4 kb long on average (standard deviation 1.8 kb), those over 7 kb dropped (3.7 kb on average, then): 40× of `Scer`, 45× of `Spar`, 40× of `Cgla`, bleed-through included; from the genome with its shared errors; bleed-through: 15 % of `Spar`'s reads are `Cgla`'s (enough to assemble some), 1 % of the others' another species' | the sequence, strand and position |
| WGS, pooled (`pool/wgs_R{1,2}.fastq.gz`) | InSilicoSeq 2.0.1, NovaSeq's errors, 100× of each species: its held-out pairs (10 %) see a genome's QV past 40 | the sequence |
| Hi-C, pooled (`pool/hic_R{1,2}.fastq.gz`) | here: Arima's DpnII and HinfI, 400,000 pairs; each pair two restriction fragments' ends ligated, a read through the junction when the molecule is short; contacts decaying as 1/s; 12 % between chromosomes, 3 % random ligations (a partner from anywhere in the pool, by abundance: ε), 5 % unligated, 5 % PCR duplicates | both ends, their sites and enzymes |
| RNA-seq, pooled (`pool/rna-<condition>-<replicate>_R{1,2}.fastq.gz`) | here: two conditions (YPD, heat) of two replicates, 500,000 pairs of 100 bp each, reverse-stranded (dUTP); 40, 35 and 25 % of the pairs `Scer`'s, `Spar`'s and `Cgla`'s; each gene's fragments (250 bp on average) drawn by its level and its length, its level log-normal, 15 % of the genes changed in heat (a log2 fold from N(0, 1.5)), each replicate's noise 15 % | the species and the gene |

Each species' mitochondrion is 5 copies per nuclear genome in every library
but the RNA-seq. `truth.json` records every number above and where each
repeat sits; `species.csv` and `libraries.csv` are the run's samplesheets.
The truth: `references/<species>.fa` (the genome, its repeats in),
`references/<species>.gtf` (its genes, in the genome's coordinates),
`references/<species>.cds.fa` (its nuclear protein-coding genes' CDS),
`truth/<species>.tpm.tsv` (each gene's true TPM in each RNA-seq library) and
`truth/<species>.ont_errors.tsv` (the nanopore reads' shared errors).

The set is deterministic: the same `GENERATOR_VERSION` writes the same reads.
It is published, with that version, as the release `benchmark-data-v<N>`, which
[the benchmarks' workflow](../.github/workflows/benchmarks.yml) downloads, or
makes and publishes when it does not exist.

```bash
# The simulators: Badread and InSilicoSeq
micromamba create -p .bench -c conda-forge -c bioconda python=3.12 badread=0.4.1 insilicoseq=2.0.1
PATH=$PWD/.bench/bin:$PATH python benchmarks/make_benchmark_data.py --outdir benchmarks/data --threads 4
```

## A run, and its checks

```bash
nextflow run . -profile benchmark,docker --outdir results
python benchmarks/evaluate.py --results results --milestone 1
```

`evaluate.py` runs BUSCO on the references, QUAST and minimap2 against them,
and plexembler's own `junctions` (`--busco`, `--quast`, `--minimap2`,
`--plexembler`: how to run them, on PATH by default) and writes
`results/benchmark/m<N>.md`; it exits 1 when a check fails (a check with
nothing to judge reads "not judged"). Each milestone's checks are its "Done
when":

- milestone 1, the kept contigs: BUSCO against the reference's, QV,
  misassemblies;
- milestone 2, the two splits (on the contigs, then on the screened contigs)
  and their epsilon, the screened contigs and the polishing: its QV, and its
  errors against the reference, QUAST's (`--hic-results`: a second run,
  polished with the Hi-C reads, `--polish_reads hic`);
- milestone 3, each genome's chromosomes against the reference's (whole, end
  to end), its misassemblies and its flagged joins; and, on its own Hi-C pairs,
  its chromosomes cut in two and laid back (each join supported), and planted
  misjoins, an inversion and a translocation (each flagged);
- milestone 4, its reference genes found among its predicted ones (a CDS
  covering 80 % of the reference's at 95 % identity), BUSCO on its proteins
  against BUSCO on its genome;
- milestone 5, Spearman's rho between its TPM and the true TPM, in every
  RNA-seq library, on its genes paired with the reference's.

## To come

`test_full` is real data: public reads of five to ten yeasts with complete
references, mixed in silico (docs/plan.md, T2). The small set's next versions:
misjoins the scaffolder itself makes (v2's are planted where they are scored),
and a species of R9.4.1 reads, for the chemistries' mix.
