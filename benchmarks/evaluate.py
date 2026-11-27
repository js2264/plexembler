"""A milestone's "Done when" (docs/plan.md), checked on the benchmark's results.

Reads a run of the pipeline on one of the benchmark pool's sets (``-profile
benchmark_easy``, ``benchmark_medium`` or ``benchmark_hard``) and the set's
truth, runs the checks the milestone names, and writes
``<results>/benchmark/m<N>.md`` and ``m<N>.json``. It exits 1 when a check
fails; a check with nothing to judge reads "not judged" and fails nothing.

    python benchmarks/evaluate.py --results results --data benchmarks/data/medium \\
        --milestone 1 --busco-db busco_downloads

Milestone 1, for each species: BUSCO's completeness of the kept contigs within
1 point of its reference's; a QV of 40 or more after polishing; no misassembly
over 10 kb against its reference (QUAST). And the pool: the pair of the genus
found as one (85 to 95 %), the distant species below.

Milestone 2: 99.5 % or more of each species' split WGS pairs its own (InSilicoSeq
names each read after its sequence: ``--save_split_reads``) and epsilon within
10 % of the simulated rate, in both splits of the pool (on the contigs, then on
the screened contigs); each species' screened contigs its own (minimap2 against
the references: the planted bleed-through moved or dropped); no QV lower after
the short-read polishing, nor more errors against the reference (QUAST's
mismatches and indels, what the QV guards); and, with ``--hic-results`` (a run with
``--polish_reads hic``), the Hi-C reads' polishing gain 80 % of the WGS pairs'
(not judged when the WGS pairs' polishing was not kept: no gain to compare).

Milestone 3, for each species: 90 % or more of its reference's chromosomes
each 90 % in one of its genome's chromosomes, and as many end to end (both ends
of the reference chromosome, within 2 kb, in that chromosome: their telomeres,
since the references are complete); a tandem array the set planted (Ty1
elements in tandem, Y' elements, the rDNA), which no read crosses, left out of
both, an array at a chromosome's end moving that end inside it; no misassembly over 10 kb against its
reference (QUAST); none of its scaffolds' joins flagged by ``plexembler
junctions``; and, on its own Hi-C pairs, its chromosomes cut in two and laid
back each join supported, a planted inversion and translocation flagged.

Milestone 4, for each species: 90 % or more of its reference's characterized
nuclear genes found (``truth/<species>.genes.tsv``; every gene in a set that
does not class them), a predicted CDS covering 80 % of the reference's at 95 %
identity (minimap2: the genome's coordinates are not the reference's); BUSCO on
its proteins within 2 points of BUSCO on its genome.

Milestone 5, for each species: Spearman's rho of 0.9 or more between its TPM
and the true TPM (``truth/<species>.tpm.tsv``) in every RNA-seq library, on its
genes paired with the reference's (each the other's best CDS match).

The tools run as ``--busco``, ``--quast`` and ``--minimap2`` say: commands on PATH by default,
or in a container (``"apptainer exec busco.sif busco"``).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import shlex
import subprocess
import sys
from pathlib import Path

BUSCO_LINE = re.compile(
    r"C:(?P<complete>[\d.]+)%\[S:(?P<single>[\d.]+)%,D:(?P<duplicated>[\d.]+)%\]"
)


def species_rows(data: Path) -> list[dict[str, str]]:
    with (data / "species.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def kept_row(results: Path, species: str) -> dict[str, str]:
    table = (
        results / "species" / species / "assembly" / f"{species}.assembly.selection.tsv"
    )
    with table.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["kept"] == "yes":
                return row
    raise ValueError(f"{table}: no candidate kept")


def reference_busco(
    busco: list[str],
    reference: Path,
    lineage: str,
    busco_db: Path | None,
    workdir: Path,
) -> float:
    """BUSCO's completeness of a reference, run once and kept in ``workdir``."""
    out = workdir / f"busco_{reference.stem}"
    summaries = list(out.glob("**/short_summary*.txt")) if out.exists() else []
    if not summaries:
        command = [
            *busco,
            "-i",
            str(reference),
            "-m",
            "genome",
            "-l",
            lineage,
            "-o",
            out.name,
            "--out_path",
            str(workdir),
            "-c",
            "4",
            "-f",
        ]
        if busco_db:
            command += ["--offline", "--download_path", str(busco_db)]
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
        summaries = list(out.glob("**/short_summary*.txt"))
    found = BUSCO_LINE.search(summaries[0].read_text())
    if not found:
        raise ValueError(f"{summaries[0]}: no BUSCO summary")
    return float(found["complete"])


def quast_report(
    quast: list[str], assembly: Path, reference: Path, out: Path
) -> dict[str, str]:
    """QUAST's report of an assembly against its reference, misassemblies
    over 10 kb."""
    command = [
        *quast, str(assembly), "-r", str(reference), "-o", str(out),
        "--extensive-mis-size", "10000", "--min-identity", "90", "-t", "4", "--fast",
    ]  # fmt: skip
    subprocess.run(
        command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    with (out / "report.tsv").open(newline="") as handle:
        return {
            row[0]: row[1] for row in csv.reader(handle, delimiter="\t") if len(row) > 1
        }


def quast_misassemblies(
    quast: list[str], assembly: Path, reference: Path, workdir: Path
) -> int:
    """QUAST's misassemblies of an assembly against its reference, over 10 kb."""
    report = quast_report(
        quast, assembly, reference, workdir / f"quast_{reference.stem}"
    )
    return int(report.get("# misassemblies", "0"))


def phred(rate: float) -> float:
    """An error rate per base as a QV, 60 for none."""
    return 60.0 if rate <= 0 else min(60.0, -10 * math.log10(rate))


def quast_errors(quast: list[str], assembly: Path, reference: Path, out: Path) -> float:
    """An assembly's mismatches and indels per 100 kb of its reference."""
    report = quast_report(quast, assembly, reference, out)
    return float(report.get("# mismatches per 100 kbp", "0")) + float(
        report.get("# indels per 100 kbp", "0")
    )


def milestone_1(
    results: Path, data: Path, busco: list[str], quast: list[str], busco_db: Path | None
) -> list[dict]:
    workdir = results / "benchmark"
    workdir.mkdir(parents=True, exist_ok=True)
    checks = []
    for row in species_rows(data):
        species = row["species"]
        reference = data / "references" / f"{species}.fa"
        kept = kept_row(results, species)
        assembly = (
            results / "species" / species / "assembly" / f"{species}.assembly.fa.gz"
        )
        complete = float(kept["busco_complete"])
        expected = reference_busco(
            busco, reference, row["busco_lineage"], busco_db, workdir
        )
        checks.append(
            {
                "species": species,
                "check": "BUSCO within 1 point of the reference's",
                "value": f"{complete:.1f} % (reference {expected:.1f} %)",
                "passed": complete >= expected - 1,
            }
        )
        # "QV >= 40 after polishing": the genome the polishing kept, against
        # its reference. The kept contigs' own QV (their candidates' table,
        # against the pool's WGS k-mers) is before the short reads' polishing
        # and the screen: it counts the errors the nanopore reads share, which
        # the second set plants for that polishing, and the bleed-through
        # contigs' (Spar's).
        before = float(kept["qv"]) if kept["qv"] not in ("NA", "") else None
        polished = results / "species" / species / "polish" / f"{species}.polish.fa.gz"
        errors = quast_errors(
            quast,
            polished if polished.exists() else assembly,
            reference,
            workdir / f"quast_polished_{species}",
        )
        qv = phred(errors / 100_000)
        checks.append(
            {
                "species": species,
                "check": "QV of 40 or more after polishing (against the reference)",
                "value": f"{qv:.1f}, {errors:.2f} errors per 100 kb"
                + (
                    f" (the kept contigs' own QV, before: {before:.1f})"
                    if before is not None
                    else ""
                ),
                "passed": qv >= 40,
            }
        )
        misassemblies = quast_misassemblies(quast, assembly, reference, workdir)
        checks.append(
            {
                "species": species,
                "check": "no misassembly over 10 kb (QUAST)",
                "value": str(misassemblies),
                "passed": misassemblies == 0,
            }
        )
    pool = species_rows(data)[0]["pool"]
    with (results / "pools" / pool / "ani" / f"{pool}.ani_pairs.tsv").open(
        newline=""
    ) as handle:
        verdicts = {
            (r["species_a"], r["species_b"]): r
            for r in csv.DictReader(handle, delimiter="\t")
        }
    genus = verdicts.get(("Scer", "Spar"))
    checks.append(
        {
            "species": "Scer, Spar",
            "check": "the pair of a genus found as one (85 to 95 %)",
            "value": f"{genus['ani']} % ({genus['verdict']})" if genus else "not found",
            "passed": bool(genus) and genus["verdict"] == "same genus",
        }
    )
    return checks


def polish_gain(results: Path, species: str) -> tuple[float, float] | None:
    """The draft's QV and the kept one's, from the polishing's table."""
    table = results / "species" / species / "polish" / f"{species}.polish.selection.tsv"
    if not table.exists():
        return None
    with table.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    draft = next(r for r in rows if r["candidate"] == "draft")
    kept = next(r for r in rows if r["kept"] == "yes")
    return float(draft["qv"]), float(kept["qv"])


def contig_species(
    minimap2: list[str], fasta: Path, references: Path, workdir: Path
) -> dict[str, str]:
    """Each contig's species: the reference it aligns most of its bases to."""
    paf = workdir / f"{fasta.name}.paf"
    with paf.open("w") as sink:
        subprocess.run(
            [*minimap2, "-x", "asm10", "--secondary=no", str(references), str(fasta)],
            check=True,
            stdout=sink,
            stderr=subprocess.DEVNULL,
        )
    best: dict[str, dict[str, int]] = {}
    for line in paf.read_text().splitlines():
        fields = line.split("\t")
        best.setdefault(fields[0], {})
        target = fields[5].split("_", 1)[0]
        best[fields[0]][target] = best[fields[0]].get(target, 0) + int(fields[9])
    return {
        contig: max(hits, key=lambda name: hits[name]) for contig, hits in best.items()
    }


def milestone_2(
    results: Path,
    data: Path,
    minimap2: list[str],
    quast: list[str],
    hic_results: Path | None,
) -> list[dict]:
    truth = json.loads((data / "truth.json").read_text())
    rows = species_rows(data)
    pool = rows[0]["pool"]
    workdir = results / "benchmark"
    workdir.mkdir(parents=True, exist_ok=True)
    checks = []
    # The pool is split twice: on the contigs, for the screen, then on the
    # screened contigs, for the polishing and the scaffolds.
    for stage, on in (("contigs", "contigs"), ("screened", "screened contigs")):
        deconvolution = results / "pools" / pool / "deconvolution" / stage
        for row in rows:
            species = row["species"]
            reads = deconvolution / "reads" / f"{species}.{pool}-wgs_R1.fastq.gz"
            own = total = 0
            with gzip.open(reads, "rt") as handle:
                for index, line in enumerate(handle):
                    if index % 4 == 0:
                        total += 1
                        own += line[1:].startswith(f"{species}_")
            share = own / total if total else 0.0
            checks.append(
                {
                    "species": species,
                    "check": f"99.5 % of its split WGS pairs its own (on the {on})",
                    "value": f"{100 * share:.2f} % of {total}",
                    "passed": share >= 0.995,
                }
            )
        numbers = json.loads(
            (deconvolution / f"{pool}.{pool}-hic.split.json").read_text()
        )
        expected = truth["hic"].get("random_ligation")
        found = numbers.get("epsilon")
        checks.append(
            {
                "species": pool,
                "check": f"epsilon within 10 % of the simulated rate (on the {on})",
                "value": f"{found:.4f} (simulated {expected})"
                if found is not None and expected
                else str(found),
                "passed": bool(
                    found and expected and abs(found - expected) <= 0.1 * expected
                ),
            }
        )
    every = workdir / "references.fa"
    every.write_text(
        "".join((data / "references" / f"{r['species']}.fa").read_text() for r in rows)
    )
    for row in rows:
        species = row["species"]
        contigs = results / "species" / species / "contigs" / f"{species}.contigs.fa.gz"
        found_species = contig_species(minimap2, contigs, every, workdir)
        wrong = sorted(
            contig for contig, other in found_species.items() if other != species
        )
        checks.append(
            {
                "species": species,
                "check": "its screened contigs all its own",
                "value": ", ".join(wrong) or "all",
                "passed": not wrong,
            }
        )
        gain = polish_gain(results, species)
        if gain is not None:
            checks.append(
                {
                    "species": species,
                    "check": "no QV lower after polishing",
                    "value": f"{gain[0]:.1f} -> {gain[1]:.1f}",
                    "passed": gain[1] >= gain[0],
                }
            )
        if hic_results is not None and gain is not None:
            hic = polish_gain(hic_results, species)
            wgs_gain = gain[1] - gain[0]
            hic_gain = (hic[1] - hic[0]) if hic else 0.0
            # A polishing is kept only when it raises the QV by more than the
            # choice's tolerance: without one from the WGS pairs, there is no
            # gain to hold the Hi-C reads' to.
            passed = None if wgs_gain <= 0 else hic_gain >= 0.8 * wgs_gain
            checks.append(
                {
                    "species": species,
                    "check": "Hi-C's polishing gain 80 % of WGS's",
                    "value": f"{hic_gain:.1f} vs {wgs_gain:.1f}",
                    "passed": passed,
                }
            )
        # What the QV guards, against the truth: the polished genome's errors
        # against its reference, the draft's (its screened contigs) before.
        reference = data / "references" / f"{species}.fa"
        runs = [("", results)] + ([(" with Hi-C", hic_results)] if hic_results else [])
        draft = None
        for how, run in runs:
            kept = run / "species" / species / "polish" / f"{species}.polish.fa.gz"
            if not kept.exists():
                continue
            if draft is None:
                draft = quast_errors(
                    quast, contigs, reference, workdir / f"quast_draft_{species}"
                )
            name = "hic" if how else "wgs"
            after = quast_errors(
                quast, kept, reference, workdir / f"quast_{name}_{species}"
            )
            checks.append(
                {
                    "species": species,
                    "check": f"no more errors against its reference after polishing{how}",
                    "value": f"{draft:.2f} -> {after:.2f} per 100 kb (QUAST)",
                    "passed": after <= draft,
                }
            )
    return checks


def merged(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    found: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if found and start <= found[-1][1]:
            found[-1] = (found[-1][0], max(found[-1][1], end))
        else:
            found.append((start, end))
    return found


def tandem_arrays(data: Path, species: str) -> dict[str, list[tuple[int, int]]]:
    """The tandem arrays the set planted in a species' genome, by chromosome:
    0-based, half-open."""
    truth = json.loads((data / "truth.json").read_text())
    arrays: dict[str, list[tuple[int, int]]] = {}
    for repeat in truth["species"][species].get("repeats", []):
        if repeat.get("tandem"):
            arrays.setdefault(repeat["chromosome"], []).append(
                (repeat["start"] - 1, repeat["end"])
            )
    return arrays


def chromosome_cover(
    minimap2: list[str],
    genome: Path,
    reference: Path,
    workdir: Path,
    arrays: dict[str, list[tuple[int, int]]] | None = None,
) -> dict[str, tuple[str, float, bool]]:
    """Each reference chromosome's best genome sequence: its name, the share of
    the reference chromosome it covers, and whether it holds both its ends.
    The ``arrays`` (0-based, half-open), which no read crosses, are left out of
    the share; one within 2 kb of an end moves that end inside it."""
    paf = workdir / f"{genome.name}.vs_reference.paf"
    with paf.open("w") as sink:
        subprocess.run(
            [*minimap2, "-x", "asm5", "--secondary=no", str(reference), str(genome)],
            check=True,
            stdout=sink,
            stderr=subprocess.DEVNULL,
        )
    lengths: dict[str, int] = {}
    covered: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for line in paf.read_text().splitlines():
        fields = line.split("\t")
        target, length = fields[5], int(fields[6])
        lengths[target] = length
        covered.setdefault((target, fields[0]), []).append(
            (int(fields[7]), int(fields[8]))
        )
    best: dict[str, tuple[str, float, bool]] = {}
    for (target, query), intervals in covered.items():
        pieces = merged(intervals)
        holes = merged((arrays or {}).get(target, []))
        length = lengths[target]
        inside = sum(
            max(0, min(end, hole_end) - max(start, hole_start))
            for start, end in pieces
            for hole_start, hole_end in holes
        )
        kept = length - sum(end - start for start, end in holes)
        share = (sum(end - start for start, end in pieces) - inside) / kept
        first, last = 0, length
        for hole_start, hole_end in holes:
            if hole_start <= first + 2000:
                first = max(first, hole_end)
        for hole_start, hole_end in reversed(holes):
            if hole_end >= last - 2000:
                last = min(last, hole_start)
        ends = pieces[0][0] <= first + 2000 and pieces[-1][1] >= last - 2000
        if target not in best or share > best[target][1]:
            best[target] = (query, share, ends)
    return best


#: A piece of a chromosome in a planted layout: its name, 0-based half-open
#: bounds, and whether it is laid reversed.
Piece = tuple[str, int, int, bool]


def planted_layouts(lengths: dict[str, int]) -> dict[str, list[list[Piece]]]:
    """Scaffolds of a genome's chromosomes, laid out three ways: each
    chromosome cut at its middle and laid back (true joins); the first one's
    second half reversed (an inversion); the first one's first half before the
    second one's second half (a translocation)."""
    names = list(lengths)
    half = {name: lengths[name] // 2 for name in names}
    first = names[0]
    layouts = {
        "true": [
            [(n, 0, half[n], False), (n, half[n], lengths[n], False)] for n in names
        ],
        "inversion": [
            [(first, 0, half[first], False), (first, half[first], lengths[first], True)]
        ],
    }
    if len(names) > 1:
        second = names[1]
        layouts["translocation"] = [
            [
                (first, 0, half[first], False),
                (second, half[second], lengths[second], False),
            ]
        ]
    return layouts


def write_layout(
    layout: list[list[Piece]], pairs: Path, agp: Path, moved: Path
) -> None:
    """The layout's AGP, and the pairs within it in its scaffolds' coordinates
    (a pair with a mate outside it left out)."""
    places: dict[str, list[tuple[int, int, str, int, bool]]] = {}
    with agp.open("w", encoding="ascii") as sink:
        for index, pieces in enumerate(layout, 1):
            scaffold = f"scaffold_{index}"
            offset = 0
            for part, (name, start, end, reverse) in enumerate(pieces, 1):
                length = end - start
                strand = "-" if reverse else "+"
                sink.write(
                    f"{scaffold}\t{offset + 1}\t{offset + length}\t{part}\tW\t"
                    f"{name}:{start + 1}-{end}\t1\t{length}\t{strand}\n"
                )
                places.setdefault(name, []).append(
                    (start, end, scaffold, offset, reverse)
                )
                offset += length

    def where(name: str, position: int) -> tuple[str, int] | None:
        """A 1-based position's scaffold and place in it."""
        for start, end, scaffold, offset, reverse in places.get(name, []):
            if start < position <= end:
                return scaffold, (
                    offset + end - position + 1
                    if reverse
                    else offset + position - start
                )
        return None

    with gzip.open(pairs, "rt") as handle, moved.open("w", encoding="ascii") as sink:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            one = where(fields[1], int(fields[2]))
            two = where(fields[3], int(fields[4]))
            if one and two:
                fields[1], fields[2] = one[0], str(one[1])
                fields[3], fields[4] = two[0], str(two[1])
                sink.write("\t".join(fields) + "\n")


def planted_verdicts(
    plexembler: list[str], results: Path, species: str, workdir: Path
) -> dict[str, list[str]]:
    """``plexembler junctions``'s verdicts on each planted layout of a species'
    chromosomes, from its own pairs."""
    genome = results / "species" / species / "genome"
    with (genome / f"{species}.names.tsv").open(newline="") as handle:
        chromosomes = [
            r["name"]
            for r in csv.DictReader(handle, delimiter="\t")
            if r["role"] == "chromosome"
        ]
    lengths = {}
    with (genome / f"{species}.chrom.sizes").open() as handle:
        for line in handle:
            name, length = line.split()[:2]
            if name in chromosomes:
                lengths[name] = int(length)
    pairs = results / "species" / species / "hic" / f"{species}.pairs.gz"
    verdicts = {}
    for kind, layout in planted_layouts(lengths).items():
        out = workdir / f"planted_{species}_{kind}"
        out.mkdir(parents=True, exist_ok=True)
        write_layout(layout, pairs, out / "layout.agp", out / "layout.pairs")
        subprocess.run(
            [
                *plexembler, "junctions", str(out / "layout.pairs"),
                "--agp", str(out / "layout.agp"), "--prefix", species,
                "--outdir", str(out), "--force", "--quiet",
            ],
            check=True, stdout=subprocess.DEVNULL,
        )  # fmt: skip
        with (out / f"{species}.junctions.tsv").open(newline="") as handle:
            verdicts[kind] = [
                r["verdict"] for r in csv.DictReader(handle, delimiter="\t")
            ]
    return verdicts


def milestone_3(
    results: Path,
    data: Path,
    minimap2: list[str],
    quast: list[str],
    plexembler: list[str],
) -> list[dict]:
    workdir = results / "benchmark"
    workdir.mkdir(parents=True, exist_ok=True)
    checks = []
    for row in species_rows(data):
        species = row["species"]
        reference = data / "references" / f"{species}.fa"
        genome = results / "species" / species / "genome" / f"{species}.fa.gz"
        names = results / "species" / species / "genome" / f"{species}.names.tsv"
        with names.open(newline="") as handle:
            roles = {
                r["name"]: r["role"] for r in csv.DictReader(handle, delimiter="\t")
            }
        cover = {
            target: found
            for target, found in chromosome_cover(
                minimap2, genome, reference, workdir, tandem_arrays(data, species)
            ).items()
            if not target.endswith("_chrM")
        }
        nuclear = [
            name
            for name in re.findall(r"^>(\S+)", reference.read_text(), re.MULTILINE)
            if not name.endswith("_chrM")
        ]
        whole = [
            t
            for t in nuclear
            if t in cover
            and cover[t][1] >= 0.9
            and roles.get(cover[t][0]) == "chromosome"
        ]
        end_to_end = [t for t in whole if cover[t][2]]
        checks.append(
            {
                "species": species,
                "check": "90 % of its chromosomes each 90 % in one of its genome's",
                "value": f"{len(whole)} of {len(nuclear)}",
                "passed": len(whole) >= 0.9 * len(nuclear),
            }
        )
        checks.append(
            {
                "species": species,
                "check": "90 % of its chromosomes end to end",
                "value": f"{len(end_to_end)} of {len(nuclear)}",
                "passed": len(end_to_end) >= 0.9 * len(nuclear),
            }
        )
        misassemblies = quast_misassemblies(quast, genome, reference, workdir)
        checks.append(
            {
                "species": species,
                "check": "no misassembly over 10 kb (QUAST)",
                "value": str(misassemblies),
                "passed": misassemblies == 0,
            }
        )
        table = (
            results / "species" / species / "genome" / "qc" / "junctions"
            / f"{species}.junctions.tsv"
        )  # fmt: skip
        if table.exists():
            with table.open(newline="") as handle:
                verdicts = [
                    r["verdict"] for r in csv.DictReader(handle, delimiter="\t")
                ]
            flagged = verdicts.count("unsupported")
            checks.append(
                {
                    "species": species,
                    "check": "none of its scaffolds' joins flagged",
                    "value": f"{flagged} of {len(verdicts)} joins",
                    "passed": flagged == 0,
                }
            )
        planted = planted_verdicts(plexembler, results, species, workdir)
        true = planted.pop("true")
        checks.append(
            {
                "species": species,
                "check": "its chromosomes, cut in two, each join supported",
                "value": f"{true.count('supported')} of {len(true)}",
                "passed": true.count("supported") == len(true),
            }
        )
        wrong = [v for verdicts in planted.values() for v in verdicts]
        checks.append(
            {
                "species": species,
                "check": f"the planted misjoins ({', '.join(planted)}) flagged",
                "value": f"{wrong.count('unsupported')} of {len(wrong)}",
                "passed": wrong.count("unsupported") == len(wrong),
            }
        )
    return checks


def busco_summary(folder: Path, pattern: str) -> float | None:
    """BUSCO's completeness in the summary of ``folder`` matching ``pattern``."""
    summaries = sorted(folder.glob(pattern)) if folder.exists() else []
    if not summaries:
        return None
    found = BUSCO_LINE.search(summaries[0].read_text())
    return float(found["complete"]) if found else None


def cds_matches(
    minimap2: list[str], predicted: Path, reference: Path, workdir: Path
) -> list[dict]:
    """Each alignment of a predicted CDS (its gene's) to a reference one: its
    identity, and how much of each it covers."""
    paf = workdir / f"{reference.stem}.genes.paf"
    with paf.open("w") as sink:
        subprocess.run(
            [
                *minimap2, "-c", "-x", "asm20", "-k", "15", "-w", "5", "-N", "10",
                str(reference), str(predicted),
            ],
            check=True, stdout=sink, stderr=subprocess.DEVNULL,
        )  # fmt: skip
    found = []
    for line in paf.read_text().splitlines():
        f = line.split("\t")
        if f[4] != "+":
            continue
        query_length, target_length = int(f[1]), int(f[6])
        found.append(
            {
                "gene": f[0].rsplit("-T", 1)[0],
                "reference": f[5],
                "matches": int(f[9]),
                "identity": int(f[9]) / int(f[10]),
                "query_cover": (int(f[3]) - int(f[2])) / query_length,
                "reference_cover": (int(f[8]) - int(f[7])) / target_length,
            }
        )
    return found


def gene_pairs(matches: list[dict]) -> dict[str, str]:
    """Each reference gene's predicted gene, where each is the other's best
    match (95 % identity, half of each covered)."""
    good = [
        m
        for m in matches
        if m["identity"] >= 0.95
        and m["query_cover"] >= 0.5
        and m["reference_cover"] >= 0.5
    ]
    best_of_gene: dict[str, dict] = {}
    best_of_reference: dict[str, dict] = {}
    for m in good:
        if m["matches"] > best_of_gene.get(m["gene"], {"matches": -1})["matches"]:
            best_of_gene[m["gene"]] = m
        if (
            m["matches"]
            > best_of_reference.get(m["reference"], {"matches": -1})["matches"]
        ):
            best_of_reference[m["reference"]] = m
    return {
        reference: m["gene"]
        for reference, m in best_of_reference.items()
        if best_of_gene[m["gene"]]["reference"] == reference
    }


def gene_classes(data: Path, species: str) -> dict[str, str]:
    """Each reference gene's class (characterized or not), from the set's
    truth; none in a set that does not class them."""
    table = data / "truth" / f"{species}.genes.tsv"
    if not table.exists():
        return {}
    with table.open(newline="") as handle:
        return {r["gene"]: r["class"] for r in csv.DictReader(handle, delimiter="\t")}


def reference_cds(data: Path, species: str) -> Path:
    return data / "references" / f"{species}.cds.fa"


def predicted_cds(results: Path, species: str) -> Path:
    return results / "species" / species / "annotation" / f"{species}.cds.fna.gz"


def milestone_4(results: Path, data: Path, minimap2: list[str]) -> list[dict]:
    workdir = results / "benchmark"
    workdir.mkdir(parents=True, exist_ok=True)
    checks = []
    for row in species_rows(data):
        species = row["species"]
        reference = reference_cds(data, species)
        every = re.findall(r"^>(\S+)", reference.read_text(), re.MULTILINE)
        matches = cds_matches(
            minimap2, predicted_cds(results, species), reference, workdir
        )
        found = {
            m["reference"]
            for m in matches
            if m["identity"] >= 0.95 and m["reference_cover"] >= 0.8
        }
        classes = gene_classes(data, species)
        counted = (
            [g for g in every if classes.get(g) == "characterized"]
            if classes
            else every
        )
        hits = sum(g in found for g in counted)
        which = "characterized reference genes" if classes else "reference genes"
        checks.append(
            {
                "species": species,
                "check": f"90 % of its {which} found (80 % of their CDS)",
                "value": f"{hits} of {len(counted)} ({100 * hits / len(counted):.1f} %)"
                + (
                    f"; all genes {sum(g in found for g in every)} of {len(every)}"
                    if classes
                    else ""
                ),
                "passed": hits >= 0.9 * len(counted),
            }
        )
        base = results / "species" / species
        genome = busco_summary(
            base / "genome" / "qc" / "busco", f"short_summary.*.{species}.fa.txt"
        )
        proteins = busco_summary(
            base / "annotation" / "qc" / "busco",
            f"short_summary.*.{species}.proteins*.txt",
        )
        checks.append(
            {
                "species": species,
                "check": "BUSCO on its proteins within 2 points of its genome's",
                "value": f"{proteins} % vs {genome} %",
                "passed": genome is not None
                and proteins is not None
                and proteins >= genome - 2,
            }
        )
    return checks


def ranks(values: list[float]) -> list[float]:
    """Each value's rank, ties sharing their average rank."""
    order = sorted(range(len(values)), key=values.__getitem__)
    found = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            found[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return found


def spearman(x: list[float], y: list[float]) -> float:
    """Spearman's rho: Pearson's correlation of the ranks."""
    rx, ry = ranks(x), ranks(y)
    n = len(rx)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    sx = sum((a - mx) ** 2 for a in rx) ** 0.5
    sy = sum((b - my) ** 2 for b in ry) ** 0.5
    return cov / (sx * sy) if sx and sy else float("nan")


def read_tpm(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    """A TPM table's libraries, and each gene's TPM in them."""
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    header = rows[0]
    first = 2 if len(header) > 1 and header[1] == "length_bp" else 1
    return header[first:], {r[0]: [float(v) for v in r[first:]] for r in rows[1:]}


def milestone_5(results: Path, data: Path, minimap2: list[str]) -> list[dict]:
    workdir = results / "benchmark"
    workdir.mkdir(parents=True, exist_ok=True)
    checks = []
    for row in species_rows(data):
        species = row["species"]
        pairs = gene_pairs(
            cds_matches(
                minimap2,
                predicted_cds(results, species),
                reference_cds(data, species),
                workdir,
            )
        )
        labels, truth = read_tpm(data / "truth" / f"{species}.tpm.tsv")
        found_labels, found = read_tpm(
            results / "species" / species / "expression" / f"{species}.tpm.tsv"
        )
        genes = sorted(g for g in truth if g in pairs and pairs[g] in found)
        rhos = {}
        for label in labels:
            if label not in found_labels:
                continue
            column = found_labels.index(label)
            rhos[label] = spearman(
                [truth[g][labels.index(label)] for g in genes],
                [found[pairs[g]][column] for g in genes],
            )
        worst = min(rhos.values()) if rhos else float("nan")
        checks.append(
            {
                "species": species,
                "check": "Spearman's rho of 0.9 against the true TPM, in every library",
                "value": f"{worst:.3f} at worst, on {len(genes)} genes ("
                + ", ".join(f"{k} {v:.3f}" for k, v in rhos.items())
                + ")",
                "passed": len(rhos) == len(labels) and worst >= 0.9,
            }
        )
    return checks


def set_name(data: Path) -> str:
    """The set's level (easy, medium, hard), from its truth; its folder's name
    for a set without one."""
    truth = json.loads((data / "truth.json").read_text())
    return truth.get("level", {}).get("name", data.name)


def write_report(
    checks: list[dict], results: Path, milestone: int, name: str = ""
) -> bool:
    out = results / "benchmark"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"m{milestone}.json").write_text(json.dumps(checks, indent=2) + "\n")
    on = f"the benchmark pool's {name} set" if name else "the benchmark pool"
    lines = [
        f"# Milestone {milestone}, on {on}",
        "",
        "| Species | Check | Value | Passed |",
        "|---|---|---|---|",
    ]
    verdicts = {True: "yes", False: "**no**", None: "not judged"}
    lines += [
        f"| {c['species']} | {c['check']} | {c['value']} | {verdicts[c['passed']]} |"
        for c in checks
    ]
    (out / f"m{milestone}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return all(c["passed"] is not False for c in checks)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--results", type=Path, required=True, help="the pipeline's --outdir"
    )
    parser.add_argument(
        "--data", type=Path, default=Path(__file__).parent / "data" / "medium"
    )
    parser.add_argument("--milestone", type=int, choices=[1, 2, 3, 4, 5], default=1)
    parser.add_argument("--busco", default="busco", help="how to run BUSCO")
    parser.add_argument("--quast", default="quast.py", help="how to run QUAST")
    parser.add_argument("--minimap2", default="minimap2", help="how to run minimap2")
    parser.add_argument(
        "--plexembler", default="plexembler", help="how to run plexembler (milestone 3)"
    )
    parser.add_argument(
        "--busco-db", type=Path, help="BUSCO's lineages, downloaded beforehand"
    )
    parser.add_argument(
        "--hic-results",
        type=Path,
        help="milestone 2: the --outdir of a run with --polish_reads hic",
    )
    arguments = parser.parse_args()
    if arguments.milestone == 1:
        checks = milestone_1(
            arguments.results,
            arguments.data,
            shlex.split(arguments.busco),
            shlex.split(arguments.quast),
            arguments.busco_db,
        )
    elif arguments.milestone == 2:
        checks = milestone_2(
            arguments.results,
            arguments.data,
            shlex.split(arguments.minimap2),
            shlex.split(arguments.quast),
            arguments.hic_results,
        )
    elif arguments.milestone == 3:
        checks = milestone_3(
            arguments.results,
            arguments.data,
            shlex.split(arguments.minimap2),
            shlex.split(arguments.quast),
            shlex.split(arguments.plexembler),
        )
    elif arguments.milestone == 4:
        checks = milestone_4(
            arguments.results, arguments.data, shlex.split(arguments.minimap2)
        )
    else:
        checks = milestone_5(
            arguments.results, arguments.data, shlex.split(arguments.minimap2)
        )
    name = set_name(arguments.data)
    return (
        0 if write_report(checks, arguments.results, arguments.milestone, name) else 1
    )


if __name__ == "__main__":
    sys.exit(main())
