"""A milestone's "Done when" (docs/plan.md), checked on the benchmark's results.

Reads a run of the pipeline on the benchmark pool's small set (``-profile
benchmark``) and the set's truth, runs the checks the milestone names, and
writes ``<results>/benchmark/m<N>.md`` and ``m<N>.json``. It exits 1 when a
check fails; a check with nothing to judge reads "not judged" and fails nothing.

    python benchmarks/evaluate.py --results results --data benchmarks/data/small \\
        --milestone 1 --busco-db busco_downloads

Milestone 1, for each species: BUSCO's completeness of the kept contigs within
1 point of its reference's; a QV of 40 or more after polishing; no misassembly
over 10 kb against its reference (QUAST). And the pool: the pair of the genus
found as one (85 to 95 %), the distant species below.

Milestone 2: 99.5 % or more of each species' split WGS pairs its own (InSilicoSeq
names each read after its sequence: ``--save_split_reads``); epsilon within 10 %
of the simulated rate; each species' screened contigs its own (minimap2 against
the references: the planted bleed-through moved or dropped); no QV lower after
the short-read polishing; and, with ``--hic-results`` (a run with
``--polish_reads hic``), the Hi-C reads' polishing gain 80 % of the WGS pairs'
(not judged when the WGS pairs' polishing was not kept: no gain to compare).

Milestone 3, for each species: 90 % or more of its reference's chromosomes
each 90 % in one of its genome's chromosomes, and as many end to end (both ends
of the reference chromosome, within 2 kb, in that chromosome: their telomeres,
since the references are complete); no misassembly over 10 kb against its
reference (QUAST); no join flagged by ``plexembler junctions`` (the set plants
no misjoin yet).

The tools run as ``--busco``, ``--quast`` and ``--minimap2`` say: commands on PATH by default,
or in a container (``"apptainer exec busco.sif busco"``).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
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


def quast_misassemblies(
    quast: list[str], assembly: Path, reference: Path, workdir: Path
) -> int:
    """QUAST's misassemblies of an assembly against its reference, over 10 kb."""
    out = workdir / f"quast_{reference.stem}"
    command = [
        *quast, str(assembly), "-r", str(reference), "-o", str(out),
        "--extensive-mis-size", "10000", "--min-identity", "90", "-t", "4", "--fast",
    ]  # fmt: skip
    subprocess.run(
        command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    with (out / "report.tsv").open(newline="") as handle:
        report = {
            row[0]: row[1] for row in csv.reader(handle, delimiter="\t") if len(row) > 1
        }
    return int(report.get("# misassemblies", "0"))


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
        qv = float(kept["qv"]) if kept["qv"] not in ("NA", "") else 0.0
        checks.append(
            {
                "species": species,
                "check": "QV of 40 or more",
                "value": f"{qv:.1f}",
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
    results: Path, data: Path, minimap2: list[str], hic_results: Path | None
) -> list[dict]:
    truth = json.loads((data / "truth.json").read_text())
    rows = species_rows(data)
    pool = rows[0]["pool"]
    workdir = results / "benchmark"
    workdir.mkdir(parents=True, exist_ok=True)
    checks = []
    deconvolution = results / "pools" / pool / "deconvolution" / "contigs"
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
                "check": "99.5 % of its split WGS pairs its own",
                "value": f"{100 * share:.2f} % of {total}",
                "passed": share >= 0.995,
            }
        )
    numbers = json.loads((deconvolution / f"{pool}.{pool}-hic.split.json").read_text())
    expected = truth["hic"].get("random_ligation")
    found = numbers.get("epsilon")
    checks.append(
        {
            "species": pool,
            "check": "epsilon within 10 % of the simulated rate",
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
    return checks


def merged(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    found: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if found and start <= found[-1][1]:
            found[-1] = (found[-1][0], max(found[-1][1], end))
        else:
            found.append((start, end))
    return found


def chromosome_cover(
    minimap2: list[str], genome: Path, reference: Path, workdir: Path
) -> dict[str, tuple[str, float, bool]]:
    """Each reference chromosome's best genome sequence: its name, the share of
    the reference chromosome it covers, and whether it holds both its ends."""
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
        share = sum(end - start for start, end in pieces) / lengths[target]
        ends = pieces[0][0] <= 2000 and pieces[-1][1] >= lengths[target] - 2000
        if target not in best or share > best[target][1]:
            best[target] = (query, share, ends)
    return best


def milestone_3(
    results: Path, data: Path, minimap2: list[str], quast: list[str]
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
                minimap2, genome, reference, workdir
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
                    "check": "no join flagged (none planted)",
                    "value": f"{flagged} of {len(verdicts)} joins",
                    "passed": flagged == 0,
                }
            )
    return checks


def write_report(checks: list[dict], results: Path, milestone: int) -> bool:
    out = results / "benchmark"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"m{milestone}.json").write_text(json.dumps(checks, indent=2) + "\n")
    lines = [
        f"# Milestone {milestone}, on the benchmark pool",
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
        "--data", type=Path, default=Path(__file__).parent / "data" / "small"
    )
    parser.add_argument("--milestone", type=int, choices=[1, 2, 3], default=1)
    parser.add_argument("--busco", default="busco", help="how to run BUSCO")
    parser.add_argument("--quast", default="quast.py", help="how to run QUAST")
    parser.add_argument("--minimap2", default="minimap2", help="how to run minimap2")
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
            arguments.hic_results,
        )
    else:
        checks = milestone_3(
            arguments.results,
            arguments.data,
            shlex.split(arguments.minimap2),
            shlex.split(arguments.quast),
        )
    return 0 if write_report(checks, arguments.results, arguments.milestone) else 1


if __name__ == "__main__":
    sys.exit(main())
