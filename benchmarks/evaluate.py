"""A milestone's "Done when" (docs/plan.md), checked on the benchmark's results.

Reads a run of the pipeline on the benchmark pool's small set (``-profile
benchmark``) and the set's truth, runs the checks the milestone names, and
writes ``<results>/benchmark/m<N>.md`` and ``m<N>.json``. It exits 1 when a
check fails.

    python benchmarks/evaluate.py --results results --data benchmarks/data/small \\
        --milestone 1 --busco-db busco_downloads

Milestone 1, for each species: BUSCO's completeness of the kept contigs within
1 point of its reference's; a QV of 40 or more after polishing; no misassembly
over 10 kb against its reference (QUAST). And the pool: the pair of the genus
found as one (85 to 95 %), the distant species below.

The tools run as ``--busco`` and ``--quast`` say: commands on PATH by default,
or in a container (``"apptainer exec busco.sif busco"``).
"""

from __future__ import annotations

import argparse
import csv
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
    lines += [
        f"| {c['species']} | {c['check']} | {c['value']} | {'yes' if c['passed'] else '**no**'} |"
        for c in checks
    ]
    (out / f"m{milestone}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return all(c["passed"] for c in checks)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--results", type=Path, required=True, help="the pipeline's --outdir"
    )
    parser.add_argument(
        "--data", type=Path, default=Path(__file__).parent / "data" / "small"
    )
    parser.add_argument("--milestone", type=int, choices=[1], default=1)
    parser.add_argument("--busco", default="busco", help="how to run BUSCO")
    parser.add_argument("--quast", default="quast.py", help="how to run QUAST")
    parser.add_argument(
        "--busco-db", type=Path, help="BUSCO's lineages, downloaded beforehand"
    )
    arguments = parser.parse_args()
    checks = milestone_1(
        arguments.results,
        arguments.data,
        shlex.split(arguments.busco),
        shlex.split(arguments.quast),
        arguments.busco_db,
    )
    return 0 if write_report(checks, arguments.results, arguments.milestone) else 1


if __name__ == "__main__":
    sys.exit(main())
