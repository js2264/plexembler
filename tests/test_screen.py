"""`plexembler screen`: each contig kept, moved to another species, or dropped."""

import gzip
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler.cli import cli
from plexembler.rename import read_fasta
from plexembler.screen import Thresholds, merged_length, screen

GC_RICH = "GCGCATATGC" * 100  # 1 kb, 60 % GC
AT_RICH = "ATATATATGC" * 100  # 1 kb, 20 % GC

#: A pool of two species. Scer_c1, Scer_c2 and Cgla_c1, Cgla_c2 are their own;
#: Scer_mito is Scer's mitochondrion; Scer_dup is Cgla's DNA, which Cgla's
#: contigs hold already; Scer_hic is Cgla's DNA too, which they lack, and
#: Scer's reads cover thinly.
CONTIGS = {
    "Scer_c1": GC_RICH * 2,
    "Scer_c2": GC_RICH * 2,
    "Scer_mito": AT_RICH,
    "Scer_dup": GC_RICH,
    "Scer_hic": GC_RICH,
    "Cgla_c1": GC_RICH * 2,
    "Cgla_c2": GC_RICH * 2,
}


def paf_line(
    query, qlen, qstart, qend, target, tlen, tstart, tend, matches=None, primary=True
):
    block = tend - tstart
    matches = block if matches is None else matches
    tags = ["tp:A:P" if primary else "tp:A:S"]
    cells = [
        query,
        qlen,
        qstart,
        qend,
        "+",
        target,
        tlen,
        tstart,
        tend,
        matches,
        block,
        60,
        *tags,
    ]
    return "\t".join(str(c) for c in cells)


def depth(contig: str, times: float, prefix: str = "r") -> list[str]:
    """Reads covering a contig ``times`` deep: one alignment of its length each."""
    length = len(CONTIGS[contig])
    reads = []
    whole = int(times)
    for i in range(whole):
        reads.append(
            paf_line(
                f"{prefix}{contig}{i}", length, 0, length, contig, length, 0, length
            )
        )
    if times > whole:
        part = int(length * (times - whole))
        reads.append(
            paf_line(f"{prefix}{contig}x", part, 0, part, contig, length, 0, part)
        )
    return reads


@pytest.fixture
def pool(tmp_path: Path) -> dict:
    fasta = tmp_path / "pool.fa"
    fasta.write_text(
        "".join(f">{n}\n{s}\n" for n, s in CONTIGS.items()), encoding="ascii"
    )
    scer = tmp_path / "Scer.paf"
    scer.write_text(
        "\n".join(
            depth("Scer_c1", 40) + depth("Scer_c2", 40) + depth("Scer_mito", 200, "m")
            + depth("Scer_dup", 2) + depth("Scer_hic", 10)
            # A secondary alignment is not depth.
            + [paf_line("s", 2000, 0, 2000, "Scer_c1", 2000, 0, 2000, primary=False)]
        )
        + "\n",
        encoding="ascii",
    )  # fmt: skip
    cgla = tmp_path / "Cgla.paf"
    cgla.write_text(
        "\n".join(depth("Cgla_c1", 30) + depth("Cgla_c2", 30) + depth("Scer_dup", 15))
        + "\n",
        encoding="ascii",
    )
    identity = tmp_path / "self.paf"
    identity.write_text(
        "\n".join(
            [
                # Scer_dup is Cgla_c1's first kilobase, 99.9 % identical.
                paf_line(
                    "Scer_dup", 1000, 0, 1000, "Cgla_c1", 2000, 0, 1000, matches=999
                ),
                # Scer_c1 is 90 % like Cgla_c2: a sister's copy, not the same DNA.
                paf_line(
                    "Scer_c1", 2000, 0, 2000, "Cgla_c2", 2000, 0, 2000, matches=1800
                ),
                # Within a species, alignments say nothing.
                paf_line("Scer_c1", 2000, 0, 2000, "Scer_c2", 2000, 0, 2000),
            ]
        )
        + "\n",
        encoding="ascii",
    )
    contacts = tmp_path / "contacts.tsv"
    contacts.write_text(
        "contig_a\tcontig_b\tpairs\n"
        "Scer_c1\tScer_c2\t500\n"
        "Cgla_c1\tCgla_c2\t400\n"
        "Cgla_c1\tScer_hic\t60\n"
        "Scer_c1\tScer_hic\t5\n"
        "Cgla_c2\tScer_c1\t3\n",
        encoding="ascii",
    )
    return {
        "fasta": fasta,
        "reads": {"Scer": scer, "Cgla": cgla},
        "identity": identity,
        "contacts": contacts,
    }


def test_merged_length():
    assert merged_length([(0, 10), (5, 20), (30, 40), (35, 38)]) == 30
    assert merged_length([]) == 0


def test_every_verdict(pool, tmp_path: Path):
    found = screen(
        pool["fasta"],
        pool["reads"],
        pool["identity"],
        [pool["contacts"]],
        tmp_path / "out",
        "poolA",
    )
    verdicts = {name: (c.verdict, c.new_name) for name, c in found.items()}
    assert verdicts == {
        "Scer_c1": ("keep", "Scer_c1"),
        "Scer_c2": ("keep", "Scer_c2"),
        "Scer_mito": ("organelle", "Scer_chrM"),
        "Scer_dup": ("drop", "Scer_dup"),
        "Scer_hic": ("move", "Cgla_Scer_hic"),
        "Cgla_c1": ("keep", "Cgla_c1"),
        "Cgla_c2": ("keep", "Cgla_c2"),
    }
    assert found["Scer_dup"].found_in == {"Cgla": 1.0}
    assert found["Scer_c1"].depth["Scer"] == 40
    assert "60 Hi-C contacts with Cgla" in found["Scer_hic"].reason
    with gzip.open(tmp_path / "out" / "Scer.screened.fa.gz", "rt") as handle:
        assert [name for name, _ in read_fasta(handle)] == [
            "Scer_c1",
            "Scer_c2",
            "Scer_chrM",
        ]
    with gzip.open(tmp_path / "out" / "Cgla.screened.fa.gz", "rt") as handle:
        assert [name for name, _ in read_fasta(handle)] == [
            "Cgla_Scer_hic",
            "Cgla_c1",
            "Cgla_c2",
        ]
    table = (tmp_path / "out" / "poolA.screen.tsv").read_text().splitlines()
    assert table[0].startswith(
        "contig\tspecies\tlength_bp\tgc_percent\tdepth_Cgla\tdepth_Scer"
    )
    multiqc = (tmp_path / "out" / "poolA.screen_mqc.tsv").read_text().splitlines()
    assert "Scer\t2\t0\t1\t1\t1" in multiqc
    assert "Cgla\t2\t1\t0\t0\t0" in multiqc


def test_moved_by_depth_without_hic(pool, tmp_path: Path):
    """Without Hi-C, a contig its own species' reads barely cover, and another
    species' reads cover, goes to that species."""
    pool["reads"]["Scer"].write_text(
        "\n".join(depth("Scer_c1", 40) + depth("Scer_c2", 40) + depth("Scer_hic", 1))
        + "\n"
    )
    pool["reads"]["Cgla"].write_text(
        "\n".join(depth("Cgla_c1", 30) + depth("Cgla_c2", 30) + depth("Scer_hic", 25))
        + "\n"
    )
    found = screen(pool["fasta"], pool["reads"], None, [], tmp_path / "out", "poolA")
    assert (found["Scer_hic"].verdict, found["Scer_hic"].new_name) == (
        "move",
        "Cgla_Scer_hic",
    )
    # Thin on its own reads, but no other species' reads cover it: kept.
    assert found["Scer_dup"].verdict == "keep"


def test_a_second_organelle_and_thresholds(pool, tmp_path: Path):
    found = screen(
        pool["fasta"],
        pool["reads"],
        pool["identity"],
        [pool["contacts"]],
        tmp_path / "out",
        "poolA",
        Thresholds(min_contacts=100),
    )
    # Too few contacts to move it on Hi-C; its own reads, thin, not thin enough
    # to move it on depth.
    assert found["Scer_hic"].verdict == "keep"


def test_hic_alone_moves_no_well_covered_contig(pool, tmp_path: Path):
    """A species of one contig has no contacts of its own; random ligation
    gives it some with the others. Its own reads cover it: it stays."""
    pool["contacts"].write_text(
        "contig_a\tcontig_b\tpairs\nCgla_c2\tScer_c1\t50\n", encoding="ascii"
    )
    found = screen(
        pool["fasta"],
        pool["reads"],
        None,
        [pool["contacts"]],
        tmp_path / "out",
        "poolA",
    )
    assert found["Cgla_c2"].contacts == {"Scer": 50}
    assert found["Cgla_c2"].verdict == "keep"
    assert found["Scer_c1"].verdict == "keep"


def test_a_moved_contig_keeps_its_name(pool, tmp_path: Path):
    """Moved, a contig is named after the species it goes to, then its own
    name: never one of that species' contigs' names."""
    pool["fasta"].write_text(
        pool["fasta"].read_text().replace(">Scer_hic", ">Scer_c3"), encoding="ascii"
    )
    for path in (*pool["reads"].values(), pool["contacts"]):
        path.write_text(path.read_text().replace("Scer_hic", "Scer_c3"))
    found = screen(
        pool["fasta"],
        pool["reads"],
        pool["identity"],
        [pool["contacts"]],
        tmp_path / "out",
        "poolA",
    )
    assert found["Scer_c3"].new_name == "Cgla_Scer_c3"


def test_the_command(pool, tmp_path: Path):
    outdir = tmp_path / "out"
    arguments = [
        "screen",
        str(pool["fasta"]),
        "-r",
        f"Scer={pool['reads']['Scer']}",
        "-r",
        f"Cgla={pool['reads']['Cgla']}",
        "-i",
        str(pool["identity"]),
        "-c",
        str(pool["contacts"]),
        "-p",
        "poolA",
        "-o",
        str(outdir),
    ]
    result = CliRunner().invoke(cli, arguments)
    assert result.exit_code == 0, result.output
    log = next(outdir.glob("poolA.plexembler_*.log")).read_text()
    assert "7 contigs: 1 drop, 4 keep, 1 move, 1 organelle" in log
    assert "Scer_hic: move (Cgla_Scer_hic)" in log
    assert CliRunner().invoke(cli, arguments).exit_code == 1
    bad = CliRunner().invoke(cli, [*arguments[:3], "Scer", *arguments[4:], "--force"])
    assert bad.exit_code == 2
    assert "not SPECIES=PAF" in bad.output
    missing = CliRunner().invoke(
        cli, [*arguments[:3], "Scer=nowhere.paf", *arguments[4:], "--force"]
    )
    assert missing.exit_code == 2


def test_an_alignment_counts_for_both_contigs(pool, tmp_path: Path):
    """minimap2's -X reports each pair of contigs once: Scer_dup found in Cgla
    though Cgla_c1 is the query."""
    pool["identity"].write_text(
        paf_line("Cgla_c1", 2000, 0, 1000, "Scer_dup", 1000, 0, 1000, matches=999)
        + "\n",
        encoding="ascii",
    )
    found = screen(
        pool["fasta"],
        pool["reads"],
        pool["identity"],
        [pool["contacts"]],
        tmp_path / "out",
        "poolA",
    )
    assert found["Scer_dup"].found_in == {"Cgla": 1.0}
    assert found["Scer_dup"].verdict == "drop"
    assert found["Cgla_c1"].found_in == {"Scer": 0.5}


def test_a_duplicate_another_species_covers_deeper(pool, tmp_path: Path):
    """Another species' organelle, assembled from bleed-through: its own reads
    at 0.4 of their depth (not thin enough on their own), the other species'
    at 2.7 of theirs."""
    pool["reads"]["Scer"].write_text(
        "\n".join(depth("Scer_c1", 40) + depth("Scer_c2", 40) + depth("Scer_dup", 16))
        + "\n"
    )
    pool["reads"]["Cgla"].write_text(
        "\n".join(depth("Cgla_c1", 30) + depth("Cgla_c2", 30) + depth("Scer_dup", 81))
        + "\n"
    )
    found = screen(
        pool["fasta"],
        pool["reads"],
        pool["identity"],
        [],
        tmp_path / "out",
        "poolA",
    )
    assert found["Scer_dup"].verdict == "drop"
    assert "Cgla's at 2.70" in found["Scer_dup"].reason
