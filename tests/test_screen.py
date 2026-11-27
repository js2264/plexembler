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
                # Within a species, two contigs of one length: neither holds
                # the other.
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


def test_the_command(pool, tmp_path: Path, plain):
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
    assert "not SPECIES=PAF" in plain(bad.output)
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


def test_a_species_most_of_whose_contigs_are_another_species(tmp_path: Path):
    """Bleed-through assembled whole, longer than the species' own contigs (the
    benchmark's easy set: all of Cgla's genome among Spar's contigs). Their
    depth is not Spar's: its own contig's is, which Spar's reads cover the most."""
    lengths = {"Spar_c1": 2000, "Spar_cg1": 3000, "Spar_cg2": 3000}
    lengths |= {"Cgla_c1": 3000, "Cgla_c2": 3000}
    fasta = tmp_path / "pool.fa"
    fasta.write_text(
        "".join(
            f">{n}\n{GC_RICH[:1000] * (n_bp // 1000)}\n" for n, n_bp in lengths.items()
        ),
        encoding="ascii",
    )

    def cover(contig: str, times: int, prefix: str) -> list[str]:
        n = lengths[contig]
        return [
            paf_line(f"{prefix}{contig}{i}", n, 0, n, contig, n, 0, n)
            for i in range(times)
        ]

    spar = tmp_path / "Spar.paf"
    spar.write_text(
        "\n".join(
            cover("Spar_c1", 38, "s")
            + cover("Spar_cg1", 3, "s")
            + cover("Spar_cg2", 3, "s")
        )
        + "\n"
    )
    cgla = tmp_path / "Cgla.paf"
    cgla.write_text(
        "\n".join(
            cover("Cgla_c1", 28, "c") + cover("Cgla_c2", 28, "c")
            + cover("Spar_cg1", 18, "c") + cover("Spar_cg2", 18, "c")
        ) + "\n"
    )  # fmt: skip
    identity = tmp_path / "self.paf"
    identity.write_text(
        "\n".join(
            paf_line(q, 3000, 0, 3000, t, 3000, 0, 3000, matches=2997)
            for q, t in (("Spar_cg1", "Cgla_c1"), ("Spar_cg2", "Cgla_c2"))
        )
        + "\n"
    )
    found = screen(
        fasta, {"Spar": spar, "Cgla": cgla}, identity, [], tmp_path / "out", "poolA"
    )
    for copy in ("Spar_cg1", "Spar_cg2"):
        assert found[copy].verdict == "drop"
        assert "Spar's reads at 0.08 of their depth" in found[copy].reason
    assert {found[c].verdict for c in ("Spar_c1", "Cgla_c1", "Cgla_c2")} == {"keep"}


def test_a_contig_its_species_holds_already(tmp_path: Path):
    """A haploid species' contig most of which a longer contig of its own
    holds, and which its reads cover less: the assembler wrote that stretch
    twice (a repeat of its graph, in a contig and alone), and the copy goes.
    Deeper than its species' contigs (a collapsed repeat), or held in part, it
    stays; and so it does for a species not known to be haploid (a diploid's
    haplotig, a hybrid's other subgenome)."""
    contigs = {
        "Scer_c1": GC_RICH * 2,
        "Scer_c2": GC_RICH * 2,
        "Scer_copy": GC_RICH,
        "Scer_repeat": GC_RICH,
        "Scer_part": GC_RICH,
    }
    fasta = tmp_path / "pool.fa"
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in contigs.items()))
    reads = []
    for name, times in [("Scer_c1", 40), ("Scer_c2", 40), ("Scer_copy", 20), ("Scer_repeat", 80), ("Scer_part", 20)]:  # fmt: skip
        length = len(contigs[name])
        reads += [
            paf_line(f"r{name}{i}", length, 0, length, name, length, 0, length)
            for i in range(times)
        ]
    paf = tmp_path / "Scer.paf"
    paf.write_text("\n".join(reads) + "\n")
    identity = tmp_path / "self.paf"
    identity.write_text(
        "\n".join(
            [
                # 80 % of Scer_copy in Scer_c1, 95 % identical; the longer
                # contig the target, then the query.
                paf_line("Scer_c1", 2000, 1000, 1500, "Scer_copy", 1000, 0, 500, matches=475),
                paf_line("Scer_copy", 1000, 500, 800, "Scer_c1", 2000, 1500, 1800, matches=285),
                paf_line("Scer_repeat", 1000, 0, 1000, "Scer_c2", 2000, 0, 1000),
                paf_line("Scer_part", 1000, 0, 300, "Scer_c2", 2000, 1000, 1300),
                # 85 % identical: not the same stretch.
                paf_line("Scer_part", 1000, 300, 1000, "Scer_c1", 2000, 0, 700, matches=595),
            ]
        )
        + "\n"
    )  # fmt: skip
    unknown = screen(fasta, {"Scer": paf}, identity, [], tmp_path / "unknown", "poolA")
    assert unknown["Scer_copy"].verdict == "keep"
    found = screen(
        fasta,
        {"Scer": paf},
        identity,
        [],
        tmp_path / "out",
        "poolA",
        haploid=frozenset({"Scer"}),
    )
    assert found["Scer_copy"].within == {"Scer_c1": pytest.approx(0.8)}
    assert found["Scer_copy"].verdict == "drop"
    assert "80 % of it in Scer_c1" in found["Scer_copy"].reason
    assert found["Scer_repeat"].verdict == "keep"
    assert found["Scer_part"].within == {"Scer_c2": pytest.approx(0.3)}
    assert found["Scer_part"].verdict == "keep"
    assert found["Scer_c1"].within == {}
    table = (tmp_path / "out" / "poolA.screen.tsv").read_text().splitlines()
    copy = next(line.split("\t") for line in table if line.startswith("Scer_copy"))
    assert copy[table[0].split("\t").index("found_in_Scer")] == "0.80"
    with gzip.open(tmp_path / "out" / "Scer.screened.fa.gz", "rt") as handle:
        assert "Scer_copy" not in [name for name, _ in read_fasta(handle)]
    # Its thresholds, from the command line.
    result = CliRunner().invoke(
        cli,
        [
            "screen", str(fasta), "-r", f"Scer={paf}", "-i", str(identity),
            "-p", "poolA", "-o", str(tmp_path / "cli"), "--haploid", "Scer",
            "--own-duplicate-depth", "0.4",
        ],
    )  # fmt: skip
    assert result.exit_code == 0, result.output
    rows = (tmp_path / "cli" / "poolA.screen.tsv").read_text()
    assert "\tkeep\tScer_copy\t" in rows
    bad = CliRunner().invoke(
        cli,
        [
            "screen", str(fasta), "-r", f"Scer={paf}", "-p", "poolA",
            "-o", str(tmp_path / "bad"), "--haploid", "S cer",
        ],
    )  # fmt: skip
    assert bad.exit_code == 2
