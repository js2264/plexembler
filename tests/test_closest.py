"""`plexembler closest`: a species' closest known genome, and its sequence.

NCBI is played by a local HTTP server, and a node without network by a port
nothing listens on: the tests never leave the machine.
"""

import gzip
import io
import socket
import threading
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest
from click.testing import CliRunner

from plexembler import closest as closest_module
from plexembler.cli import cli
from plexembler.closest import (
    Match,
    download_url,
    fetch_genome,
    read_gather,
    read_lineages,
)

#: sourmash gather's columns (4.9), those the command reads among them.
COLUMNS = (
    "intersect_bp,f_orig_query,f_match,name,md5,gather_result_rank,query_name,"
    "query_containment_ani,match_containment_ani,average_containment_ani,max_containment_ani"
)
SCER = "2712000,0.80,0.24,GCA_000146045.2 Saccharomyces cerevisiae S288C,8243f891,0,Q,0.993,0.955,0.974,0.993"
SPAR = "1200000,0.30,0.10,GCF_000166955.1 Saccharomyces paradoxus,a1b2c3d4,1,Q,0.91,0.88,0.895,0.91"
CUSTOM = "900000,0.20,0.20,my own yeast,00000000,0,Q,0.90,0.90,0.90,0.90"

GENOME = ">BK006935.2 TPA_inf: Saccharomyces cerevisiae S288C chromosome I, complete sequence\nACGTACGTAC\nGTACGT\n"


def gather(tmp_path: Path, rows: list[str], name: str = "Scer.csv.gz") -> Path:
    path = tmp_path / name
    with gzip.open(path, "wt", encoding="utf-8") as sink:
        sink.write("\n".join([COLUMNS, *rows]) + "\n")
    return path


def package(accession: str, fasta: str | None = GENOME) -> bytes:
    """NCBI Datasets' package of an assembly: its sequence, and the rest."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as bundle:
        bundle.writestr("README.md", "NCBI Datasets\n")
        bundle.writestr("ncbi_dataset/data/assembly_data_report.jsonl", "{}\n")
        if fasta is not None:
            bundle.writestr(
                f"ncbi_dataset/data/{accession}/{accession}_R64_genomic.fna", fasta
            )
    return buffer.getvalue()


#: What the fake NCBI answers, by accession: a status and a body.
ANSWERS = {
    "GCA_000146045.2": (200, package("GCA_000146045.2")),
    "GCF_000146045.2": (200, package("GCF_000146045.2")),
    "GCA_000000001.1": (200, b"PK\x03\x04 not a whole zip"),
    "GCA_000000002.1": (200, package("GCA_000000002.1", fasta=None)),
    "GCA_000000003.1": (404, b"not found"),
    "GCA_000000004.1": (503, b"busy"),
    "GCA_000000005.1": (403, b"forbidden"),
}


class FakeNCBI(BaseHTTPRequestHandler):
    def do_GET(self):
        self.server.paths.append(self.path)  # type: ignore[attr-defined]
        accession = self.path.split("/genome/accession/", 1)[-1].split("/", 1)[0]
        status, body = ANSWERS.get(accession, (404, b"unknown"))
        self.send_response(status)
        self.send_header("Content-Type", "application/zip")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *_args):
        pass


@pytest.fixture
def ncbi(monkeypatch):
    """The fake NCBI's API, reached directly, whatever proxy the machine has."""
    monkeypatch.setenv("no_proxy", "*")
    monkeypatch.setenv("NO_PROXY", "*")
    monkeypatch.setattr(closest_module.time, "sleep", lambda _seconds: None)
    server = ThreadingHTTPServer(("127.0.0.1", 0), FakeNCBI)
    server.paths = []  # type: ignore[attr-defined]
    thread = threading.Thread(
        target=server.serve_forever, kwargs={"poll_interval": 0.02}, daemon=True
    )
    thread.start()
    yield server, f"http://127.0.0.1:{server.server_address[1]}/datasets/v2"
    server.shutdown()
    server.server_close()


@pytest.fixture
def offline(monkeypatch):
    """An API where nothing answers: a node without network."""
    monkeypatch.setenv("no_proxy", "*")
    monkeypatch.setenv("NO_PROXY", "*")
    monkeypatch.setattr(closest_module.time, "sleep", lambda _seconds: None)
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    return f"http://127.0.0.1:{port}/datasets/v2"


def run(*arguments: str):
    return CliRunner().invoke(cli, ["closest", *arguments])


def test_read_gather_best_first(tmp_path: Path):
    matches = read_gather(gather(tmp_path, [SPAR, SCER]))
    assert [m.accession for m in matches] == ["GCA_000146045.2", "GCF_000166955.1"]
    best = matches[0]
    assert best.rank == 1
    assert best.ani == pytest.approx(99.3)
    assert best.species_containment == pytest.approx(80.0)
    assert best.match_containment == pytest.approx(24.0)
    assert best.shared_bp == 2712000
    assert best.title == "Saccharomyces cerevisiae S288C"


def test_read_gather_without_matches_or_columns(tmp_path: Path):
    empty = tmp_path / "empty.csv.gz"
    with gzip.open(empty, "wt") as sink:
        sink.write("\n")  # sourmash gather's --create-empty-results, or a stub's
    assert read_gather(empty) == []
    bad = tmp_path / "bad.csv"
    bad.write_text("name,md5\nx,y\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no column intersect_bp"):
        read_gather(bad)


def test_a_name_without_accession(tmp_path: Path):
    (match,) = read_gather(gather(tmp_path, [CUSTOM]))
    assert match.accession is None
    assert match.title == "my own yeast"
    assert Match(1, "GCA_1", None, None, None, None, 0).title == "GCA_1"


def test_read_lineages_by_ranks_or_lineage(tmp_path: Path):
    ranks = tmp_path / "lineages.csv.gz"
    with gzip.open(ranks, "wt", encoding="utf-8") as sink:
        sink.write(
            "ident,superkingdom,phylum,class,order,family,genus,species,strain\n"
            "GCA_000146045,Eukaryota,Ascomycota,Saccharomycetes,Saccharomycetales,"
            "Saccharomycetaceae,Saccharomyces,Saccharomyces cerevisiae,\n"
            "GCA_000000009.1,Eukaryota,,,,,,,\n"
        )
    found = read_lineages(ranks, {"GCA_000146045.2"})
    assert found == {
        "GCA_000146045.2": "Eukaryota;Ascomycota;Saccharomycetes;Saccharomycetales;"
        "Saccharomycetaceae;Saccharomyces;Saccharomyces cerevisiae"
    }
    joined = tmp_path / "lineages.csv"
    joined.write_text(
        "accession,lineage\nGCA_000146045.2,d__Eukaryota; p__Ascomycota\n",
        encoding="utf-8",
    )
    assert read_lineages(joined, {"GCA_000146045.2"}) == {
        "GCA_000146045.2": "d__Eukaryota;p__Ascomycota"
    }
    wrong = tmp_path / "wrong.csv"
    wrong.write_text("name,genus\nx,y\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not a sourmash taxonomy"):
        read_lineages(wrong, {"GCA_000146045.2"})


def test_the_download_url():
    assert download_url("https://api/v2/", "GCA_1.1") == (
        "https://api/v2/genome/accession/GCA_1.1/download"
        "?include_annotation_type=GENOME_FASTA"
    )


def test_fetch_a_genome(ncbi, tmp_path: Path):
    server, api = ncbi
    output = tmp_path / "Scer.closest.fa.gz"
    fetched = fetch_genome("GCA_000146045.2", output, api=api)
    assert fetched.status == "downloaded"
    assert "GCA_000146045.2_R64_genomic.fna" in fetched.detail
    with gzip.open(output, "rt") as handle:
        assert handle.read() == GENOME
    assert server.paths == [
        "/datasets/v2/genome/accession/GCA_000146045.2/download"
        "?include_annotation_type=GENOME_FASTA"
    ]


@pytest.mark.parametrize(
    ("accession", "status", "detail", "tries"),
    [
        ("GCA_000000001.1", "not_found", "NCBI sent no genome", 1),
        ("GCA_000000002.1", "not_found", "has no sequence", 1),
        ("GCA_000000003.1", "not_found", "HTTP 404", 1),
        ("GCA_000000004.1", "failed", "HTTP 503", 3),
        ("GCA_000000005.1", "failed", "HTTP 403", 1),
    ],
)
def test_what_ncbi_may_answer(ncbi, tmp_path: Path, accession, status, detail, tries):
    server, api = ncbi
    output = tmp_path / "x.fa.gz"
    fetched = fetch_genome(accession, output, api=api, attempts=3)
    assert fetched.status == status
    assert detail in fetched.detail
    assert len(server.paths) == tries
    assert not output.exists()


def test_a_broken_package(tmp_path: Path, monkeypatch):
    def broken(*_args, **_kwargs):
        raise zipfile.BadZipFile("Bad CRC-32")

    monkeypatch.setattr(closest_module.shutil, "copyfileobj", broken)
    output = tmp_path / "x.fa.gz"
    fetched = closest_module.unpack(
        io.BytesIO(package("GCA_000146045.2")), "GCA_000146045.2", output
    )
    assert fetched.status == "failed"
    assert "is broken" in fetched.detail
    assert not output.exists()
    assert not list(tmp_path.iterdir())


def test_fetch_without_network(offline, tmp_path: Path):
    fetched = fetch_genome(
        "GCA_000146045.2", tmp_path / "x.fa.gz", api=offline, attempts=2
    )
    assert fetched.status == "offline"
    assert "could not be reached" in fetched.detail


def test_fetch_when_the_connection_drops(ncbi, tmp_path: Path, monkeypatch):
    _server, api = ncbi

    def dropped(*_args, **_kwargs):
        raise TimeoutError("timed out")

    monkeypatch.setattr(closest_module, "_receive", dropped)
    fetched = fetch_genome("GCA_000146045.2", tmp_path / "x.fa.gz", api=api, attempts=1)
    assert fetched.status == "offline"
    assert "timed out" in fetched.detail


def test_the_command_downloads_the_closest_genome(ncbi, tmp_path: Path):
    _server, api = ncbi
    lineages = tmp_path / "lineages.csv"
    lineages.write_text(
        "ident,superkingdom,genus,species\n"
        "GCA_000146045.2,Eukaryota,Saccharomyces,Saccharomyces cerevisiae\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "out"
    result = run(
        str(gather(tmp_path, [SCER, SPAR])),
        "-p", "Scer", "-o", str(outdir), "-l", str(lineages), "--api", api,
    )  # fmt: skip
    assert result.exit_code == 0, result.output
    table = (outdir / "Scer.closest.tsv").read_text().splitlines()
    assert table[0] == (
        "rank\taccession\tname\tani\tspecies_containment\tmatch_containment"
        "\tshared_bp\tlineage\tsequence"
    )
    assert table[1] == (
        "1\tGCA_000146045.2\tGCA_000146045.2 Saccharomyces cerevisiae S288C\t99.30"
        "\t80.00\t24.00\t2712000\tEukaryota;Saccharomyces;Saccharomyces cerevisiae"
        "\tdownloaded"
    )
    assert table[2].endswith("\tNA\tNA")
    assert (outdir / "Scer.closest.fa.gz").exists()
    multiqc = (outdir / "Scer.closest_mqc.tsv").read_text()
    assert "# id: 'plexembler_closest'" in multiqc
    assert multiqc.splitlines()[-1].startswith(
        "Scer\tSaccharomyces cerevisiae S288C\tGCA_000146045.2\t99.30"
    )
    (log,) = outdir.glob("Scer.plexembler_*.log")
    assert "its closest genome is GCA_000146045.2" in log.read_text()


def test_the_command_without_network_warns_and_succeeds(offline, tmp_path: Path, plain):
    outdir = tmp_path / "out"
    result = run(
        str(gather(tmp_path, [SCER])),
        "-p", "Scer", "-o", str(outdir), "--api", offline, "--attempts", "1",
    )  # fmt: skip
    assert result.exit_code == 0, result.output
    assert not (outdir / "Scer.closest.fa.gz").exists()
    assert (
        (outdir / "Scer.closest.tsv").read_text().splitlines()[1].endswith("\toffline")
    )
    assert "Scer has no dotplot; the rest of the run goes on" in plain(result.output)
    (log,) = outdir.glob("Scer.plexembler_*.log")
    assert (
        ":: WARNING :: Scer: its closest genome's sequence was not downloaded"
        in log.read_text()
    )


def test_the_command_without_a_match(tmp_path: Path, plain):
    outdir = tmp_path / "out"
    result = run(str(gather(tmp_path, [])), "-p", "Scer", "-o", str(outdir))
    assert result.exit_code == 0, result.output
    assert (outdir / "Scer.closest.tsv").read_text().count("\n") == 1
    assert (
        (outdir / "Scer.closest_mqc.tsv")
        .read_text()
        .splitlines()[-1]
        .startswith("Scer\tno match\tNA")
    )
    assert "no genome of the database shares enough of its k-mers" in plain(
        result.output
    )


def test_the_command_without_an_accession_or_a_download(tmp_path: Path, plain):
    outdir = tmp_path / "custom"
    result = run(str(gather(tmp_path, [CUSTOM])), "-p", "Scer", "-o", str(outdir))
    assert result.exit_code == 0, result.output
    assert (
        (outdir / "Scer.closest.tsv").read_text().splitlines()[1].endswith("\tnot_ncbi")
    )
    assert "has no NCBI assembly accession" in plain(result.output)

    outdir = tmp_path / "skipped"
    result = run(
        str(gather(tmp_path, [SCER])), "-p", "Scer", "-o", str(outdir), "--no-download"
    )
    assert result.exit_code == 0, result.output
    assert (
        (outdir / "Scer.closest.tsv")
        .read_text()
        .splitlines()[1]
        .endswith("\tnot_downloaded")
    )


def test_the_command_refuses_to_overwrite_and_bad_input(tmp_path: Path, plain):
    outdir = tmp_path / "out"
    table = gather(tmp_path, [SCER])
    assert (
        run(str(table), "-p", "Scer", "-o", str(outdir), "--no-download").exit_code == 0
    )
    stale = outdir / "Scer.closest.fa.gz"
    stale.write_bytes(b"an old run's")
    result = run(str(table), "-p", "Scer", "-o", str(outdir), "--no-download")
    assert result.exit_code != 0
    assert "pass --force" in plain(result.output)
    result = run(
        str(table), "-p", "Scer", "-o", str(outdir), "--no-download", "--force"
    )
    assert result.exit_code == 0, result.output
    assert not stale.exists()

    bad = tmp_path / "bad.csv"
    bad.write_text("name\nx\n", encoding="utf-8")
    result = run(str(bad), "-p", "Scer", "-o", str(tmp_path / "bad"))
    assert result.exit_code != 0
    assert "no column intersect_bp" in plain(result.output)
    result = run(str(table), "-p", "S_cer", "-o", str(tmp_path / "x"))
    assert result.exit_code == 2
