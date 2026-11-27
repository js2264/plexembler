"""`plexembler run`: the pipeline, run with Nextflow at the package's version.

The command hands its process over to Nextflow (`os.execve`): the tests catch
that call, and check the command and the environment it would run.
"""

import json
import os
import re
import shutil
import stat
import sysconfig
from importlib.metadata import requires
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from plexembler import __version__, run
from plexembler.cli import cli
from plexembler.run import (
    REPOSITORY,
    default_pipeline,
    environment,
    find_nextflow,
    installed_from,
    pipeline_to_run,
    split_arguments,
)

ROOT = Path(__file__).resolve().parents[1]
RELEASE = f"v{__version__}"


class Launched(Exception):
    """Raised by the fake exec in place of Nextflow."""


@pytest.fixture
def nextflow(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A Nextflow on PATH, and none beside plexembler's Python."""
    folder = tmp_path / "path"
    folder.mkdir()
    path = folder / "nextflow"
    path.write_text("#!/bin/sh\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", str(folder))
    empty = tmp_path / "scripts"
    empty.mkdir()
    monkeypatch.setattr(sysconfig, "get_path", lambda name: str(empty))
    return path


@pytest.fixture
def released(monkeypatch: pytest.MonkeyPatch) -> None:
    """plexembler as a release installs it: from an index, with no direct_url.json."""
    monkeypatch.setattr(run, "installed_from", lambda: None)


@pytest.fixture
def launched(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, list[str], dict]]:
    """The calls to os.execve, which raise Launched instead of running anything."""
    calls = []

    def execve(path, argv, env):
        calls.append((path, argv, env))
        raise Launched

    monkeypatch.setattr(os, "execve", execve)
    return calls


def invoke(*arguments: str):
    return CliRunner().invoke(cli, ["run", *arguments])


def origin(monkeypatch: pytest.MonkeyPatch, record: dict | None) -> None:
    """Install plexembler, as far as default_pipeline knows, from `record`."""
    monkeypatch.setattr(run, "installed_from", lambda: record)


# The command line


@pytest.mark.usefixtures("released")
def test_runs_the_release_of_its_version(nextflow: Path, launched):
    result = invoke("-profile", "docker", "--input", "species.csv", "--outdir", "out")
    assert isinstance(result.exception, Launched), result.output
    path, argv, _ = launched[0]
    assert path == str(nextflow)
    assert argv == [
        str(nextflow),
        "run",
        REPOSITORY,
        "-r",
        RELEASE,
        "-profile",
        "docker",
        "--input",
        "species.csv",
        "--outdir",
        "out",
    ]
    assert f"Running {nextflow} run {REPOSITORY} -r {RELEASE}" in result.stderr


@pytest.mark.usefixtures("released")
def test_nextflow_options_go_through_as_they_are(nextflow: Path, launched):
    # click would read -with-docker as -w -i -t -h…: -h would print the help.
    arguments = [
        "-profile",
        "test,docker",
        "-with-docker",
        "-queue-size",
        "2",
        "-resume",
        "-c",
        "my.config",
        "-params-file",
        "params.yaml",
        "--input=species.csv",
        "--min_chromosome",
        "-1",
        "--multiqc_title",
        "two words",
    ]
    invoke(*arguments)
    assert launched[0][1][5:] == arguments


@pytest.mark.usefixtures("released")
def test_after_a_double_dash_every_argument_is_nextflows(nextflow: Path, launched):
    invoke("--dry-run", "--", "--help", "--pipeline", "elsewhere", "-h", "--")
    assert not launched
    invoke("--", "--help", "--pipeline", "elsewhere", "-h", "--dry-run")
    argv = launched[0][1]
    assert argv[:3] == [str(nextflow), "run", REPOSITORY]
    assert argv[5:] == ["--help", "--pipeline", "elsewhere", "-h", "--dry-run"]


@pytest.mark.usefixtures("released", "nextflow")
def test_a_revision_given_is_the_one_run(launched):
    invoke("-r", "main", "-profile", "docker")
    assert launched[0][1][2:] == [REPOSITORY, "-r", "main", "-profile", "docker"]
    invoke("-revision", "0123abc")
    assert launched[1][1][2:] == [REPOSITORY, "-revision", "0123abc"]


@pytest.mark.usefixtures("released", "nextflow")
def test_a_pipeline_folder_runs_as_it_is(tmp_path: Path, launched):
    invoke("--pipeline", str(tmp_path), "-profile", "test")
    assert launched[0][1][2:] == [str(tmp_path), "-profile", "test"]


@pytest.mark.usefixtures("released", "nextflow")
def test_a_pipeline_repository_runs_at_the_version(launched):
    invoke("--pipeline=someone/plexembler", "-resume")
    assert launched[0][1][2:] == ["someone/plexembler", "-r", RELEASE, "-resume"]


@pytest.mark.usefixtures("released")
def test_a_dry_run_prints_the_command_and_runs_nothing(nextflow: Path, launched):
    result = invoke("--dry-run", "-profile", "docker", "--input", "my species.csv")
    assert result.exit_code == 0, result.output
    assert not launched
    assert result.stdout == (
        f"{nextflow} run {REPOSITORY} -r {RELEASE} "
        "-profile docker --input 'my species.csv'\n"
    )


@pytest.mark.usefixtures("released")
def test_without_nextflow_it_fails_and_says_how_to_get_it(
    nextflow: Path, launched, plain
):
    nextflow.unlink()
    result = invoke("-profile", "docker")
    assert result.exit_code == 1
    assert not launched
    assert "Nextflow was not found" in plain(result.output)
    assert "uv tool install plexembler" in plain(result.output)
    # A dry run prints the command all the same, with a warning.
    result = invoke("--dry-run", "-profile", "docker")
    assert result.exit_code == 0
    assert result.stdout == f"nextflow run {REPOSITORY} -r {RELEASE} -profile docker\n"
    assert "WARNING :: Nextflow was not found" in result.stderr


@pytest.mark.usefixtures("released")
def test_a_nextflow_that_cannot_run_is_an_error(
    nextflow: Path, monkeypatch: pytest.MonkeyPatch, plain
):
    def execve(path, argv, env):
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr(os, "execve", execve)
    result = invoke("-profile", "docker")
    assert result.exit_code == 1
    # The path may break over the error panel's lines: its end does not.
    assert "nextflow: [Errno 13] Permission denied" in plain(result.output)


def test_help_and_usage_errors(plain):
    for arguments in (["-h"], ["--help"], []):
        result = invoke(*arguments)
        assert "--pipeline" in result.output
        assert "--dry-run" in result.output
    assert invoke("-h").exit_code == 0
    result = invoke("-profile", "docker", "--pipeline")
    assert result.exit_code == 2
    assert "Option '--pipeline' requires an argument." in plain(result.output)
    assert invoke("--pipeline", "--", "-profile", "docker").exit_code == 2


# Which pipeline, and which revision


def test_split_arguments():
    flags, valued = {"-h", "--dry-run"}, {"--pipeline"}
    assert split_arguments(
        ["-profile", "docker", "--pipeline", "p", "-hub", "--dry-run", "--x=1"],
        flags,
        valued,
    ) == (["--pipeline", "p", "--dry-run"], ["-profile", "docker", "-hub", "--x=1"])
    assert split_arguments(["--pipeline=p", "--", "-h"], flags, valued) == (
        ["--pipeline=p"],
        ["-h"],
    )
    with pytest.raises(ValueError, match="requires an argument"):
        split_arguments(["--pipeline"], flags, valued)


def test_a_release_runs_its_tag(monkeypatch: pytest.MonkeyPatch, tmp_path: Path):
    origin(monkeypatch, None)
    assert default_pipeline() == (REPOSITORY, RELEASE)
    # An archive (a wheel's file), or the folder of a conda package's build,
    # now gone: the release too.
    origin(monkeypatch, {"url": "file:///tmp/dist/x.whl", "archive_info": {}})
    assert default_pipeline() == (REPOSITORY, RELEASE)
    gone = (tmp_path / "work").as_uri()
    origin(monkeypatch, {"url": gone, "dir_info": {}})
    assert default_pipeline() == (REPOSITORY, RELEASE)


@pytest.mark.parametrize("editable", [True, False])
def test_an_install_from_a_folder_runs_that_folder(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, editable: bool
):
    checkout = tmp_path / "my clone"
    checkout.mkdir()
    (checkout / "main.nf").write_text("workflow {}\n", encoding="utf-8")
    origin(monkeypatch, {"url": checkout.as_uri(), "dir_info": {"editable": editable}})
    assert default_pipeline() == (str(checkout), None)
    assert pipeline_to_run(None, ["-profile", "docker"]) == (str(checkout), None)
    # A revision is a repository's.
    assert pipeline_to_run(None, ["-r", "main"]) == (REPOSITORY, None)


@pytest.mark.parametrize(
    ("url", "repository"),
    [
        ("https://github.com/js2264/plexembler", REPOSITORY),
        ("https://github.com/someone/plexembler.git", "someone/plexembler"),
        ("ssh://git@github.com/js2264/plexembler.git", REPOSITORY),
        (
            "https://gitlab.com/someone/plexembler.git",
            "https://gitlab.com/someone/plexembler.git",
        ),
    ],
)
def test_an_install_from_git_runs_its_commit(
    monkeypatch: pytest.MonkeyPatch, url: str, repository: str
):
    commit = "7921be1537eac1e97bc40179a57f0349c2aabfcd"
    vcs_info = {"vcs": "git", "requested_revision": "main", "commit_id": commit}
    origin(monkeypatch, {"url": url, "vcs_info": vcs_info})
    assert default_pipeline() == (repository, commit)
    assert pipeline_to_run(None, []) == (repository, commit)
    assert pipeline_to_run(None, ["-r", "dev"]) == (repository, None)


def test_installed_from_reads_pep_610(monkeypatch: pytest.MonkeyPatch):
    class Distribution:
        def __init__(self, text):
            self.text = text

        def read_text(self, name):
            assert name == "direct_url.json"
            return self.text

    record = {"url": "file:///src/plexembler", "dir_info": {"editable": True}}
    for text, expected in [
        (json.dumps(record), record),
        (None, None),
        ("", None),
        ("{not json", None),
        ("[]", None),
    ]:
        monkeypatch.setattr(
            run, "distribution", lambda _, text=text: Distribution(text)
        )
        assert installed_from() == expected


def test_this_installation_is_found():
    # Whatever installed the tests' plexembler, the answer is one of the three.
    pipeline, revision = default_pipeline()
    assert (Path(pipeline) / "main.nf").is_file() or revision is not None


# Nextflow, and its environment


def test_nextflow_beside_python_comes_first(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    scripts, elsewhere = tmp_path / "env" / "bin", tmp_path / "usr" / "bin"
    for folder in (scripts, elsewhere):
        folder.mkdir(parents=True)
        executable = folder / "nextflow"
        executable.write_text("#!/bin/sh\n", encoding="utf-8")
        executable.chmod(0o755)
    monkeypatch.setattr(sysconfig, "get_path", lambda name: str(scripts))
    monkeypatch.setenv("PATH", str(elsewhere))
    assert find_nextflow() == str(scripts / "nextflow")
    # Its environment's executables come first: its Java is there.
    env = environment(str(scripts / "nextflow"))
    assert env["PATH"].split(os.pathsep) == [str(scripts), str(elsewhere)]
    # The Nextflow of PATH runs in the environment as it is.
    assert environment(str(elsewhere / "nextflow"))["PATH"] == str(elsewhere)
    (scripts / "nextflow").unlink()
    assert find_nextflow() == str(elsewhere / "nextflow")
    monkeypatch.setenv("PATH", "")
    assert find_nextflow() is None


@pytest.mark.usefixtures("released")
def test_nextflow_runs_with_its_environment(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, launched
):
    scripts = tmp_path / "bin"
    scripts.mkdir()
    (scripts / "nextflow").write_text("#!/bin/sh\n", encoding="utf-8")
    (scripts / "nextflow").chmod(0o755)
    monkeypatch.setattr(sysconfig, "get_path", lambda name: str(scripts))
    monkeypatch.setenv("PATH", "/usr/bin")
    monkeypatch.setenv("NXF_VER", "26.04.0")
    invoke("-profile", "docker")
    _, _, env = launched[0]
    assert env["PATH"] == os.pathsep.join([str(scripts), "/usr/bin"])
    assert env["NXF_VER"] == "26.04.0"


# What plexembler run and the pipeline share


def _camel(name: str) -> str:
    head, *rest = name.split("-")
    return head + "".join(word.title() for word in rest)


def test_no_option_of_plexembler_run_hides_a_parameter_of_the_pipeline():
    schema_path = ROOT / "nextflow_schema.json"
    if not schema_path.exists():
        pytest.skip("the pipeline is not in the sdist")
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    parameters = set(schema.get("properties", {}))
    for group in schema.get("$defs", {}).values():
        parameters |= set(group.get("properties", {}))
    assert {"input", "outdir", "help"} <= parameters
    command = cli.commands["run"]
    names = {
        name.removeprefix("--")
        for param in command.params
        if isinstance(param, click.Option)
        for name in [*param.opts, *param.secondary_opts]
        if name.startswith("--")
    }
    # The pipeline's --help is reached after --: plexembler run -- --help.
    names.discard("help")
    # Nextflow reads --dry-run as params.'dry-run' and params.dryRun.
    hidden = (names | {_camel(name) for name in names}) & parameters
    assert not hidden, f"plexembler run's options hide the pipeline's {hidden}"


# Nextflow, which comes with the package: PyPI's nextflow, Nextflow's launcher


def _numbers(version: str) -> tuple[int, ...]:
    """A version's numbers, as PEP 440 compares them: 26.04.0 is 26.4."""
    numbers = [int(part) for part in version.split(".")]
    while numbers and numbers[-1] == 0:
        numbers.pop()
    return tuple(numbers)


def test_nextflow_comes_with_the_package_at_the_pipelines_floor():
    floors = {}
    for requirement in requires("plexembler") or []:
        found = re.fullmatch(r"([A-Za-z0-9._-]+)>=([0-9.]+)", requirement)
        if found:
            floors[found[1]] = found[2]
    assert "nextflow" in floors, "PyPI's nextflow is a dependency of plexembler"
    config = ROOT / "nextflow.config"
    if not config.exists():
        pytest.skip("the pipeline is not in the sdist")
    manifest = re.search(
        r"nextflowVersion\s*=\s*'!?>=(.+?)'", config.read_text(encoding="utf-8")
    )
    assert manifest is not None
    assert _numbers(floors["nextflow"]) == _numbers(manifest[1])


def test_this_installation_has_nextflow_beside_python():
    # Installed with its dependencies, plexembler has Nextflow's launcher
    # beside its Python, which it runs before any Nextflow on PATH.
    beside = shutil.which("nextflow", path=sysconfig.get_path("scripts"))
    assert beside is not None
    assert find_nextflow() == beside
