"""The command line as a whole: its version, its help and its commands."""

from click.testing import CliRunner

from plexembler import __version__
from plexembler.cli import cli


def test_version():
    result = CliRunner().invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_help_lists_the_commands():
    result = CliRunner().invoke(cli, ["-h"])
    assert result.exit_code == 0
    assert "rename" in result.output


def test_an_unknown_command_is_a_usage_error():
    result = CliRunner().invoke(cli, ["assemble"])
    assert result.exit_code == 2
