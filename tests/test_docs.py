"""The command line's reference, docs/cli.md, follows the command line.

docs/cli.md is written by hand: a section per command, headed
``### `command` { #command }``, or ``#### `group command` `` for a command of
a group, with a table whose first column lists the command's options, short
form first (``-o, --outdir``; ``--balance/--no-balance``). These tests fail
when the page misses a command or an option, or documents one that the command
line does not have.

- A command without a section of its own is documented in its group's; a group
  needs no section when each of its commands has one, and it has no option.
- Commands that share their options may share a section, whose heading names
  each of them: ``#### `normalize normcc` and `normalize hiczin` ``.
- A section runs to the next heading of its level or above: the tables of its
  sub-sections (``#### Examples``) document its commands too.
- A table outside every command's section, such as the page's conventions,
  documents its options for every command that has them.
- Code blocks are not read: a comment in an example is not a heading.
- -h/--help and --version are left out; a hidden option, or a hidden command
  (an alias), may be documented, and need not be.

chromokit's template (CONVENTIONS.md § 8) writes this file, and `copier
update` brings its fixes.
"""

import re
from pathlib import Path

import click
import pytest

from plexembler.cli import cli

CLI_MD = Path(__file__).resolve().parents[1] / "docs" / "cli.md"

#: A heading of the page, its level and its text.
_HEADING = re.compile(r"^(#{1,6}) +(.*)$", re.MULTILINE)
#: A fenced code block, which holds neither headings nor tables.
_FENCE = re.compile(r"^[ \t]*(`{3,}|~{3,}).*?^[ \t]*\1[ \t]*$", re.M | re.S)
#: A command's name in a heading: ### `pipeline` { #pipeline }, #### `create pairs`.
_COMMAND = re.compile(r"`([a-z0-9][a-z0-9_ -]*)`")
#: What joins the commands of a heading that names several: `a` and `b`, `c`.
_JOINED = re.compile(r"[\s,/&]*(?:(?:and|or)[\s,/&]+)*")
#: The first cell of a table's row, and an option in it: | `-o, --outdir` | …
_FIRST_CELL = re.compile(r"^\|([^|\n]*)\|", re.MULTILINE)
_OPTION = re.compile(r"`(-[^`]+)`")
#: Options every command has, which the page leaves out.
_LEFT_OUT = {"-h", "--help", "-V", "--version"}


def _options(text):
    """The options of the tables of text: the first cell's, in backticks."""
    found = set()
    for cell in _FIRST_CELL.findall(text):
        for span in _OPTION.findall(cell):
            for name in re.split(r"[,/]", span):
                if name.strip():
                    found.add(name.split()[0].split("=")[0])
    return found - _LEFT_OUT


def _heading_commands(heading):
    """The commands a heading names: ``### `pipeline` { #pipeline }``; or, for a
    section that several share, each of them: ``#### `normcc` and `hiczin` ``."""
    title = heading.split("{")[0].strip()
    names = _COMMAND.findall(title)
    if not names or not title.startswith("`"):
        return []
    if _JOINED.fullmatch(_COMMAND.sub("", title)):
        return names
    return names[:1]


def documented():
    """The page's command sections, each the commands its heading names and the
    options its tables document, its sub-sections' included; and the options of
    the tables outside every command's section, which document every command."""
    text = _FENCE.sub("", CLI_MD.read_text(encoding="utf-8"))
    headings = list(_HEADING.finditer(text))
    sections = []
    shared = _options(text[: headings[0].start()] if headings else text)
    # The command sections that the heading is in, the innermost last.
    within = []
    for heading, following in zip(headings, [*headings[1:], None], strict=True):
        end = following.start() if following else len(text)
        found = _options(text[heading.end() : end])
        level = len(heading.group(1))
        names = _heading_commands(heading.group(2)) if level in (3, 4) else []
        while within and within[-1][0] >= level:
            within.pop()
        if names:
            sections.append((names, found))
            within.append((level, found))
        elif within:
            within[-1][1].update(found)
        else:
            shared |= found
    return sections, shared


def _by_name(sections):
    """The options that the sections of each command document, by its name."""
    found = {}
    for names, documented_here in sections:
        for name in names:
            found.setdefault(name, set()).update(documented_here)
    return found


def commands():
    """Each command of the command line, by its name: a command of a group by
    the group's name and its own ("create pairs")."""
    found = {}

    def walk(command, path):
        found[" ".join(path)] = command
        if isinstance(command, click.Group):
            context = click.Context(command)
            for name in command.list_commands(context):
                walk(command.get_command(context, name), [*path, name])

    context = click.Context(cli)
    for name in cli.list_commands(context):
        walk(cli.get_command(context, name), [name])
    return found


def options(path):
    """The options of a command, and its hidden ones apart."""
    shown, hidden = set(), set()
    for param in commands()[path].params:
        if isinstance(param, click.Option):
            names = param.opts + param.secondary_opts
            (hidden if param.hidden else shown).update(names)
    return shown - _LEFT_OUT, hidden - _LEFT_OUT


def _groups(path):
    """Whether a command is a group."""
    return isinstance(commands()[path], click.Group)


def section_of(path, sections):
    """The section that documents a command: its own, else its group's."""
    words = path.split()
    while words:
        if " ".join(words) in sections:
            return " ".join(words)
        words.pop()
    return None


def _needs_a_section(path):
    """A command needs a section, but a group without options of its own."""
    return not (_groups(path) and not options(path)[0])


def _hidden(path):
    """Whether a command is hidden, itself or as one of a hidden group (an alias)."""
    words = path.split()
    paths = [" ".join(words[:end]) for end in range(1, len(words) + 1)]
    return any(getattr(commands()[each], "hidden", False) for each in paths)


def test_every_command_has_a_section():
    named = _by_name(documented()[0])
    shown = [path for path in commands() if not _hidden(path)]
    undocumented = [path for path in shown if section_of(path, named) is None]
    missing = [path for path in undocumented if _needs_a_section(path)]
    assert not missing, f"commands that docs/cli.md misses: {missing}"


def test_every_section_is_a_command():
    named = _by_name(documented()[0])
    unknown = sorted(set(named) - set(commands()))
    assert not unknown, f"docs/cli.md documents commands that do not exist: {unknown}"


@pytest.mark.parametrize("path", sorted(p for p in commands() if not _hidden(p)))
def test_every_option_is_documented(path):
    sections, shared = documented()
    named = _by_name(sections)
    shown, _ = options(path)
    documented_here = named.get(section_of(path, named), set()) | shared
    missing = sorted(shown - documented_here)
    assert not missing, f"options of `{path}` that docs/cli.md misses: {missing}"


def test_every_documented_option_exists():
    sections, shared = documented()
    named = _by_name(sections)
    known = {path: set().union(*options(path)) for path in commands()}
    unknown = []
    for names, documented_here in sections:
        here = [known[path] for path in known if section_of(path, named) in names]
        extra = sorted(documented_here - set().union(*here))
        unknown += [f"`{' / '.join(names)}` {name}" for name in extra]
    every = set().union(*known.values())
    unknown += [f"(every command) {name}" for name in sorted(shared - every)]
    assert not unknown, f"docs/cli.md documents options that do not exist: {unknown}"
