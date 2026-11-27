"""What the tests share: a command's output read as plain text, and sequences."""

import random
import re

import pytest

#: A terminal's colour codes, which rich-click writes where it sees a terminal (CI).
ANSI = re.compile(r"\x1b\[[0-9;]*m")


def plain_text(output: str) -> str:
    """The output without rich-click's panel: its colours, borders and wrapping.

    rich-click wraps an error in a panel as wide as the terminal; a message that
    names a long path (macOS' temporary folders) breaks over its lines.
    """
    text = ANSI.sub("", output)
    return " ".join(text.translate(str.maketrans("", "", "│╭╮╰╯─")).split())


@pytest.fixture
def plain():
    """``plain(result.output)``: a command's output as one line of text."""
    return plain_text


def random_sequence(length: int, seed: int, gc: float = 0.4) -> str:
    """A random sequence of ``length`` bases, the same for the same seed."""
    rng = random.Random(seed)
    at, cg = (1 - gc) / 2, gc / 2
    return "".join(rng.choices("ACGT", weights=[at, cg, cg, at], k=length))


@pytest.fixture
def sequence():
    """``sequence(length, seed, gc=0.4)``: a random sequence, the same for the
    same seed."""
    return random_sequence
