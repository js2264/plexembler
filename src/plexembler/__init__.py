"""plexembler: genome references for every species of a pool.

The Python side of the pipeline: the commands its steps run on their data
(``plexembler <command>``), each tested in-process, and the report of a run.
"""

from importlib.metadata import version

__version__ = version("plexembler")

__all__ = ["__version__"]
