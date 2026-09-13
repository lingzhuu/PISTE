"""Command-line entry point for PathwayLGM."""

import argparse
from collections.abc import Sequence

from pathwaylgm import __version__
from pathwaylgm.config import RuntimeConfig, configure_logging


def create_parser() -> argparse.ArgumentParser:
    """Create the top-level parser without starting an analysis."""
    parser = argparse.ArgumentParser(
        prog="pathwaylgm",
        description="Pathway-level modeling of raw or normalized spatial expression data.",
        epilog="Package initialization only. Analysis commands will be added incrementally.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", "-v", action="version", version=f"PathwayLGM {__version__}")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the package entry point.

    Parameters
    ----------
    argv : sequence of str, optional
        Command-line arguments. Defaults to the process arguments.

    Returns
    -------
    int
        Zero after displaying usage. Argument parsing handles help, version
        and invalid arguments using standard argparse exit codes.
    """
    parser = create_parser()
    parser.parse_args(argv)
    configure_logging(RuntimeConfig())
    parser.print_help()
    return 0
