#!/usr/bin/env python
# coding: utf-8

"""
Command-line Interface for Genome Diagram Generation (gbdraw)

This script provides the command-line interface for `gbdraw`, a tool for generating 
genome diagrams in circular or linear formats from GenBank/EMBL/DDBJ-format annotated 
genomes. It supports various customization options, including color schemes, feature 
inclusions, and skewness calculations.

Functions:
    main() - Parses command-line arguments and redirects to appropriate functions for diagram generation.
    circular_main(args) - Handles the generation of circular genome diagrams.
    linear_main(args) - Handles the generation of linear genome diagrams.
"""


import sys
from typing import NoReturn
from .circular import circular_main
from .linear import linear_main
from .version import __version__


def print_version() -> None:
    print(f"gbdraw version {__version__}")


def print_help_message() -> NoReturn:
    print(f"gbdraw v. {__version__}: A diagram generator for small genomes")
    print("")
    print("Usage:")
    print("  gbdraw <subcommand> [options]")
    print("")
    print("Subcommands:")
    print("  circular  Generate a circular genome diagram")
    print("  linear    Generate a linear genome diagram")
    print("")
    print("For each subcommand, you can get additional help by running:")
    print("  gbdraw <subcommand> --help")
    print("")
    print("Examples:")
    print("  gbdraw circular -i input.gb")
    print("  gbdraw linear -i input.gb")
    print("")
    print("Options (examples):")
    print("  -i, --input          Input GenBank file(s) (required)")
    print("  -o, --output         Output file prefix (optional)")
    print("  -t, --table          Color table file (optional)")
    print("  -b, --blast          BLAST result file in tab-separated format (-outfmt 6 or 7) (optional; currently implemented for linear mode only)")
    print("")
    print("Additional Information:")
    print("  - For full documentation, visit: https://github.com/satoshikawato/gbdraw/")
    print("  - For issues and source code, visit the GitHub repository: https://github.com/satoshikawato/gbdraw/")
    print("  - For support, contact: kawato[at]kaiyodai.ac.jp")
    sys.exit(0)


def main() -> None:
    """
    Parses command-line arguments and directs to the appropriate genome diagram generation function.

    This function serves as the entry point for the command-line interface. It parses the arguments
    provided by the user and determines whether to generate a circular or linear genome diagram.
    Based on the user's choice, it calls either `circular_main` or `linear_main` functions with the
    parsed arguments.

    Args:
        None - Command-line arguments are obtained from sys.argv by argparse.

    This function does not return any values. Instead, it initiates the genome diagram generation process.
    """
    if len(sys.argv) <= 1 or sys.argv[1] in ['-h', '--help']:
        print_help_message()
        sys.exit(0)
    elif '--version' in sys.argv:
        print_version()
        sys.exit(0)

    command: str = sys.argv[1].lower()
    args: list[str] = sys.argv[2:]

    if command == "circular":
        circular_main(args)
    elif command == "linear":
        linear_main(args)
    else:
        print("Oops! It seems like you entered an invalid command.")
        print("Please use 'circular' or 'linear' followed by the respective options.")
        print("For example:")
        print("  gbdraw circular -i input.gb")
        print("  gbdraw linear -i input.gb")
        sys.exit(1)


if __name__ == "__main__":
    main()
