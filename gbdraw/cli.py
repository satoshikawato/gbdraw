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
import subprocess
from typing import NoReturn
from .circular import circular_main
from .linear import linear_main
from importlib import resources
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
    print("  gui       Launch the graphical user interface (requires Streamlit)")
    print("")
    print("For each subcommand, you can get additional help by running:")
    print("  gbdraw <subcommand> --help")
    print("")
    print("Examples:")
    print("  gbdraw circular --gbk input.gb")
    print("  gbdraw linear --gbk input.gb")
    print("  gbdraw gui")
    print("")
    print("Options (examples):")
    print("  --gbk                Input GenBank file(s)")
    print("  --gff               Input GFF# file(s) (rquires --fasta; mutually exclusive with --gbk)")
    print("  --fasta             Input FASTA file(s) (required with --gff; mutually exclusive with --gbk)")
    print("  -o, --output         Output file prefix (optional)")
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
    elif command == "gui":
        # Check if streamlit is installed
        try:
            import streamlit
        except ImportError:
            print("Error: Streamlit is not installed.", file=sys.stderr)
            print("Please install it to use the GUI:", file=sys.stderr)
            print("Examples:", file=sys.stderr)
            print("mamba install -c conda-forge streamlit", file=sys.stderr)
            print("micromamba install -c conda-forge streamlit", file=sys.stderr)
            sys.exit(1)

        # Get the path to the app.py file within the package
        app_path = resources.files('gbdraw').joinpath('app.py')
        
        # Launch the Streamlit app
        print("Launching gbdraw GUI...")
        print("If the GUI does not open automatically, please open your web browser and navigate to one of the URLs provided below.")
        print("When you are finished using the GUI, return to this terminal and press Ctrl+C to stop the server.")
        subprocess.run([sys.executable, "-m", "streamlit", "run", str(app_path)])
        
    else:
        print("Oops! It seems like you entered an invalid command.")
        print("Please use 'circular' or 'linear' followed by the respective options.")
        print("For example:")
        print("  gbdraw circular --gbk input.gb")
        print("  gbdraw linear --gbk input.gb")
        print("  gbdraw gui")        
        sys.exit(1)


if __name__ == "__main__":
    main()
