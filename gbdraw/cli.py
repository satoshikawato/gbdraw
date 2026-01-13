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
import webbrowser
import http.server
import socketserver
import threading
import socket
from typing import NoReturn
from importlib import resources
from functools import partial

from .circular import circular_main
from .linear import linear_main
from .version import __version__

def print_version() -> None:
    print(f"gbdraw version {__version__}")


def find_free_port():
    """Find a free port on localhost"""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('127.0.0.1', 0)) 
        return s.getsockname()[1]

def start_local_server(directory: str):
    """
    Starts a local HTTP server to serve the gbdraw web UI.
    """

    port = find_free_port()
    
    handler = partial(http.server.SimpleHTTPRequestHandler, directory=directory)
    
    try:
        httpd = socketserver.TCPServer(("127.0.0.1", port), handler)
    except OSError as e:
        print(f"Error starting server: {e}", file=sys.stderr)
        return

    thread = threading.Thread(target=httpd.serve_forever)
    thread.daemon = True
    thread.start()

    url = f"http://localhost:{port}"
    print(f"✅ Local server started at: {url}")
    print(f"📂 Serving UI from: {directory}")
    print("Press Ctrl+C to stop the server.")
    
    webbrowser.open(url)

    try:
        while True:
            pass
    except KeyboardInterrupt:
        print("\nStopping server...")
        httpd.shutdown()
        sys.exit(0)

def print_help_message() -> NoReturn:
    print(f"gbdraw v. {__version__}: A diagram generator for small genomes")
    print("")
    print("Usage:")
    print("  gbdraw <subcommand> [options]")
    print("")
    print("Subcommands:")
    print("  circular  Generate a circular genome diagram")
    print("  linear    Generate a linear genome diagram")
    print("  gui       Launch the local web UI (no Streamlit required)")
    print("")
    print("For each subcommand, you can get additional help by running:")
    print("  gbdraw <subcommand> --help")
    print("")
    print("Examples:")
    print("  gbdraw circular --gbk input.gb")
    print("  gbdraw circular --gff input.gff --fasta input.fna")
    print("  gbdraw linear --gbk input.gb")
    print("  gbdraw linear --gff input.gff --fasta input.fna")
    print("  gbdraw linear --gbk input1.gb input2.gb input3.gb -b input1_input2.blast.outfmt7.txt input2_input3.blast.outfmt7.txt")
    print("  gbdraw linear --gff input1.gff input2.gff input3.gff --fasta input1.fna input2.fna input3.fna -b input1_input2.blast.outfmt7.txt input2_input3.blast.outfmt7.txt")
    print("  gbdraw gui")
    print("")
    print("Options (examples):")
    print("  --gbk                Input GenBank file(s)")
    print("  --gff                Input GFF3 file(s) (requires --fasta; mutually exclusive with --gbk)")
    print("  --fasta              Input FASTA file(s) (required with --gff; mutually exclusive with --gbk)")
    print("  -o, --output         Output file prefix (optional)")
    print("  -b, --blast          BLAST result file in tab-separated format (-outfmt 6 or 7) (optional; implemented for linear mode only)")

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
            try:
                # Python 3.9+
                web_dir = resources.files('gbdraw').joinpath('web')
            except AttributeError:
                import os
                web_dir = os.path.join(os.path.dirname(__file__), 'web')

            if not str(web_dir).endswith("index.html") and not str(web_dir).endswith("web"):
                pass
            
            print("Launching gbdraw GUI (Local Mode)...")
            start_local_server(str(web_dir))
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
