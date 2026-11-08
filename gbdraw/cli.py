#!/usr/bin/env python
# coding: utf-8

"""
The gbdraw command: it dispatches to the circular, linear, render, setup-losat and gui subcommands.

This script provides the command-line interface for `gbdraw`, a tool for generating
genome diagrams in circular or linear formats from GenBank/EMBL/DDBJ-format annotated
genomes. It supports various customization options, including color schemes, feature
inclusions, and skewness calculations.

Functions:
    main() - Parses command-line arguments and redirects to appropriate functions for diagram generation.
    circular_main(args) - Handles the generation of circular genome diagrams.
    linear_main(args) - Handles the generation of linear genome diagrams.
"""


import argparse
import hashlib
import os
import socket
import sys
import warnings
import webbrowser
import http.server
import socketserver
from http import HTTPStatus
from collections.abc import Callable
from typing import TYPE_CHECKING, BinaryIO, NoReturn, TextIO
from importlib import resources
from functools import partial

from .api.request_render import _without_derived_losat_entries
from .circular import circular_main
from .linear import linear_main
from .exceptions import LabelLeadersCrossWarning, LabelsMayNotFitWarning, GbdrawError
from .render.output_paths import CLI_OVERWRITE_OPTION, OutputExistsError
from .version import __version_display__

if TYPE_CHECKING:
    if sys.version_info >= (3, 11):
        from importlib.resources.abc import Traversable
    else:
        from importlib.abc import Traversable


class IsolatedSimpleHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    """Serve the browser UI with the headers required for SharedArrayBuffer.

    The browser revalidates every file, and a cached file is reused only when its
    ETag matches exactly: the ETag names the file's path, size and modification
    time and the gbdraw version, so another installation or an older gbdraw at the
    same address never runs from the cache. A modification date alone never
    answers 304: a downgrade or another environment can serve older dates."""

    _etag: str | None = None

    def _file_etag(self) -> str | None:
        path = self.translate_path(self.path)
        if os.path.isdir(path):
            path = os.path.join(path, "index.html")
        try:
            stat = os.stat(path)
        except OSError:
            return None
        if not os.path.isfile(path):
            return None
        identity = f"{os.path.realpath(path)}\0{stat.st_size}\0{stat.st_mtime_ns}\0{__version_display__}"
        return '"' + hashlib.sha256(identity.encode("utf-8", "surrogateescape")).hexdigest()[:32] + '"'

    def send_head(self) -> BinaryIO | None:
        del self.headers["If-Modified-Since"]
        self._etag = self._file_etag()
        matches = {tag.strip() for tag in self.headers.get("If-None-Match", "").split(",")}
        if self._etag is not None and self._etag in matches:
            self.send_response(HTTPStatus.NOT_MODIFIED)
            self.end_headers()
            return None
        return super().send_head()

    def end_headers(self) -> None:
        if self._etag is not None:
            self.send_header("ETag", self._etag)
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        self.send_header("Cross-Origin-Resource-Policy", "same-origin")
        self.send_header("Content-Security-Policy", "frame-ancestors 'none'")
        super().end_headers()

def print_version() -> None:
    print(f"gbdraw version {__version_display__}")


GUI_PORT = 24372
"""The port ``gbdraw gui`` serves the Web app on. One address across runs keeps the
browser's storage for it: the layout choice, and the temporary input copies that
the next run deletes when a run ended without deleting them (PD-OI-121)."""


class _GuiServer(socketserver.TCPServer):
    # A port the previous run just closed can be bound again at once; Windows
    # would let a second server bind a port in use with this option, so not there.
    allow_reuse_address = sys.platform != "win32"


def _port_answers(port: int) -> bool:
    """Whether a server already accepts connections on ``port`` of 127.0.0.1 or ::1,
    either of which ``localhost`` can reach. It also finds a server on every
    address, which the reuse option would let a bind on macOS share."""
    for family, host in ((socket.AF_INET, "127.0.0.1"), (socket.AF_INET6, "::1")):
        try:
            with socket.socket(family, socket.SOCK_STREAM) as probe:
                probe.settimeout(0.5)
                if probe.connect_ex((host, port)) == 0:
                    return True
        except OSError:
            continue
    return False


def bind_gui_server(handler: Callable[..., http.server.BaseHTTPRequestHandler]) -> tuple[socketserver.TCPServer, bool]:
    """The server on ``GUI_PORT`` of 127.0.0.1, or, when a server answers on that port
    or it cannot be bound, on a free port; the flag says whether it fell back."""
    if not _port_answers(GUI_PORT):
        try:
            return _GuiServer(("127.0.0.1", GUI_PORT), handler), False
        except OSError:
            pass
    return _GuiServer(("127.0.0.1", 0), handler), True


def start_local_server(directory: str):
    """
    Starts a local HTTP server to serve the gbdraw web UI.
    """

    handler = partial(IsolatedSimpleHTTPRequestHandler, directory=directory)

    try:
        httpd, fell_back = bind_gui_server(handler)
    except OSError as e:
        print(f"Error starting server: {e}", file=sys.stderr)
        return

    port = httpd.server_address[1]
    url = f"http://localhost:{port}"
    if fell_back:
        print(f"Port {GUI_PORT} is in use, so the UI is served at another address.")
        print("The layout choice made there is not kept, and its temporary input copies are")
        print(f"deleted only when port {port} is used again.")
    print(f"✅ Local server started at: {url}")
    print(f"📂 Serving UI from: {directory}")
    print("Press Ctrl+C to stop the server.")

    webbrowser.open(url)

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server...")
    finally:
        httpd.shutdown()
        httpd.server_close()

def print_help_message() -> NoReturn:
    print(f"gbdraw v. {__version_display__}: A diagram generator for small genomes")
    print("")
    print("Usage:")
    print("  gbdraw <subcommand> [options]")
    print("")
    print("Subcommands:")
    print("  circular  Generate a circular genome diagram")
    print("  linear    Generate a linear genome diagram")
    print("  render    Render the drawings of a saved Session")
    print("  setup-losat  Install the pinned native LOSAT release")
    print("  gui       Launch the local web UI in your browser")
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
    print("  gbdraw linear --gbk input1.gb input2.gb --losat losatp --losatp_mode pairwise")
    print("  gbdraw linear --gff input1.gff input2.gff input3.gff --fasta input1.fna input2.fna input3.fna -b input1_input2.blast.outfmt7.txt input2_input3.blast.outfmt7.txt")
    print("  gbdraw render --session project.gbdraw-session.json")
    print("  gbdraw gui")
    print("")
    print("Options (examples):")
    print("  --gbk                Input GenBank file(s)")
    print("  --gff                Input GFF3 file(s) (requires --fasta; mutually exclusive with --gbk)")
    print("  --fasta              Input FASTA file(s) (with --gff, its sequences; alone, records without features; mutually exclusive with --gbk)")
    print("  -o, --output         Output file prefix (optional)")
    print("  --overwrite          Replace existing output files (default: refuse)")
    print("  -b, --blast          BLAST result file in tab-separated format (-outfmt 6 or 7) (optional; implemented for linear mode only)")
    print("  --losat losatp       Run a LOSATP protein comparison for linear diagrams")
    print("  --losatp_mode {similarity_groups,collinear,pairwise}  LOSATP display (default: similarity_groups)")
    print("  --losat_threads N    Pass thread count to each LOSAT job")

    print("")
    print("Additional Information:")
    print("  - For full documentation, visit: https://github.com/satoshikawato/gbdraw/")
    print("  - For issues and source code, visit the GitHub repository: https://github.com/satoshikawato/gbdraw/")
    print("  - For support, contact: kawato[at]kaiyodai.ac.jp")
    sys.exit(0)


def _cli_error_message(error: GbdrawError) -> str:
    """The error's message, with an existing-output refusal worded for the CLI.

    A wrapper keeps the refusal as its cause (``raise ... from``).
    """

    cause: BaseException | None = error
    seen: set[int] = set()
    while cause is not None and id(cause) not in seen:
        seen.add(id(cause))
        if isinstance(cause, OutputExistsError):
            return cause.message(CLI_OVERWRITE_OPTION)
        cause = cause.__cause__
    return str(error)


def _gui_parser() -> argparse.ArgumentParser:
    """The ``gbdraw gui`` parser: it takes no options, so ``--help`` and unknown arguments stop before the server starts."""

    return argparse.ArgumentParser(
        prog="gbdraw gui",
        description="Serve the gbdraw web UI from this installation on localhost and open it in your browser.",
    )


class _CliWarningPrinter:
    """Print gbdraw's label fit warnings like its errors; pass other warnings on."""

    def __init__(self, show_other: Callable[..., None]) -> None:
        self._show_other = show_other

    def __call__(
        self,
        message: Warning | str,
        category: type[Warning],
        filename: str,
        lineno: int,
        file: TextIO | None = None,
        line: str | None = None,
    ) -> None:
        if issubclass(category, (LabelLeadersCrossWarning, LabelsMayNotFitWarning)):
            print(f"WARNING: {message}", file=sys.stderr)
            return
        self._show_other(message, category, filename, lineno, file, line)


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
    elif sys.argv[1] == '--version':
        print_version()
        sys.exit(0)

    command: str = sys.argv[1].lower()
    args: list[str] = sys.argv[2:]

    # A CLI render never reads the derived LOSAT entries, so it does not build them.
    with warnings.catch_warnings(), _without_derived_losat_entries():
        # Crossing outer leaders draw the figure and print their report (LP-D46).
        warnings.simplefilter("always", LabelLeadersCrossWarning)
        # Labels near the room around the map warn before they are placed (CAPACITY-WARN).
        warnings.simplefilter("always", LabelsMayNotFitWarning)
        warnings.showwarning = _CliWarningPrinter(warnings.showwarning)
        try:
            if command == "circular":
                circular_main(args)
            elif command == "linear":
                linear_main(args)
            elif command == "render":
                from .cli_utils.session import render_main
                render_main(args)
            elif command == "setup-losat":
                from .losat_setup import setup_main
                setup_main(args)
            elif command == "gui":
                _gui_parser().parse_args(args)
                web_dir: Traversable | str
                try:
                    # Python 3.9+
                    web_dir = resources.files('gbdraw').joinpath('web')
                except AttributeError:
                    import os
                    web_dir = os.path.join(os.path.dirname(__file__), 'web')

                print("Launching gbdraw GUI (Local Mode)...")
                start_local_server(str(web_dir))
            else:
                print(
                    f"Unknown command {sys.argv[1]!r}. Choose circular, linear, render, setup-losat, or gui; "
                    "run gbdraw --help for usage.",
                    file=sys.stderr,
                )
                sys.exit(1)
        except GbdrawError as exc:
            print(f"ERROR: {_cli_error_message(exc)}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
