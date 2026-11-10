"""Run the ``gbdraw`` command in the test's own process.

A test whose subject is a figure, a Session or a message, not the process, runs the command here instead of
starting a Python interpreter for it: the start alone (``import gbdraw.cli``) costs about a second. The call sees
what a fresh ``python -m gbdraw.cli`` sees where its output depends on it (the argument list and its program name,
the working folder, the interpreter's warning filters with warnings printed to standard error, a root logger
without handlers so the command installs its own on standard error, and the config file not yet read, so its
once-per-process log line appears), returns the exit status a process would, and
leaves the test's process as it found it.

Tests about the process itself (the entry point, imports in a fresh interpreter, the installed package, equality
between the CLI and Python surfaces) keep their subprocesses.
"""

from __future__ import annotations

import contextlib
import io
import logging
import os
import sys
import traceback
import warnings
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import NamedTuple

from gbdraw import cli
from gbdraw.config import toml as config_toml
from gbdraw.core import fonts

# ``python -m gbdraw.cli`` puts the module's file in ``sys.argv[0]``; argparse names the program after it.
_PROGRAM = cli.__file__

# The warning filters of an interpreter started without ``-W``, ``-X dev`` or ``PYTHONWARNINGS``.
_INTERPRETER_WARNING_FILTERS: tuple[tuple[str, type[Warning], str], ...] = (
    ("default", DeprecationWarning, "__main__"),
    ("ignore", DeprecationWarning, ""),
    ("ignore", PendingDeprecationWarning, ""),
    ("ignore", ImportWarning, ""),
    ("ignore", ResourceWarning, ""),
)


class CliResult(NamedTuple):
    """What a finished ``gbdraw`` command returned and printed."""

    returncode: int
    stdout: str
    stderr: str


def _print_warning(
    message: Warning | str,
    category: type[Warning],
    filename: str,
    lineno: int,
    file: object = None,
    line: str | None = None,
) -> None:
    """Print a warning to standard error as the interpreter does when nothing records it."""

    sys.stderr.write(warnings.formatwarning(message, category, filename, lineno, line))


@contextlib.contextmanager
def _fresh_process_state(cwd: Path, argv: list[str]) -> Iterator[None]:
    """Give the command a fresh process's argv, folder, warnings and logging; restore the test's afterwards."""

    root = logging.getLogger()
    saved_handlers = list(root.handlers)
    saved_level = root.level
    saved_argv = sys.argv
    saved_cwd = os.getcwd()
    saved_fonts = {registered.face.path for registered in fonts.registered_fonts()}
    for handler in saved_handlers:
        root.removeHandler(handler)
    config_toml._config_text.cache_clear()
    try:
        with warnings.catch_warnings():
            warnings.resetwarnings()
            for action, category, module in reversed(_INTERPRETER_WARNING_FILTERS):
                warnings.filterwarnings(action, category=category, module=module)
            warnings.showwarning = _print_warning
            sys.argv = argv
            os.chdir(cwd)
            yield
    finally:
        os.chdir(saved_cwd)
        sys.argv = saved_argv
        for handler in list(root.handlers):
            root.removeHandler(handler)
        for handler in saved_handlers:
            root.addHandler(handler)
        root.setLevel(saved_level)
        for registered in fonts.registered_fonts():
            if registered.face.path not in saved_fonts:
                fonts.unregister_font_file(registered.face.path)


def run_gbdraw(args: Sequence[str], *, cwd: Path) -> CliResult:
    """Run ``gbdraw <args>`` in ``cwd`` in this process and return its exit status, stdout and stderr.

    The exit status follows the interpreter: ``SystemExit(None)`` is 0, an integer is itself, any other value is
    printed to stderr and gives 1, and an exception that escapes the command prints its traceback and gives 1.
    """

    stdout = io.StringIO()
    stderr = io.StringIO()
    returncode = 0
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        with _fresh_process_state(Path(cwd), [_PROGRAM, *map(str, args)]):
            try:
                cli.main()
            except SystemExit as exit_:
                if exit_.code is None or isinstance(exit_.code, int):
                    returncode = exit_.code or 0
                else:
                    print(exit_.code, file=sys.stderr)
                    returncode = 1
            except Exception:
                traceback.print_exc(file=sys.stderr)
                returncode = 1
    return CliResult(returncode, stdout.getvalue(), stderr.getvalue())
