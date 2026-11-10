"""Run the ``gbdraw`` command in the test's own process.

A test whose subject is a figure, a Session or a message, not the process, runs the command here instead of starting
a Python interpreter for it: the start alone (``import gbdraw.cli``) costs about a second. The call sees what a
fresh ``python -m gbdraw.cli`` sees where its output depends on it (the argument list and ``__main__.__spec__``,
from which argparse names the program, the working folder, the warning filters a fresh interpreter holds once
``gbdraw.cli`` is imported (its defaults plus those that numpy and the other imports add) with warnings printed to
standard error, a root logger without handlers so the command installs its own on standard error, and the config
file not yet read, so its once-per-process log line appears), returns the exit status a process would, and leaves
the test's process as it found it.

Tests about the process itself (the entry point, imports in a fresh interpreter, the installed package, equality
between the CLI and Python surfaces) keep their subprocesses.
"""

from __future__ import annotations

import contextlib
import importlib
import inspect
import io
import json
import logging
import os
import subprocess
import sys
import traceback
import warnings
from collections.abc import Iterator, Sequence
from functools import cache
from pathlib import Path
from typing import NamedTuple

from gbdraw import cli
from gbdraw.config import toml as config_toml
from gbdraw.core import fonts

# ``python -m gbdraw.cli`` puts the module's file in ``sys.argv[0]`` and the module's spec in ``__main__.__spec__``.
# argparse names the program after the file before Python 3.14 (``cli.py``) and after the spec from 3.14 on
# (``python -m gbdraw.cli``), so a run here sets both.
_PROGRAM = cli.__file__
_PROGRAM_SPEC = cli.__spec__

# One warning filter as text: action, message pattern, category ("module:qualname"), module pattern, line.
WarningFilter = tuple[str, str, str, str, int]


def describe_filters(filters: Sequence[tuple[object, ...]]) -> list[WarningFilter]:
    """``warnings.filters`` as ``WarningFilter`` rows (a fresh interpreter runs this function's source too)."""

    return [
        (action, getattr(message, "pattern", message or ""), f"{category.__module__}:{category.__qualname__}",
         getattr(module, "pattern", module or ""), lineno)
        for action, message, category, module, lineno in filters
    ]


@cache
def fresh_process_warning_filters() -> tuple[WarningFilter, ...]:
    """The warning filters of a fresh interpreter once ``gbdraw.cli`` is imported, read once per test process.

    The test process imported numpy and the rest long ago, under pytest's own filters, so it cannot rebuild the
    filters those imports add (numpy's ``RankWarning`` ``always`` among them) except by asking a fresh one.
    """

    script = (
        "from __future__ import annotations\n" + inspect.getsource(describe_filters)
        + "import json, warnings, gbdraw.cli\nprint(json.dumps(describe_filters(warnings.filters)))"
    )
    printed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, check=True).stdout
    return tuple(tuple(row) for row in json.loads(printed))


def _install_warning_filters(filters: Sequence[WarningFilter]) -> None:
    warnings.resetwarnings()
    for action, message, category, module, lineno in filters:
        category_module, qualname = category.split(":")
        found: object = importlib.import_module(category_module)
        for name in qualname.split("."):
            found = getattr(found, name)
        assert isinstance(found, type) and issubclass(found, Warning), category
        warnings.filterwarnings(action, message, found, module, lineno, append=True)


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
    main = sys.modules["__main__"]
    saved_main_spec = main.__spec__
    saved_cwd = os.getcwd()
    saved_fonts = {registered.face.path for registered in fonts.registered_fonts()}
    for handler in saved_handlers:
        root.removeHandler(handler)
    config_toml._config_text.cache_clear()
    try:
        with warnings.catch_warnings():
            _install_warning_filters(fresh_process_warning_filters())
            warnings.showwarning = _print_warning
            sys.argv = argv
            main.__spec__ = _PROGRAM_SPEC
            os.chdir(cwd)
            yield
    finally:
        os.chdir(saved_cwd)
        sys.argv = saved_argv
        main.__spec__ = saved_main_spec
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
