"""Console entry point of the installed ``gbdraw`` command.

It lives outside the ``gbdraw`` package because ``import gbdraw`` imports NumPy, and OpenBLAS reads
``OPENBLAS_NUM_THREADS`` when NumPy loads. Setting the variable here, before the package is imported,
keeps the one-shot command from starting a BLAS thread pool it does not use. A value the user already
set is kept, and the Python API (``import gbdraw``) is not affected.
"""

from __future__ import annotations

import os


def main() -> None:
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    from gbdraw.cli import main as cli_main

    cli_main()
