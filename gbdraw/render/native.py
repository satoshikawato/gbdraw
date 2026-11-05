"""The optional Rust drawing engine, ``gbdraw_core``.

Native Python imports it as an extension (``rust/python``); the Web Worker
registers its Wasm module under the same name (``workers/gbdraw-core.js``). A
group the engine does not draw is drawn in Python, and so is every group when
the engine is missing or ``GBDRAW_NATIVE=0``.
"""

from __future__ import annotations

import importlib
import os
from types import ModuleType


def native_engine() -> ModuleType | None:
    """``gbdraw_core``, or ``None`` when Python draws everything."""

    if os.environ.get("GBDRAW_NATIVE") == "0":
        return None
    try:
        return importlib.import_module("gbdraw_core")
    except ImportError:
        return None


__all__ = ["native_engine"]
