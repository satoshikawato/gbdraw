"""The optional Rust drawing engine, ``gbdraw_core``.

Native Python imports it as an extension (``rust/python``); the Web Worker
registers its Wasm module under the same name (``workers/gbdraw-core.js``). A
group the engine does not draw is drawn in Python, and so is every group when
the engine is missing, its contract is not ``ENGINE_CONTRACT``, or
``GBDRAW_NATIVE=0``.
"""

from __future__ import annotations

import functools
import importlib
import logging
import os
from types import ModuleType

logger = logging.getLogger(__name__)

# The engine's CONTRACT_VERSION (rust/engine/src/lib.rs) this package draws with.
ENGINE_CONTRACT = 1


@functools.cache
def _loaded_engine() -> ModuleType | None:
    """``gbdraw_core`` when it is present and has this package's contract."""

    try:
        engine = importlib.import_module("gbdraw_core")
    except ImportError:
        return None
    contract = getattr(engine, "CONTRACT_VERSION", None)
    if contract != ENGINE_CONTRACT:
        logger.warning(
            "WARNING: gbdraw_core has contract %s, not %s: Python draws every group. Rebuild or reinstall gbdraw_core.",
            contract,
            ENGINE_CONTRACT,
        )
        return None
    return engine


def native_engine() -> ModuleType | None:
    """``gbdraw_core``, or ``None`` when Python draws everything."""

    if os.environ.get("GBDRAW_NATIVE") == "0":
        return None
    return _loaded_engine()


__all__ = ["ENGINE_CONTRACT", "native_engine"]
