"""Add the bytecode of every module to the browser wheel.

Pyodide compiles each module of the gbdraw wheel from source on every Worker
boot, because the wheel carries no ``.pyc`` and its Emscripten file system does
not keep one between page loads. ``add_bytecode`` writes, next to each ``.py``,
an unchecked-hash ``.pyc`` (PEP 552) compiled by the running interpreter, and
lists it in RECORD. A bytecode file is valid only for the interpreter that
wrote it, so ``tools/compile_browser_wheel_bytecode.mjs`` runs this module in
the vendored Pyodide, whatever Python builds the wheel. The ``.py`` sources stay:
tracebacks and ``inspect`` read them.

Usage inside an interpreter that should own the bytecode:
    add_bytecode(source_wheel, target_wheel)
"""

from __future__ import annotations

import base64
import hashlib
import importlib.util
import sys
import zipfile
from importlib._bootstrap_external import _code_to_hash_pyc


def _record_hash(data: bytes) -> str:
    digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode("ascii")
    return f"sha256={digest}"


def _bytecode_path(source_path: str) -> str:
    folder, _, name = source_path.rpartition("/")
    prefix = f"{folder}/" if folder else ""
    return f"{prefix}__pycache__/{name[:-3]}.{sys.implementation.cache_tag}.pyc"


def add_bytecode(source_wheel: str, target_wheel: str) -> int:
    """Copy ``source_wheel`` to ``target_wheel`` with a ``.pyc`` per ``.py``; return their count.

    Entries keep their order, timestamps and attributes; each ``.pyc`` follows
    its source and takes its timestamp, so the same wheel gives the same bytes.
    """

    added: list[str] = []
    with zipfile.ZipFile(source_wheel) as source, zipfile.ZipFile(target_wheel, "w") as target:
        record_name = next(name for name in source.namelist() if name.endswith(".dist-info/RECORD"))
        for info in source.infolist():
            if info.filename == record_name:
                continue
            data = source.read(info)
            target.writestr(info, data)
            if not info.filename.endswith(".py"):
                continue
            code = compile(data, info.filename, "exec", dont_inherit=True)
            bytecode = bytes(_code_to_hash_pyc(code, importlib.util.source_hash(data), checked=False))
            entry = zipfile.ZipInfo(_bytecode_path(info.filename), info.date_time)
            entry.external_attr = info.external_attr
            entry.compress_type = zipfile.ZIP_DEFLATED
            target.writestr(entry, bytecode)
            added.append(f"{entry.filename},{_record_hash(bytecode)},{len(bytecode)}")
        record = source.read(record_name).decode("utf-8").splitlines()
        rows = [row for row in record if row and not row.startswith(f"{record_name},")]
        target.writestr(source.getinfo(record_name), "\n".join([*rows, *added, f"{record_name},,"]) + "\n")
    return len(added)
