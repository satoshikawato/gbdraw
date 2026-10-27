#!/usr/bin/env python
# coding: utf-8

"""Native LOSAT and NCBI BLAST+ runtime owner for every comparison program.

LOSATN, TLOSATX, and LOSATP share one resolution order, one argv builder, one
subprocess path, and one raw-cache store. Programs differ only by the data in
``LOSAT_PROGRAMS``. Program identity (protein manifests, nucleotide hashes)
and result parsing stay with their analysis owners. Installation and the
managed cache stay in ``gbdraw.losat_setup``.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
import logging
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from types import MappingProxyType
from typing import Callable, Literal, Mapping, Sequence, TypedDict, cast

from gbdraw.comparisons.losat_jobs import losat_entry_filename
from gbdraw.exceptions import ValidationError

logger = logging.getLogger(__name__)

# The string "losat" is the automatic-selection token of the executable option.
AUTOMATIC_LOSAT_BIN = "losat"
# Current CLI spellings of the runtime overrides, used in diagnostics.
LOSAT_BIN_OPTION = "--losat_bin"
NCBI_BLAST_BIN_OPTION = "--ncbi_blast_bin"
_PROBE_TIMEOUT_SECONDS = 30
_UNAVAILABLE: Mapping[str, str] = MappingProxyType(
    {"code": "LOSAT_RUNTIME", "reason": "UNAVAILABLE"}
)

LosatProgramName = Literal["losatn", "tlosatx", "losatp"]
LosatRuntimeKind = Literal["losat", "ncbi-blast"]
# Sessions written by gbdraw 0.14.2 or earlier may also record "bundled" (a
# development build once tracked in the source tree). Python never reads the
# record back; the Web shows its source and path as stored.
LosatRuntimeSource = Literal["explicit", "conda", "managed", "path"]


class LosatRuntimeRecord(TypedDict):
    """The ``runtime`` record a Session stores beside a raw LOSAT entry.

    It is not part of the raw key. Python copies stored records as they are;
    the Web run information reads ``source`` and ``path``.
    """

    kind: LosatRuntimeKind
    # None when the ``--version`` probe failed.
    version: str | None
    source: LosatRuntimeSource
    # The executable's file name only (a shared Session carries no directories).
    path: str
    # The LOSAT subcommand or NCBI BLAST+ executable name.
    program: str


LosatRuntimeCallback = Callable[[LosatRuntimeRecord], None]


@dataclass(frozen=True)
class LosatOptionFlags:
    """One search option spelled for the runtime CLIs and for the raw-cache key."""

    # LOSAT CLI v2 (released v0.1.0) and NCBI BLAST+ share this spelling.
    search: str
    cache: str


# Cache args keep the Web v1 form so CLI and Web raw keys match.
LOSAT_OPTION_FLAGS: Mapping[str, LosatOptionFlags] = MappingProxyType(
    {
        "task": LosatOptionFlags("-task", "--task"),
        "query_gencode": LosatOptionFlags("-query_gencode", "--query-gencode"),
        "db_gencode": LosatOptionFlags("-db_gencode", "--db-gencode"),
        "max_hsps": LosatOptionFlags("-max_hsps", "--max-hsps-per-subject"),
        "max_target_seqs": LosatOptionFlags("-max_target_seqs", "--max-target-seqs"),
    }
)


@dataclass(frozen=True)
class LosatProgram:
    """One comparison program as data."""

    name: LosatProgramName
    # LOSAT subcommand, NCBI BLAST+ executable, and raw-cache ``program``.
    search: str
    description: str
    activity: str
    options: tuple[str, ...]
    required_options: frozenset[str]
    fasta_suffix: str


LOSAT_PROGRAMS: Mapping[str, LosatProgram] = MappingProxyType(
    {
        "losatn": LosatProgram(
            "losatn",
            "blastn",
            "nucleotide BLASTN",
            "nucleotide comparison",
            ("task",),
            frozenset({"task"}),
            ".fna",
        ),
        "tlosatx": LosatProgram(
            "tlosatx",
            "tblastx",
            "translated TBLASTX",
            "translated nucleotide comparison",
            ("query_gencode", "db_gencode"),
            # An unset table is the runtime default (1), as in the Web (D17).
            frozenset(),
            ".fna",
        ),
        "losatp": LosatProgram(
            "losatp",
            "blastp",
            "protein BLASTP",
            "protein colinearity",
            ("max_hsps", "max_target_seqs"),
            frozenset(),
            ".faa",
        ),
    }
)


@dataclass(frozen=True)
class LosatSearchArgs:
    """Program-specific search options; unset options are omitted."""

    task: str | None = None
    query_gencode: int | None = None
    db_gencode: int | None = None
    max_hsps: int | None = None
    max_target_seqs: int | None = None


@dataclass(frozen=True)
class LosatRuntime:
    """Resolved executable for a comparison search."""

    kind: LosatRuntimeKind
    executable: str
    source: LosatRuntimeSource


def losat_program(program: str | LosatProgram) -> LosatProgram:
    if isinstance(program, LosatProgram):
        return program
    spec = LOSAT_PROGRAMS.get(str(program))
    if spec is None:
        raise ValidationError(
            f"Unknown LOSAT program {program!r}; expected one of: "
            + ", ".join(LOSAT_PROGRAMS),
            diagnostic={"code": "COMPARISON_INPUT"},
        )
    return spec


def _option_values(spec: LosatProgram, options: LosatSearchArgs) -> list[tuple[str, str]]:
    for name in LOSAT_OPTION_FLAGS:
        if name not in spec.options and getattr(options, name) is not None:
            raise ValidationError(
                f"LOSAT option {name} does not apply to {spec.name}.",
                diagnostic={"code": "COMPARISON_INPUT"},
            )
    values: list[tuple[str, str]] = []
    for name in spec.options:
        value = getattr(options, name)
        if value is None:
            if name in spec.required_options:
                raise ValidationError(
                    f"{spec.name} requires the LOSAT option {name}.",
                    diagnostic={"code": "COMPARISON_INPUT", "reason": "REQUIRED"},
                )
            continue
        values.append((name, str(value) if name == "task" else str(int(value))))
    return values


def losat_cache_args(
    program: str | LosatProgram,
    options: LosatSearchArgs,
) -> list[str]:
    """Return raw-cache key args in the Web v1 form."""

    spec = losat_program(program)
    args: list[str] = []
    for name, value in _option_values(spec, options):
        args.extend([LOSAT_OPTION_FLAGS[name].cache, value])
    return args


# Runtime resolution.


def _path_executable(name: str) -> str | None:
    return shutil.which(str(name).strip())


def _conda_losat_runtime() -> LosatRuntime | None:
    prefix = Path(sys.prefix)
    if not (prefix / "conda-meta").is_dir():
        return None

    candidate = (prefix / "bin" / "losat").absolute()
    try:
        mode = candidate.stat().st_mode
    except FileNotFoundError:
        if candidate.is_symlink():
            raise ValidationError(
                f"Conda LOSAT candidate is a broken symbolic link: {candidate}",
                diagnostic=_UNAVAILABLE,
            )
        return None
    except OSError as exc:
        raise ValidationError(
            f"Could not inspect the conda LOSAT candidate at {candidate}: {exc}",
            diagnostic=_UNAVAILABLE,
        ) from exc

    if not stat.S_ISREG(mode):
        raise ValidationError(
            f"Conda LOSAT candidate is not a regular file: {candidate}",
            diagnostic=_UNAVAILABLE,
        )
    if not os.access(candidate, os.X_OK):
        raise ValidationError(
            f"Conda LOSAT candidate is not executable: {candidate}",
            diagnostic=_UNAVAILABLE,
        )

    return LosatRuntime("losat", str(candidate), "conda")


def _runtime_label(runtime: LosatRuntime, spec: LosatProgram) -> str:
    if runtime.kind == "ncbi-blast":
        return f"NCBI BLAST+ {spec.search}"
    return f"LOSAT {spec.search}"


def _runtime_unavailable_error(spec: LosatProgram) -> ValidationError:
    return ValidationError(
        f"{spec.description[:1].upper()}{spec.description[1:]} comparison needs "
        "LOSAT or NCBI BLAST+. "
        "No managed LOSAT is installed, "
        "`losat` was not found on PATH, and "
        f"`{spec.search}` was not found on PATH. "
        "Run gbdraw setup-losat to install the pinned release where it is supported. "
        f"Install NCBI BLAST+ and make `{spec.search}` available on PATH, or pass a "
        f"native LOSAT executable with {LOSAT_BIN_OPTION}, or pass an NCBI BLAST+ "
        f"executable with {NCBI_BLAST_BIN_OPTION}.",
        diagnostic={**_UNAVAILABLE, "program": spec.name},
    )


def resolve_losat_runtime(
    program: str | LosatProgram,
    *,
    losat_bin: str | None = None,
    ncbi_blast_bin: str | None = None,
) -> LosatRuntime:
    """Resolve the runtime for one program in the shared order.

    Explicit LOSAT, explicit NCBI BLAST+, conda LOSAT, managed LOSAT, PATH
    ``losat``, then PATH NCBI ``<program>``. ``run_losat_search`` logs the
    chosen runtime.
    """

    spec = losat_program(program)
    requested_bin = str(losat_bin or AUTOMATIC_LOSAT_BIN).strip() or AUTOMATIC_LOSAT_BIN
    requested_ncbi_bin = str(ncbi_blast_bin or "").strip() or None
    if requested_bin != AUTOMATIC_LOSAT_BIN and requested_ncbi_bin is not None:
        raise ValidationError(
            f"Pass either {LOSAT_BIN_OPTION} or {NCBI_BLAST_BIN_OPTION} for "
            f"{spec.description} comparisons, not both.",
            diagnostic={"code": "COMPARISON_INPUT"},
        )
    if requested_bin != AUTOMATIC_LOSAT_BIN:
        return LosatRuntime("losat", requested_bin, "explicit")
    if requested_ncbi_bin is not None:
        return LosatRuntime("ncbi-blast", requested_ncbi_bin, "explicit")

    conda_runtime = _conda_losat_runtime()
    if conda_runtime is not None:
        return conda_runtime

    from gbdraw.losat_setup import managed_losat

    managed_path = managed_losat()
    if managed_path is not None:
        return LosatRuntime("losat", str(managed_path), "managed")

    path_losat = _path_executable(AUTOMATIC_LOSAT_BIN)
    if path_losat is not None:
        return LosatRuntime("losat", path_losat, "path")

    path_ncbi = _path_executable(spec.search)
    if path_ncbi is not None:
        return LosatRuntime("ncbi-blast", path_ncbi, "path")

    raise _runtime_unavailable_error(spec)


# Runtime identity: version, probed once per executable.

_RUNTIME_VERSIONS: dict[tuple[str, str], str | None] = {}


def _probe(command: list[str]) -> subprocess.CompletedProcess[str] | None:
    try:
        return subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=_PROBE_TIMEOUT_SECONDS,
        )
    except (OSError, subprocess.SubprocessError):
        return None


_TASK_VALUES: dict[str, frozenset[str] | None] = {}
_POSSIBLE_VALUES = re.compile(r"\[possible values:\s*([^\]]+)\]")


def losat_task_values(executable: str) -> frozenset[str] | None:
    """Return the ``blastn`` task values a LOSAT runtime lists, or ``None``.

    Released LOSAT v0.1.0 lists ``[possible values: megablast, blastn]`` under
    its task option; a runtime that lists none is not checked. Probed once per
    executable.
    """

    key = str(executable)
    if key in _TASK_VALUES:
        return _TASK_VALUES[key]
    completed = _probe([key, "blastn", "--help"])
    values: frozenset[str] | None = None
    if completed is not None:
        lines = f"{completed.stdout}\n{completed.stderr}".splitlines()
        for index, line in enumerate(lines):
            if re.match(r"^\s*-{1,2}task\b", line):
                block = "\n".join(lines[index:index + 8])
                match = _POSSIBLE_VALUES.search(block)
                if match:
                    values = frozenset(
                        item.strip() for item in match.group(1).split(",") if item.strip()
                    )
                break
    _TASK_VALUES[key] = values
    return values


def require_losat_task_support(
    task: str,
    *,
    losat_bin: str | None = None,
    ncbi_blast_bin: str | None = None,
) -> None:
    """Fail before searching when the resolved LOSAT runtime lacks ``task`` (D11).

    NCBI BLAST+ ``blastn`` accepts every LOSATN task.
    """

    spec = losat_program("losatn")
    runtime = resolve_losat_runtime(spec, losat_bin=losat_bin, ncbi_blast_bin=ncbi_blast_bin)
    if runtime.kind != "losat":
        return
    values = losat_task_values(runtime.executable)
    if values is None or task in values:
        return
    version = _runtime_version(runtime) or "unknown version"
    raise ValidationError(
        f"LOSATN task {task!r} is not supported by the LOSAT runtime at "
        f"{runtime.executable} ({version}; supported: "
        f"{', '.join(sorted(values))}). Choose another --losatn_task, or pass a "
        f"runtime that supports it with {LOSAT_BIN_OPTION} or {NCBI_BLAST_BIN_OPTION}.",
        diagnostic={
            "code": "COMPARISON_INPUT",
            "reason": "LOSAT_TASK",
            "field": "losatn_task",
            "program": spec.name,
        },
    )


def _runtime_version(runtime: LosatRuntime) -> str | None:
    key = (runtime.kind, runtime.executable)
    if key in _RUNTIME_VERSIONS:
        return _RUNTIME_VERSIONS[key]
    flag = "--version" if runtime.kind == "losat" else "-version"
    completed = _probe([runtime.executable, flag])
    if completed is None:
        return None
    version: str | None = None
    if completed.returncode == 0:
        first_line = next(
            (line.strip() for line in completed.stdout.splitlines() if line.strip()),
            "",
        )
        # "losat 0.1.0" or "blastn: 2.16.0+"
        match = re.match(r"^\S+?:?\s+(\S+)", first_line)
        version = match.group(1) if match else (first_line[:200] or None)
    _RUNTIME_VERSIONS[key] = version
    return version


def _recorded_runtime_path(runtime: LosatRuntime) -> str:
    # Sessions are shared, so a runtime records only its executable's name;
    # source and version identify it.
    return Path(runtime.executable).name or runtime.executable


def losat_runtime_record(
    runtime: LosatRuntime,
    program: str | LosatProgram,
) -> LosatRuntimeRecord:
    """Return the non-key runtime identity stored beside a raw cache entry.

    Records written by gbdraw 0.14.2 or earlier may carry a ``cli`` dialect
    field; it is kept as stored and nothing reads it.
    """

    spec = losat_program(program)
    return {
        "kind": runtime.kind,
        "version": _runtime_version(runtime),
        "source": runtime.source,
        "path": _recorded_runtime_path(runtime),
        "program": spec.search,
    }


# Argv and execution.


def build_losat_command(
    runtime: LosatRuntime,
    program: str | LosatProgram,
    *,
    query_path: Path,
    subject_path: Path,
    options: LosatSearchArgs = LosatSearchArgs(),
    threads: int | None = None,
) -> list[str]:
    """Build the outfmt 6 argv for one search."""

    spec = losat_program(program)
    values = _option_values(spec, options)
    if runtime.kind == "ncbi-blast":
        command = [str(runtime.executable)]
    else:
        command = [str(runtime.executable), spec.search]
    command.extend(
        [
            "-query",
            str(query_path),
            "-subject",
            str(subject_path),
            "-outfmt",
            "6",
        ]
    )
    for name, value in values:
        command.extend([LOSAT_OPTION_FLAGS[name].search, value])
    if threads is not None:
        command.extend(["-num_threads", str(int(threads))])
    return command


def _run_losat_subprocess(
    command: list[str],
    *,
    runtime_label: str,
) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise ValidationError(
            f"{runtime_label} executable not found: {command[0]}",
            diagnostic=_UNAVAILABLE,
        ) from exc
    except PermissionError as exc:
        raise ValidationError(
            f"{runtime_label} executable is not executable: {command[0]}",
            diagnostic=_UNAVAILABLE,
        ) from exc
    except OSError as exc:
        raise ValidationError(
            f"{runtime_label} executable could not be started at {command[0]}: {exc}",
            diagnostic={"code": "LOSAT_RUNTIME", "reason": "FAILED"},
        ) from exc


def run_losat_search(
    program: str | LosatProgram,
    query_fasta: str,
    subject_fasta: str,
    *,
    options: LosatSearchArgs = LosatSearchArgs(),
    losat_bin: str | None = None,
    ncbi_blast_bin: str | None = None,
    threads: int | None = None,
    runtime_callback: LosatRuntimeCallback | None = None,
) -> str:
    """Run one search and return its raw outfmt 6 text.

    ``runtime_callback`` receives the runtime record after a successful search.
    """

    spec = losat_program(program)
    runtime_record: LosatRuntimeRecord | None = None
    runtime = resolve_losat_runtime(spec, losat_bin=losat_bin, ncbi_blast_bin=ncbi_blast_bin)
    with tempfile.TemporaryDirectory(prefix="gbdraw_losat_") as temp_dir:
        temp_path = Path(temp_dir)
        query_path = temp_path / f"query{spec.fasta_suffix}"
        subject_path = temp_path / f"subject{spec.fasta_suffix}"
        query_path.write_text(query_fasta, encoding="utf-8")
        subject_path.write_text(subject_fasta, encoding="utf-8")

        command = build_losat_command(
            runtime,
            spec,
            query_path=query_path,
            subject_path=subject_path,
            options=options,
            threads=threads,
        )
        runtime_label = _runtime_label(runtime, spec)
        # The one log line that names the backend. The version probe runs
        # once per executable; the runtime record reuses it.
        version = _runtime_version(runtime)
        logger.info(
            "INFO: Running %s for %s (backend %s, source %s, executable %s%s).",
            runtime_label,
            spec.activity,
            runtime.kind,
            runtime.source,
            runtime.executable,
            f", version {version}" if version else "",
        )
        completed = _run_losat_subprocess(command, runtime_label=runtime_label)
        if completed.returncode == 0 and runtime_callback is not None:
            runtime_record = losat_runtime_record(runtime, spec)

    if completed.returncode != 0:
        stderr = completed.stderr.strip()
        detail = f": {stderr}" if stderr else ""
        raise ValidationError(
            f"{runtime_label} failed at {runtime.executable} "
            f"with exit code {completed.returncode}{detail}",
            diagnostic={
                "code": "LOSAT_RUNTIME",
                "reason": "FAILED",
                "program": spec.name,
                "exitCode": completed.returncode,
            },
        )
    if runtime_callback is not None and runtime_record is not None:
        runtime_callback(runtime_record)
    return completed.stdout


# Raw cache store.


LosatCacheWriteKind = Literal["entry", "display", "retire-legacy"]


@dataclass(frozen=True)
class LosatCacheWrite:
    """One write a lookup or a search made into a loaded raw cache.

    ``entry`` stores ``entry`` under ``key``; ``display`` shows ``key`` under
    ``filename``; ``retire-legacy`` drops legacy candidate ``index`` (protein
    caches). Replayed in order on a cache in the same state, the writes give
    the same Session entries (``LosatRawCache.apply_writes``).
    """

    kind: LosatCacheWriteKind
    key: str = ""
    entry: Mapping[str, object] | None = None
    filename: str = ""
    index: int = -1


class LosatRawCache:
    """Raw LOSAT entries by raw key, kept in Session display order.

    Program identity owners derive keys and validate entries. This store keeps
    the entries, their display filenames, and the non-key ``runtime`` record.
    Every write after loading is also kept in ``writes``.
    """

    def __init__(self) -> None:
        self._entries_by_key: dict[str, dict[str, object]] = {}
        self._display_order: list[str] = []
        self._display_info: dict[str, tuple[str, bool]] = {}
        self._writes: list[LosatCacheWrite] = []

    @property
    def writes(self) -> tuple[LosatCacheWrite, ...]:
        """The writes made since loading, in order."""

        return tuple(self._writes)

    def apply_writes(self, writes: Sequence[LosatCacheWrite]) -> None:
        """Replay writes another cache in this cache's state made."""

        for write in writes:
            self._apply_write(write)

    def _apply_write(self, write: LosatCacheWrite) -> None:
        if write.kind == "entry" and write.entry is not None:
            self._store_entry(write.key, dict(write.entry))
        elif write.kind == "display":
            self._mark_display(write.key, write.filename)
        else:
            raise ValueError(f"{type(self).__name__} cannot apply a {write.kind!r} write.")

    def _raw_state(
        self,
    ) -> tuple[tuple[tuple[str, Mapping[str, object]], ...], tuple[tuple[str, str], ...]]:
        """The stored entries and the displayed keys with their filenames, in order."""

        return (
            tuple(self._entries_by_key.items()),
            tuple((key, self._display_info[key][0]) for key in self._display_order),
        )

    def _restore_raw_state(
        self,
        entries: Sequence[tuple[str, Mapping[str, object]]],
        display: Sequence[tuple[str, str]],
    ) -> None:
        # Entries are shared, not copied: a write stores a new dict.
        self._entries_by_key = {key: cast(dict, entry) for key, entry in entries}
        self._display_order = [key for key, _filename in display]
        self._display_info = {key: (filename, True) for key, filename in display}

    @property
    def has_entries(self) -> bool:
        return bool(self._entries_by_key)

    def _cached_entry(self, key: str) -> dict[str, object] | None:
        return self._entries_by_key.get(key)

    def cached_entry(self, key: str) -> dict[str, object] | None:
        """The entry stored under a raw key, if any."""

        return self._cached_entry(key)

    def store_search_entry(
        self,
        key: str,
        entry: dict[str, object],
        *,
        runtime: Mapping[str, object] | None,
        filename: str,
    ) -> None:
        """Store a displayed entry that a search produced in this run."""

        self._add_search_entry(key, entry, runtime=runtime, display=True, filename=filename)

    def _add_loaded_entry(
        self,
        key: str,
        entry: dict[str, object],
        source: Mapping[str, object],
    ) -> None:
        """Store a Session entry; the first displayed copy keeps its filename."""

        runtime = source.get("runtime")
        if isinstance(runtime, Mapping):
            entry["runtime"] = copy.deepcopy(dict(runtime))
        self._entries_by_key[key] = entry
        if entry.get("display") is not False and key not in self._display_info:
            self._display_order.append(key)
            self._display_info[key] = (losat_entry_filename(entry), True)

    def _add_search_entry(
        self,
        key: str,
        entry: dict[str, object],
        *,
        runtime: Mapping[str, object] | None = None,
        display: bool = False,
        filename: str = "",
    ) -> None:
        if runtime is not None:
            entry["runtime"] = dict(runtime)
        self._store_entry(key, entry)
        if display:
            self._mark_display(key, filename)

    def _store_entry(self, key: str, entry: dict[str, object]) -> None:
        self._entries_by_key[key] = entry
        self._writes.append(LosatCacheWrite("entry", key, entry=dict(entry)))

    def _retain_entries(self, keep: Callable[[dict[str, object]], bool]) -> None:
        self._entries_by_key = {
            key: entry
            for key, entry in self._entries_by_key.items()
            if keep(entry)
        }
        self._display_order = [
            key for key in self._display_order if key in self._entries_by_key
        ]
        self._display_info = {
            key: value
            for key, value in self._display_info.items()
            if key in self._entries_by_key
        }

    def _mark_display(self, key: str, filename: str) -> None:
        if key not in self._display_info:
            self._display_order.append(key)
        self._display_info[key] = (str(filename or ""), True)
        self._writes.append(LosatCacheWrite("display", key, filename=str(filename or "")))

    def session_entries(self) -> tuple[dict[str, object], ...]:
        result: list[dict[str, object]] = []
        seen: set[str] = set()
        for key in self._display_order:
            entry = self._entries_by_key.get(key)
            if entry is None:
                continue
            filename, display = self._display_info.get(key, ("", True))
            rendered = dict(entry)
            rendered["filename"] = filename
            rendered["display"] = display
            result.append(rendered)
            seen.add(key)
        for key, entry in self._entries_by_key.items():
            if key in seen:
                continue
            rendered = dict(entry)
            rendered["filename"] = ""
            rendered["display"] = False
            result.append(rendered)
        return tuple(result)


__all__ = [
    "AUTOMATIC_LOSAT_BIN",
    "LOSAT_BIN_OPTION",
    "LOSAT_OPTION_FLAGS",
    "LOSAT_PROGRAMS",
    "LosatCacheWrite",
    "LosatCacheWriteKind",
    "LosatProgram",
    "LosatRawCache",
    "LosatRuntime",
    "LosatRuntimeCallback",
    "LosatRuntimeRecord",
    "LosatSearchArgs",
    "NCBI_BLAST_BIN_OPTION",
    "build_losat_command",
    "losat_cache_args",
    "losat_program",
    "losat_runtime_record",
    "losat_task_values",
    "require_losat_task_support",
    "resolve_losat_runtime",
    "run_losat_search",
]
