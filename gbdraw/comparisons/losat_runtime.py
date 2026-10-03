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
from contextlib import ExitStack
from dataclasses import dataclass
from importlib import resources
import logging
import os
import platform
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from types import MappingProxyType
from typing import Callable, Literal, Mapping

from gbdraw.exceptions import ValidationError

logger = logging.getLogger(__name__)

# The string "losat" is the automatic-selection token of the executable option.
AUTOMATIC_LOSAT_BIN = "losat"
# Current CLI spellings of the runtime overrides, used in diagnostics.
LOSAT_BIN_OPTION = "--losatp_bin"
NCBI_BLAST_BIN_OPTION = "--ncbi_blastp_bin"
_BUNDLED_LOSAT_DIR = "bin"
_PROBE_TIMEOUT_SECONDS = 30

LosatProgramName = Literal["losatn", "tlosatx", "losatp"]
LosatRuntimeKind = Literal["losat", "ncbi-blast"]
LosatRuntimeSource = Literal["explicit", "conda", "managed", "bundled", "path"]
LosatCliDialect = Literal["v1", "v2"]
LosatRuntimeCallback = Callable[[dict[str, object]], None]


@dataclass(frozen=True)
class LosatOptionFlags:
    """One search option spelled for each runtime CLI and for the raw-cache key."""

    losat_v2: str
    losat_v1: str
    ncbi: str
    cache: str


# LOSAT CLI v2 is the released v0.1.0 spelling. CLI v1 is the tracked
# gbdraw/bin development build, which spells genetic codes with double dashes.
# Cache args keep the Web v1 form so CLI and Web raw keys match.
LOSAT_OPTION_FLAGS: Mapping[str, LosatOptionFlags] = MappingProxyType(
    {
        "task": LosatOptionFlags("-task", "-task", "-task", "--task"),
        "query_gencode": LosatOptionFlags(
            "-query_gencode", "--query-gencode", "-query_gencode", "--query-gencode"
        ),
        "db_gencode": LosatOptionFlags(
            "-db_gencode", "--db-gencode", "-db_gencode", "--db-gencode"
        ),
        "max_hsps": LosatOptionFlags(
            "-max_hsps", "-max_hsps", "-max_hsps", "--max-hsps-per-subject"
        ),
        "max_target_seqs": LosatOptionFlags(
            "-max_target_seqs", "-max_target_seqs", "-max_target_seqs", "--max-target-seqs"
        ),
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
            frozenset({"query_gencode", "db_gencode"}),
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
class LosatSearchOptions:
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
            + ", ".join(LOSAT_PROGRAMS)
        )
    return spec


def _option_values(spec: LosatProgram, options: LosatSearchOptions) -> list[tuple[str, str]]:
    for name in LOSAT_OPTION_FLAGS:
        if name not in spec.options and getattr(options, name) is not None:
            raise ValidationError(f"LOSAT option {name} does not apply to {spec.name}.")
    values: list[tuple[str, str]] = []
    for name in spec.options:
        value = getattr(options, name)
        if value is None:
            if name in spec.required_options:
                raise ValidationError(f"{spec.name} requires the LOSAT option {name}.")
            continue
        values.append((name, str(value) if name == "task" else str(int(value))))
    return values


def losat_cache_args(
    program: str | LosatProgram,
    options: LosatSearchOptions,
) -> list[str]:
    """Return raw-cache key args in the Web v1 form."""

    spec = losat_program(program)
    args: list[str] = []
    for name, value in _option_values(spec, options):
        args.extend([LOSAT_OPTION_FLAGS[name].cache, value])
    return args


# Platform detection and runtime resolution.


def _normalize_machine(machine: str | None = None) -> str:
    normalized = str(machine or platform.machine()).strip().lower()
    aliases = {
        "amd64": "x86_64",
        "x64": "x86_64",
        "arm64": "aarch64",
    }
    return aliases.get(normalized, normalized)


def _bundled_platform_dir() -> str | None:
    machine = _normalize_machine()
    if sys.platform.startswith("linux"):
        if machine == "x86_64":
            return "linux-x86_64"
        if machine == "aarch64":
            return "linux-aarch64"
    if sys.platform == "darwin":
        if machine == "x86_64":
            return "macos-x86_64"
        if machine == "aarch64":
            return "macos-arm64"
    if os.name == "nt":
        if machine == "x86_64":
            return "windows-x86_64"
        if machine == "aarch64":
            return "windows-arm64"
    return None


def _bundled_filename() -> str:
    return "losat.exe" if os.name == "nt" else "losat"


def _bundled_losat_resource():
    platform_dir = _bundled_platform_dir()
    if platform_dir is None:
        return None
    try:
        package_root = resources.files("gbdraw")
    except (ModuleNotFoundError, FileNotFoundError):
        return None
    candidate = (
        package_root
        .joinpath(_BUNDLED_LOSAT_DIR)
        .joinpath(platform_dir)
        .joinpath(_bundled_filename())
    )
    if not candidate.is_file():
        return None
    return candidate


def _ensure_executable(path: Path) -> None:
    if os.name == "nt":
        return
    try:
        mode = path.stat().st_mode
    except OSError:
        return
    if mode & stat.S_IXUSR:
        return
    try:
        path.chmod(mode | stat.S_IXUSR)
    except OSError:
        logger.debug("Could not mark bundled LOSAT binary executable: %s", path)


def bundled_losat_runtime(*, stack: ExitStack) -> LosatRuntime | None:
    """Return the source-checkout LOSAT for this platform, marked executable."""

    bundled_resource = _bundled_losat_resource()
    if bundled_resource is None:
        return None
    bundled_path = stack.enter_context(resources.as_file(bundled_resource))
    _ensure_executable(bundled_path)
    return LosatRuntime("losat", str(bundled_path), "bundled")


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
                f"Conda LOSAT candidate is a broken symbolic link: {candidate}"
            )
        return None
    except OSError as exc:
        raise ValidationError(
            f"Could not inspect the conda LOSAT candidate at {candidate}: {exc}"
        ) from exc

    if not stat.S_ISREG(mode):
        raise ValidationError(
            f"Conda LOSAT candidate is not a regular file: {candidate}"
        )
    if not os.access(candidate, os.X_OK):
        raise ValidationError(
            f"Conda LOSAT candidate is not executable: {candidate}"
        )

    runtime = LosatRuntime("losat", str(candidate), "conda")
    logger.debug("Using conda LOSAT runtime: %s", runtime.executable)
    return runtime


def _runtime_label(runtime: LosatRuntime, spec: LosatProgram) -> str:
    if runtime.kind == "ncbi-blast":
        return f"NCBI BLAST+ {spec.search}"
    return f"LOSAT {spec.search}"


def _runtime_unavailable_error(spec: LosatProgram, *, platform_dir: str | None) -> ValidationError:
    platform_name = platform_dir or "this platform"
    platform_note = ""
    if platform_dir and (platform_dir.startswith("macos-") or platform_dir.startswith("windows-")):
        platform_note = (
            " gbdraw does not currently include a LOSAT binary for this platform."
        )
    return ValidationError(
        f"{spec.description[:1].upper()}{spec.description[1:]} comparison needs "
        "LOSAT or NCBI BLAST+. "
        f"no bundled LOSAT binary was found for {platform_name}, "
        "`losat` was not found on PATH, and "
        f"`{spec.search}` was not found on PATH. "
        "Run gbdraw setup-losat to install the pinned release when available."
        f"{platform_note} "
        f"Install NCBI BLAST+ and make `{spec.search}` available on PATH, or pass a "
        f"native LOSAT executable with {LOSAT_BIN_OPTION}, or pass an NCBI BLAST+ "
        f"executable with {NCBI_BLAST_BIN_OPTION}."
    )


def resolve_losat_runtime(
    program: str | LosatProgram,
    *,
    losat_bin: str | None = None,
    ncbi_blast_bin: str | None = None,
    stack: ExitStack,
) -> LosatRuntime:
    """Resolve the runtime for one program in the shared order.

    Explicit LOSAT, explicit NCBI BLAST+, conda LOSAT, managed LOSAT, bundled
    LOSAT (source checkout), PATH ``losat``, then PATH NCBI ``<program>``.
    """

    spec = losat_program(program)
    requested_bin = str(losat_bin or AUTOMATIC_LOSAT_BIN).strip() or AUTOMATIC_LOSAT_BIN
    requested_ncbi_bin = str(ncbi_blast_bin or "").strip() or None
    if requested_bin != AUTOMATIC_LOSAT_BIN and requested_ncbi_bin is not None:
        raise ValidationError(
            f"Pass either {LOSAT_BIN_OPTION} or {NCBI_BLAST_BIN_OPTION} for "
            f"{spec.description} comparisons, not both."
        )
    if requested_bin != AUTOMATIC_LOSAT_BIN:
        runtime = LosatRuntime("losat", requested_bin, "explicit")
        logger.debug("Using explicit LOSAT runtime: %s", runtime.executable)
        return runtime
    if requested_ncbi_bin is not None:
        runtime = LosatRuntime("ncbi-blast", requested_ncbi_bin, "explicit")
        logger.debug("Using explicit NCBI BLAST+ runtime: %s", runtime.executable)
        return runtime

    conda_runtime = _conda_losat_runtime()
    if conda_runtime is not None:
        return conda_runtime

    from gbdraw.losat_setup import managed_losat

    managed_path = managed_losat()
    if managed_path is not None:
        return LosatRuntime("losat", str(managed_path), "managed")

    bundled = bundled_losat_runtime(stack=stack)
    if bundled is not None:
        logger.debug("Using bundled LOSAT runtime: %s", bundled.executable)
        return bundled

    path_losat = _path_executable(AUTOMATIC_LOSAT_BIN)
    if path_losat is not None:
        runtime = LosatRuntime("losat", path_losat, "path")
        logger.debug("Using PATH LOSAT runtime: %s", runtime.executable)
        return runtime

    path_ncbi = _path_executable(spec.search)
    if path_ncbi is not None:
        runtime = LosatRuntime("ncbi-blast", path_ncbi, "path")
        logger.debug("Using PATH NCBI BLAST+ %s runtime: %s", spec.search, runtime.executable)
        return runtime

    raise _runtime_unavailable_error(spec, platform_dir=_bundled_platform_dir())


# Runtime identity: CLI dialect and version, probed once per executable.

_CLI_DIALECTS: dict[str, LosatCliDialect] = {}
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


def detect_losat_cli_dialect(executable: str) -> LosatCliDialect:
    """Return ``v1`` (double-dash gencode flags) or ``v2`` (released v0.1.0)."""

    key = str(executable)
    cached = _CLI_DIALECTS.get(key)
    if cached is not None:
        return cached
    completed = _probe([key, "tblastx", "--help"])
    if completed is None:
        # The search reports the start failure with its own diagnostic.
        return "v2"
    help_text = f"{completed.stdout}\n{completed.stderr}"
    flags = LOSAT_OPTION_FLAGS["query_gencode"]
    dialect: LosatCliDialect = (
        "v1" if flags.losat_v1 in help_text and flags.losat_v2 not in help_text else "v2"
    )
    _CLI_DIALECTS[key] = dialect
    return dialect


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
    if runtime.source == "bundled":
        platform_dir = _bundled_platform_dir()
        if platform_dir is not None:
            return "/".join(("gbdraw", _BUNDLED_LOSAT_DIR, platform_dir, _bundled_filename()))
    return runtime.executable


def losat_runtime_record(
    runtime: LosatRuntime,
    program: str | LosatProgram,
) -> dict[str, object]:
    """Return the non-key runtime identity stored beside a raw cache entry."""

    spec = losat_program(program)
    record: dict[str, object] = {
        "kind": runtime.kind,
        "version": _runtime_version(runtime),
        "source": runtime.source,
        "path": _recorded_runtime_path(runtime),
        "program": spec.search,
    }
    if runtime.kind == "losat":
        record["cli"] = detect_losat_cli_dialect(runtime.executable)
    return record


# Argv and execution.


def build_losat_command(
    runtime: LosatRuntime,
    program: str | LosatProgram,
    *,
    query_path: Path,
    subject_path: Path,
    options: LosatSearchOptions = LosatSearchOptions(),
    threads: int | None = None,
    dialect: LosatCliDialect | None = None,
) -> list[str]:
    """Build the outfmt 6 argv for one search.

    LOSAT CLI dialects differ only in flags whose spellings differ; the dialect
    is probed only when such a flag is present.
    """

    spec = losat_program(program)
    values = _option_values(spec, options)
    if runtime.kind == "ncbi-blast":
        command = [str(runtime.executable)]
        spelling = "ncbi"
    else:
        command = [str(runtime.executable), spec.search]
        if dialect is None and any(
            LOSAT_OPTION_FLAGS[name].losat_v1 != LOSAT_OPTION_FLAGS[name].losat_v2
            for name, _value in values
        ):
            dialect = detect_losat_cli_dialect(runtime.executable)
        if dialect not in {None, "v1", "v2"}:
            raise ValidationError(f"Unknown LOSAT CLI dialect: {dialect!r}")
        spelling = "losat_v1" if dialect == "v1" else "losat_v2"
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
        command.extend([getattr(LOSAT_OPTION_FLAGS[name], spelling), value])
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
        raise ValidationError(f"{runtime_label} executable not found: {command[0]}") from exc
    except PermissionError as exc:
        raise ValidationError(f"{runtime_label} executable is not executable: {command[0]}") from exc
    except OSError as exc:
        raise ValidationError(
            f"{runtime_label} executable could not be started at {command[0]}: {exc}"
        ) from exc


def run_losat_search(
    program: str | LosatProgram,
    query_fasta: str,
    subject_fasta: str,
    *,
    options: LosatSearchOptions = LosatSearchOptions(),
    losat_bin: str | None = None,
    ncbi_blast_bin: str | None = None,
    threads: int | None = None,
    runtime_callback: LosatRuntimeCallback | None = None,
) -> str:
    """Run one search and return its raw outfmt 6 text.

    ``runtime_callback`` receives the runtime record after a successful search.
    """

    spec = losat_program(program)
    runtime_record: dict[str, object] | None = None
    with ExitStack() as stack:
        runtime = resolve_losat_runtime(
            spec,
            losat_bin=losat_bin,
            ncbi_blast_bin=ncbi_blast_bin,
            stack=stack,
        )
        temp_dir = stack.enter_context(tempfile.TemporaryDirectory(prefix="gbdraw_losat_"))
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
        logger.info("INFO: Running %s for %s.", runtime_label, spec.activity)
        completed = _run_losat_subprocess(command, runtime_label=runtime_label)
        if completed.returncode == 0 and runtime_callback is not None:
            runtime_record = losat_runtime_record(runtime, spec)

    if completed.returncode != 0:
        stderr = completed.stderr.strip()
        detail = f": {stderr}" if stderr else ""
        raise ValidationError(
            f"{runtime_label} failed at {runtime.executable} "
            f"with exit code {completed.returncode}{detail}"
        )
    if runtime_callback is not None and runtime_record is not None:
        runtime_callback(runtime_record)
    return completed.stdout


# Raw cache store.


class LosatRawCache:
    """Raw LOSAT entries by raw key, kept in Session display order.

    Program identity owners derive keys and validate entries. This store keeps
    the entries, their display filenames, and the non-key ``runtime`` record.
    """

    def __init__(self) -> None:
        self._entries_by_key: dict[str, dict[str, object]] = {}
        self._display_order: list[str] = []
        self._display_info: dict[str, tuple[str, bool]] = {}

    @property
    def has_entries(self) -> bool:
        return bool(self._entries_by_key)

    def _cached_entry(self, key: str) -> dict[str, object] | None:
        return self._entries_by_key.get(key)

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
            self._display_info[key] = (str(entry.get("filename") or ""), True)

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
        self._entries_by_key[key] = entry
        if display:
            self._mark_display(key, filename)

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
    "LosatCliDialect",
    "LosatProgram",
    "LosatRawCache",
    "LosatRuntime",
    "LosatSearchOptions",
    "NCBI_BLAST_BIN_OPTION",
    "build_losat_command",
    "bundled_losat_runtime",
    "detect_losat_cli_dialect",
    "losat_cache_args",
    "losat_program",
    "losat_runtime_record",
    "resolve_losat_runtime",
    "run_losat_search",
]
